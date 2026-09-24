"""
engine.py
=========
Faithful, vectorized port of the C++ kMC engine (calculate_reactions.cpp,
create_lattice.cpp population step, count_sites.cpp, initialize.cpp).

Algorithmic equivalences (documented on purpose)
-----------------------------------------------
* First-Reaction-Method  ==  Direct/BKL.  The C++ drew an independent waiting
  time t = -ln(u)/rate for every candidate event and kept the global minimum.
  That is mathematically identical to: total rate W = sum_j rate_j ; advance
  the clock by dt = -ln(u)/W ; choose the firing event with probability
  rate_j/W.  We use the BKL form because it is O(candidates) instead of
  O(candidates) random draws and is far faster, while reproducing the same
  statistics.  When several identical candidates exist (e.g. one decomposition
  per qualifying neighbour x per rate-channel) we fold them into a multiplicity.

* Placement (LiPlating / packOC / packBA) used rejection sampling over uniform
  random sites.  We instead sample uniformly among the eligible sites, which is
  the same stationary distribution, just without the wasted draws.

* RNG: we use numpy's PCG64 (seeded).  This is NOT bit-compatible with GSL's
  `taus`, so a single trajectory will differ realization-by-realization from
  the C++ run; the physics, observables and ensemble statistics are preserved.

The two interleaved loops (time-scale separation) are kept exactly: a SEI/
reaction call is followed by MaxToStop = randint(1,15) diffusion-only calls,
and diffusion events do not advance the kMC step counter.
"""

from __future__ import annotations
import numpy as np

from . import lattice as L
from .lattice import BCB, BCO, BA, OC, TE

# physical constants (match declarations.cpp)
JOULES_TO_EV = 6.24150974e18
KB_EV = 1.380649e-23 * JOULES_TO_EV           # Boltzmann in eV/K
AVOGADRO = 6.02214076e23
MOLAR_TO_ATOMS_PER_A3 = AVOGADRO / 1.0e27


class SpeciesTable:
    """Bidirectional species <-> integer-code map, code 0 == empty."""
    def __init__(self):
        self.sym2code = {"": 0}
        self.code2sym = {0: ""}
        for s in ("Li", "FSI", "SOL", "ETH", "F", "O", "S", "N", "SFO", "F5D"):
            self.add(s)

    def add(self, sym: str) -> int:
        if sym not in self.sym2code:
            c = len(self.sym2code)
            self.sym2code[sym] = c
            self.code2sym[c] = sym
        return self.sym2code[sym]

    def __getitem__(self, sym: str) -> int:
        return self.sym2code[sym]


class Engine:
    def __init__(self, lat: L.Lattice, params, geom, decomp, mob, mech, output):
        self.lat = lat
        self.p = params
        self.g = geom
        self.dec = decomp
        self.mob = mob
        self.mech = mech
        self.out = output
        self.rng = np.random.default_rng(params.seed)

        self.N = lat.n
        self.spt = SpeciesTable()
        # register any extra species appearing in mechanism / cathode.
        # Reservoir (dissolved) species are NOT lattice species: they live in
        # self.n_dis and are never given a site code.
        _RES_OPS = ("RES", "DISSOLVE")
        self.res_species = []
        for r in mech.reactions.values():
            if r.trigger and r.trigger != "EMPTY":
                self.spt.add(r.trigger)
            for (spc, _n) in r.requires:
                if spc not in self.res_species:
                    self.res_species.append(spc)
            if r.rate_scale and r.rate_scale not in self.res_species:
                self.res_species.append(r.rate_scale)
            for ch in r.channels.values():
                for (_op, spc, _q, _n) in ch:
                    if _op in _RES_OPS:
                        if spc not in self.res_species:
                            self.res_species.append(spc)
                    elif spc:
                        self.spt.add(spc)
        if geom.cathode_enabled:
            self.spt.add(geom.cathode_species)
        if geom.anode_species and geom.anode_species != "Li":
            self.spt.add(geom.anode_species)
        if geom.anode_vacancy_species:
            self.spt.add(geom.anode_vacancy_species)

        # state arrays
        self.occ = np.zeros(self.N, dtype=np.int8)
        self.spc = np.zeros(self.N, dtype=np.int32)      # species code
        self.chg = np.zeros(self.N, dtype=np.int32)

        # species code shortcuts
        S = self.spt
        self.cLi, self.cFSI, self.cSOL, self.cETH = S["Li"], S["FSI"], S["SOL"], S["ETH"]
        self.cF, self.cO = S["F"], S["O"]

        self._build_edges()
        self._build_interaction_matrix()
        self._build_minR()

        # scalars
        self.kT = KB_EV * params.temperature
        self.presentV = params.BeginV
        if str(getattr(params, "first_half", "begin")).lower() == "end":
            self.presentV = params.EndV          # start with a discharge
        self.idle_limit = int(getattr(params, "end_half_when_idle", 0))
        self.idle_events = 0                     # consecutive non-Li-transfer events
        self.cycleNumber = 0
        self.currentStep = 0
        self.currentTime = 0.0
        self.timeCurrentV = 0.0
        self.stepCurrentV = 0
        self.stop = 0
        self.MaxToStop = 0
        self.failedAttempts = 0
        self.electrolyteVol = 0.0

        # reaction counters
        self.R = dict(Plating=0, PlatingSEI=0, FSI=0, SOL=0, SFO=0, F5D=0,
                      SOL2=0, LiStripping=0, LiSurface=0)
        # cathode / conversion reactions get their own counters
        self.conv_rxns = [r for r in mech.reactions.values() if r.is_conversion]
        for r in self.conv_rxns:
            self.R.setdefault(r.name, 0)
        # CATHODE-AWARE PATCH: cathode material must never be classified as
        # anode SEI, otherwise wrap()/updateLiMetal() plate Li around the
        # cathode slab (validated artifact). Collect every cathode species:
        # the geometry species + conversion triggers + CONVERT products.
        _cath = set()
        if getattr(geom, "cathode_enabled", False) and geom.cathode_species:
            _cath.add(self.spt.add(geom.cathode_species))
        for r in self.conv_rxns:
            # only CATHODE-region reactions define cathode material; the
            # anode-surface shuttle triggers on Li and must not put Li here.
            if r.region != "cathode":
                continue
            if r.trigger != "EMPTY":
                _cath.add(self.spt.add(r.trigger))
            for ch in r.channels.values():
                for (op, sname, _q, _n) in ch:
                    if op in ("CONVERT", "PRECIP"):
                        _cath.add(self.spt.add(sname))
        self.cath_code_arr = np.array(sorted(_cath), dtype=self.spc.dtype)
        # DEPOSIT products (e.g. Li2S2_an): inert passivating deposits on the
        # anode. They are excluded from the SEI class so the anode framework
        # (wrap / updateLiMetal) neither grows Li around them nor removes Li
        # because of them. Empty for any mechanism without DEPOSIT -> legacy
        # masks are untouched.
        _dep = set()
        _cei = set()
        for r in mech.reactions.values():
            for ch in r.channels.values():
                for (op, sname, _q, _n) in ch:
                    if op == "DEPOSIT":
                        _dep.add(self.spt.add(sname))
                    elif op == "CEI":
                        # CEI film species: same inert treatment as deposits
                        _dep.add(self.spt.add(sname))
                        _cei.add(self.spt.add(sname))
        self.dep_code_arr = np.array(sorted(_dep), dtype=self.spc.dtype)
        self.cei_code_arr = np.array(sorted(_cei), dtype=self.spc.dtype)
        self.cei_active = bool(_cei)
        self.s_cei = 0            # S atoms sequestered in the CEI (S_LOSS)
        self.li_cei = 0           # Li trapped in the CEI (LI_LOSS, reporting)
        self.deposit_unplaced = 0 # DEPOSIT sites that found no room (audit)
        # ---- structural brakes (opt-in) ----
        self.li_pool_max = int(getattr(params, "li_pool_max", -1))
        if str(getattr(params, "li_pool_mode", "fixed")).lower() != "shared":
            self.li_pool_max = -1          # meaningless without a pool
        self.stop_elec_frac = float(getattr(params, "stop_electrolyte_fraction", 0.0))
        self.n_elec0 = 0                   # electrolyte sites at t=0
        self.stop_anode_halves = int(getattr(params, "stop_anode_inactive_halves", 0))
        # ---- physically based anode kinetics (opt-in) ----
        self.bv = str(getattr(params, "anode_kinetics", "legacy")).lower() == "bv"
        self.bv_k0 = float(getattr(params, "bv_k0_site", 74.0))
        self.bv_alpha = float(getattr(params, "bv_alpha", 0.5))
        self.eta_c = float(getattr(params, "eta_charge", -0.03))
        self.eta_d = float(getattr(params, "eta_discharge", 0.03))
        self.kT_V = KB_EV * params.temperature      # k_B T / e in volts
        self.cath_bv = str(getattr(params, "cathode_kinetics", "legacy")).lower() == "bv"
        self.cath_V_c = float(getattr(params, "cathode_V_charge", 2.45))
        self.cath_V_d = float(getattr(params, "cathode_V_discharge", 1.9))
        self.anode_dead = False
        self.dried_out = False
        self.stalled = False      # set by step() when no event can fire
        self._last_W = 0.0        # total non-diffusion rate of the last scan
        self.pool_full_blocks = 0          # audit: candidate sets pruned by the cap
        self.li_consumed = 0      # Li+ consumed by cathode reduction (bookkeeping)
        self.li_released = 0      # Li+ released by cathode oxidation (bookkeeping)
        # ---- shared Li+ pool (li_pool_mode shared) + well-mixed reservoir ----
        self.pool_shared = (str(getattr(params, "li_pool_mode", "fixed")).lower()
                            == "shared")
        self.li_pool = 0          # Li+ available in the electrolyte (shared mode)
        self.li_pool0 = 0         # initial pool (normalises the plating rate)
        self.li_shuttled = 0      # Li0 corroded from the anode by the shuttle
        self.li_deposit = 0       # insoluble Li2Sx sites deposited on the anode
        self.li_plated_pool = 0   # Li0 placed through the pool (shared mode)
        # anode bulk-metal reservoir (shared mode): the framework draws from it
        self.li_bulk_init = int(getattr(params, "li_bulk_init", -1))
        self.li_bulk = max(self.li_bulk_init, 0)
        self.li_bulk_drawn = 0    # Li created by wrap()/fill (from the foil)
        self.li_bulk_returned = 0 # floating Li removed by updateLiMetal (to foil)
        self.n_dis = {spc: 0 for spc in self.res_species}   # dissolved counts
        self.reservoir_sites = int(getattr(params, "reservoir_sites", -1))
        self.cath_access = None   # static accessibility mask (None = all)
        # Cathode lattice sites (BC sites inside the cathode z-band) are NOT
        # anode framework: updateLiMetal / wrap must never fill them with Li
        # and createEther must never fill them with ETH. Without this, a
        # dissolved S8 site was refilled with foil Li on the next step (Li
        # metal inside the cathode, re-precipitation impossible: found
        # 2026-09-21). None when the cathode is disabled -> legacy path.
        self.cath_bc_band = None
        if getattr(geom, "cathode_enabled", False):
            _z = lat.frac[:, 2]
            self.cath_bc_band = lat.is_BC() & self._zband(
                _z, geom.cathode_center, geom.cathode_thickness / 2.0)
        # ---- dynamic Li2S passivation of the cathode (opt-in, point 3) ----
        # Active only when cathode_passivation_nmin > 0 AND the species list is
        # non-empty. Species codes are looked up lazily in _passivation_mask()
        # so the species table is NOT touched when the feature is off.
        self.pass_nmin = int(getattr(geom, "cathode_passivation_nmin", 0))
        _pl = str(getattr(geom, "cathode_passivating_species", "") or "")
        self.pass_species = [t for t in _pl.replace(";", ",").split(",") if t.strip()]
        self.pass_species = [t.strip() for t in self.pass_species]
        if self.pass_nmin > 0 and not self.pass_species:
            self.pass_nmin = 0     # nothing can block -> feature off
        self.pass_code_arr = None  # filled on first use (codes are stable)
        self.n_cathode_blocked = 0 # last evaluated number of blocked sites
        # per-reaction flag: does any channel transfer Li+ (electron transfer)?
        # Only these are gated by the passivation mask.
        self._li_transfer = {}
        for r in self.conv_rxns:
            self._li_transfer[r.name] = any(
                op in ("CONSUME_LI", "RELEASE_LI")
                for ch in r.channels.values() for (op, _s, _q, _n) in ch)
        # per-reaction Li+ requirement (max CONSUME_LI over channels) for gating
        self._li_need = {}
        # per-reaction deposit requirement (max DEPOSIT count over channels):
        # a shuttle reaction is a candidate only if that many eligible BA
        # surface sites exist, otherwise the deposit could not be placed and
        # its sulfur would vanish (found with s_total, 2026-09-20). Physically:
        # a Li surface saturated with Li2S2 is passivated toward the shuttle.
        self._dep_need = {}
        self._li_release = {}
        for r in self.conv_rxns:
            need = 0
            for ch in r.channels.values():
                need = max(need, sum(n for (op, _s, _q, n) in ch if op == "CONSUME_LI"))
            self._li_need[r.name] = need
            dneed = 0
            for ch in r.channels.values():
                dneed = max(dneed, sum(n for (op, _s, _q, n) in ch if op == "DEPOSIT"))
            self._dep_need[r.name] = dneed
            self._li_release[r.name] = max(
                (sum(n for (op, _s, _q, n) in ch if op == "RELEASE_LI")
                 for ch in r.channels.values()), default=0)
        self.LiMetal = self.LiIon = self.LiMetalS = self.LiIonS = 0
        self.nO = self.nF = 0

        # per-half-cycle ML ledger (cycle_stats.csv)
        from .stats import CycleLedger
        self.ledger = CycleLedger()
        if geom.anode_thickness_is_fraction:
            _h = geom.anode_thickness / 2.0
        else:
            _h = (geom.anode_thickness_A / self.lat.box[2, 2]) / 2.0
        _Lz = self.lat.box[2, 2]
        self._anode_band = ((geom.anode_center - _h) * _Lz,
                            (geom.anode_center + _h) * _Lz)
        # last fired reaction (for logging)
        self.last = dict(type="NoRxn", Ospcs="", Ea=0.0, expValue=0.0,
                         reactRate=0.0, time=0.0)

    # ------------------------------------------------------------------ setup
    def _build_edges(self):
        """Combined and grouped directed edge arrays + per-edge distances."""
        lat = self.lat
        cart = lat.frac @ lat.box
        Lv = np.array([lat.box[0, 0], lat.box[1, 1], lat.box[2, 2]])

        src_all, dst_all = [], []
        rel_of_edge = []
        for r in (BCB, BCO, BA, OC, TE):
            ptr, idx = lat.nbr_indptr[r], lat.nbr_indices[r]
            counts = np.diff(ptr)
            s = np.repeat(np.arange(self.N), counts)
            src_all.append(s); dst_all.append(idx)
            rel_of_edge.append(np.full(idx.size, r, dtype=np.int8))
        self.src = np.concatenate(src_all)
        self.dst = np.concatenate(dst_all)
        rel = np.concatenate(rel_of_edge)

        d = cart[self.src] - cart[self.dst]
        d -= Lv * np.round(d / Lv)
        self.edist = np.sqrt((d * d).sum(axis=1))

        # grouped edge masks for vectorized neighbour counts
        self.m_g4 = np.isin(rel, (BCB, BA, BCO, OC))   # for nLi/nSOL/nETH/nSALT
        self.m_te = rel == TE                           # nLi also counts TE
        self.m_ba2 = np.isin(rel, (BCB, BA))            # nBA_SEI
        self.m_oc2 = np.isin(rel, (BCO, OC))            # nOC_SEI
        self.rel = rel

        # PERFORMANCE: the group masks are static, so pre-filter the edge lists
        # once. _count then avoids re-masking the full edge arrays every call
        # (this was ~80% of total runtime).
        self.eg4 = (self.src[self.m_g4].astype(np.int64),
                    self.dst[self.m_g4].astype(np.int64))
        self.ete = (self.src[self.m_te].astype(np.int64),
                    self.dst[self.m_te].astype(np.int64))
        self.eba2 = (self.src[self.m_ba2].astype(np.int64),
                     self.dst[self.m_ba2].astype(np.int64))
        self.eoc2 = (self.src[self.m_oc2].astype(np.int64),
                     self.dst[self.m_oc2].astype(np.int64))

    def _build_interaction_matrix(self):
        n = len(self.spt.sym2code)
        M = np.zeros((n, n))
        for a, row in self.mob.interaction.items():
            if a not in self.spt.sym2code:
                continue
            ia = self.spt.sym2code[a]
            for b, eb in row.items():
                if b in self.spt.sym2code:
                    M[ia, self.spt.sym2code[b]] = eb
        self.INT = M

    def _build_minR(self):
        """Per-site coordination radius used by Eact / charge assignment.
        For OC-type SEI sites: distance to nearest BA neighbour.
        For BA-type SEI sites: distance to nearest TE neighbour.
        Non-SEI sites get a large value (no restriction)."""
        lat = self.lat
        minR = np.full(self.N, 1000.0)
        # nearest BA neighbour distance per site
        for r, mask_site in ((BA, lat.is_OCsite()), (TE, lat.is_BAsite())):
            ptr, idx = lat.nbr_indptr[r], lat.nbr_indices[r]
            for i in np.where(mask_site)[0]:
                a, b = ptr[i], ptr[i + 1]
                if b > a:
                    # distance of first such neighbour (regular lattice => uniform)
                    e0 = self._edge_dist(i, idx[a])
                    if r == BA:                 # OC-type uses BA distance
                        minR[i] = min(minR[i], e0) if minR[i] < 1000 else e0
                    else:                       # BA-type uses TE distance
                        minR[i] = min(minR[i], e0) if minR[i] < 1000 else e0
        self.minR_geom = minR   # geometric; SEI gating applied at use time

    def _edge_dist(self, i, j):
        cart = self.lat.frac @ self.lat.box
        Lv = np.array([self.lat.box[0, 0], self.lat.box[1, 1], self.lat.box[2, 2]])
        d = cart[i] - cart[j]
        d -= Lv * np.round(d / Lv)
        return float(np.sqrt((d * d).sum()))

    # --------------------------------------------------------------- masks
    def _sei_mask(self):
        m = (self.occ == 1) & (self.spc != self.cLi) & (self.spc != self.cSOL) \
            & (self.spc != self.cFSI) & (self.spc != self.cETH)
        if self.cath_code_arr.size:
            m &= ~np.isin(self.spc, self.cath_code_arr)
        if self.dep_code_arr.size:
            m &= ~np.isin(self.spc, self.dep_code_arr)
        return m

    def _nonsei_mask(self):
        """Sites the Li framework may overwrite (wrap): empty, Li or
        electrolyte, and never a cathode lattice site."""
        s = self.spc
        m = (s == 0) | (s == self.cLi) | (s == self.cSOL) | (s == self.cFSI) \
            | (s == self.cETH)
        if self.cath_bc_band is not None:
            m &= ~self.cath_bc_band
        return m

    # site-class codes for the combined neighbour-count pass
    _CLS_NONE, _CLS_LI, _CLS_ETH, _CLS_SOL, _CLS_FSI, _CLS_SEI, _CLS_CATH = range(7)
    _NCLS = 7

    def _count(self, edges, dst_value_mask):
        """Count, per src site, edges (src_pre, dst_pre) whose dst satisfies
        dst_value_mask. Pure-integer bincount on pre-filtered edge lists."""
        src_pre, dst_pre = edges
        return np.bincount(src_pre[dst_value_mask[dst_pre]], minlength=self.N)

    def _class_of_sites(self):
        """Per-site class label: 0 empty, 1 Li, 2 ETH, 3 SOL, 4 FSI, 5 SEI."""
        cls = np.zeros(self.N, dtype=np.int64)
        occ1 = self.occ == 1
        s = self.spc
        cls[occ1 & (s == self.cLi)] = self._CLS_LI
        cls[occ1 & (s == self.cETH)] = self._CLS_ETH
        cls[occ1 & (s == self.cSOL)] = self._CLS_SOL
        cls[occ1 & (s == self.cFSI)] = self._CLS_FSI
        cls[occ1 & (s != self.cLi) & (s != self.cETH)
            & (s != self.cSOL) & (s != self.cFSI)] = self._CLS_SEI
        if self.cath_code_arr.size:
            cls[occ1 & np.isin(s, self.cath_code_arr)] = self._CLS_CATH
        if self.dep_code_arr.size:
            cls[occ1 & np.isin(s, self.dep_code_arr)] = self._CLS_CATH
        return cls

    def _count_all_classes(self, edges, cls):
        """One bincount that yields, per src site, the neighbour count of every
        class simultaneously: returns (N, NCLS) int array."""
        src_pre, dst_pre = edges
        key = src_pre * self._NCLS + cls[dst_pre]
        return np.bincount(key, minlength=self.N * self._NCLS)\
                 .reshape(self.N, self._NCLS)

    def _refresh_counts(self):
        cls = self._class_of_sites()
        g4 = self._count_all_classes(self.eg4, cls)
        te = self._count_all_classes(self.ete, cls)
        self.nLi = g4[:, self._CLS_LI] + te[:, self._CLS_LI]
        self.nETH = g4[:, self._CLS_ETH]
        self.nSOL = g4[:, self._CLS_SOL]
        self.nSALT = g4[:, self._CLS_FSI]
        seiM = cls == self._CLS_SEI
        self.nBA_SEI = self._count(self.eba2, seiM)
        self.nOC_SEI = self._count(self.eoc2, seiM)
        self.seiM = seiM

    # --------------------------------------------------------- Eact / rates
    def _eact_all(self):
        """Environment-dependent activation energy for every site (vectorized)."""
        occ1 = self.occ == 1
        is_li_dst = self.spc[self.dst] == self.cLi
        # SEI sites: restrict to minR; non-SEI: no restriction (minR huge)
        minR = np.where(self.seiM, self.minR_geom, 1000.0)
        within = np.where(is_li_dst, self.edist <= minR[self.src], True)
        Iv = self.INT[self.spc[self.src], self.spc[self.dst]]
        valid = occ1[self.dst] & within & (Iv > 0) & (Iv < 1)
        return np.bincount(self.src, weights=Iv * valid, minlength=self.N)

    def _arr_rate(self, sigma, k0, Ea, alpha, E0, V=None):
        """rate = globalA sigma k0 exp(-(Ea - alpha (V - E0)) / kT). V defaults
        to the cell voltage label (legacy); the cathode BV mode passes the
        cathode potential of the half-cycle instead."""
        Vx = self.presentV if V is None else V
        expv = (Ea - alpha * (Vx - E0)) / self.kT
        return self.p.globalA * sigma * k0 * np.exp(-expv), expv

    def _decomp_sumrate(self, rtype, V=None):
        """Sum of Arrhenius rates over all (non-zero sigma) rate channels."""
        tot = 0.0
        for c in self.dec.channels.get(rtype, []):
            if c["sigma"] != 0.0:
                r, _ = self._arr_rate(c["sigma"], c["k0"], c["Ea"],
                                      c["alpha"], c["E0"], V)
                tot += r
        return tot

    def _cathode_V(self):
        """Cathode potential (V vs Li/Li+) of the current half-cycle."""
        return self.cath_V_c if self.presentV == self.p.BeginV else self.cath_V_d

    # --------------------------------------------------------- population ops
    def _eligible_BC(self):
        return (self.occ != 1) & self.lat.is_BC()

    def LiPlating(self, n):
        placed = 0
        for _ in range(n):
            self._refresh_counts()
            elig = self._eligible_BC() & (self.nETH > 0) & (self.nLi >= 2) & (self.nLi <= 5)
            idx = np.where(elig)[0]
            if idx.size == 0:
                break
            i = idx[self.rng.integers(idx.size)]
            self.occ[i] = 1; self.spc[i] = self.cLi; self.chg[i] = 0
            placed += 1
        if self.pool_shared:
            self.li_pool -= placed
            self.li_plated_pool += placed
        return placed

    def packOC(self, sym, charge, n):
        code = self.spt.add(sym)
        for _ in range(n):
            self._refresh_counts()
            elig = (self.occ != 1) & self.lat.is_OCsite() & (self.nETH >= 1) \
                & (self.nBA_SEI == 0) & (self.nLi >= 2)
            idx = np.where(elig)[0]
            if idx.size == 0:
                return
            i = idx[self.rng.integers(idx.size)]
            self.occ[i] = 1; self.spc[i] = code; self.chg[i] = charge

    def packBA(self, sym, charge, n):
        code = self.spt.add(sym)
        for _ in range(n):
            self._refresh_counts()
            elig = (self.occ != 1) & self.lat.is_BAsite() & (self.nETH >= 1) \
                & (self.nOC_SEI == 0) & (self.nLi >= 2)
            idx = np.where(elig)[0]
            if idx.size == 0:
                return
            i = idx[self.rng.integers(idx.size)]
            self.occ[i] = 1; self.spc[i] = code; self.chg[i] = charge

    def _apply_channel(self, actions, site=None):
        for (op, spc, q, count) in actions:
            if op == "PLATE_LI":
                self.LiPlating(1)
            elif op == "PACK_OC":
                self.packOC(spc, q, count)
            elif op == "PACK_BA":
                self.packBA(spc, q, count)
            elif op == "CONVERT" and site is not None:
                code = self.spt.add(spc)
                self.occ[site] = 1; self.spc[site] = code; self.chg[site] = q
            elif op == "CONSUME_LI":
                self.li_consumed += count
                if self.pool_shared:
                    self.li_pool -= count
            elif op == "RELEASE_LI":
                self.li_released += count
                if self.pool_shared:
                    self.li_pool += count
            # ---- well-mixed reservoir / shuttle actions -------------------
            elif op == "RES":
                self.n_dis[spc] = self.n_dis.get(spc, 0) + count
            elif op == "DISSOLVE" and site is not None:
                self._empty_one(site)
                self.n_dis[spc] = self.n_dis.get(spc, 0) + 1
            elif op == "PRECIP" and site is not None:
                code = self.spt.add(spc)
                self.occ[site] = 1; self.spc[site] = code; self.chg[site] = q
            elif op == "CEI" and site is not None:
                code = self.spt.add(spc)
                self.occ[site] = 1; self.spc[site] = code; self.chg[site] = q
            elif op == "S_LOSS":
                self.s_cei += count
            elif op == "LI_LOSS":
                self.li_cei += count
            elif op == "STRIP_LI0":
                self._strip_li0(site, count)
            elif op == "DEPOSIT":
                placed = self.packBA_count(spc, q, count)
                self.li_deposit += placed
                if placed < count:
                    # should not happen after the candidate gate; keep the
                    # balance auditable rather than silent
                    self.deposit_unplaced += count - placed
                    self.out.log(f"WARNING: DEPOSIT {spc} placed {placed}/{count} "
                                 f"at step {self.currentStep}; S balance affected")

    def _surface_li0_mask(self):
        """Li0 sites exposed to electrolyte (same criterion as LiStripping)."""
        occ1 = self.occ == 1
        li0 = occ1 & (self.spc == self.cLi) & (self.chg == 0)
        return li0 & (self.nLi <= 4) & (self.nETH >= 1)

    def _strip_li0(self, site, n):
        """Anode corrosion by the shuttle: remove the trigger Li0 site and up to
        n-1 further surface Li0 sites. The Li ends up bound in polysulfide
        (reservoir stoichiometry carries it), so the Li+ pool is untouched."""
        removed = 0
        if site is not None and self.occ[site] == 1 and self.spc[site] == self.cLi:
            self._empty_one(site); removed += 1
        while removed < n:
            self._refresh_counts()
            idx = np.where(self._surface_li0_mask())[0]
            if idx.size == 0:
                break
            i = idx[self.rng.integers(idx.size)]
            self._empty_one(i); removed += 1
        self.li_shuttled += removed
        return removed

    def _pool_full(self):
        """True when the shared Li+ pool is at its cap (li_pool_max)."""
        return self.li_pool_max > 0 and self.li_pool >= self.li_pool_max

    def _electrolyte_count(self):
        occ1 = self.occ == 1
        return int((occ1 & ((self.spc == self.cETH) | (self.spc == self.cSOL)
                            | (self.spc == self.cFSI))).sum())

    def _deposit_eligible_count(self):
        """Number of BA sites where packBA_count could place a deposit now
        (same criterion as packBA_count; counts are refreshed by sei_step)."""
        elig = (self.occ != 1) & self.lat.is_BAsite() & (self.nETH >= 1) \
            & (self.nOC_SEI == 0) & (self.nLi >= 2)
        return int(elig.sum())

    def packBA_count(self, sym, charge, n):
        """packBA that reports how many sites were actually placed."""
        code = self.spt.add(sym)
        placed = 0
        for _ in range(n):
            self._refresh_counts()
            elig = (self.occ != 1) & self.lat.is_BAsite() & (self.nETH >= 1) \
                & (self.nOC_SEI == 0) & (self.nLi >= 2)
            idx = np.where(elig)[0]
            if idx.size == 0:
                break
            i = idx[self.rng.integers(idx.size)]
            self.occ[i] = 1; self.spc[i] = code; self.chg[i] = charge
            placed += 1
        return placed

    def _apply_conversion(self, reaction, site):
        """Apply a cathode conversion reaction in place at `site`."""
        self._apply_channel(reaction.pick_channel(self.rng), site=site)

    def electrolyteDecom(self, rtype):
        if rtype in self.mech:
            self._apply_channel(self.mech.get(rtype).pick_channel(self.rng))
        # counters
        if rtype in self.R:
            self.R[rtype] += 1

    def flushElectrolyte(self):
        m = (self.occ == 1) & ((self.spc == self.cSOL) | (self.spc == self.cFSI))
        self._empty(m)

    def _empty(self, mask):
        self.occ[mask] = 0; self.spc[mask] = 0; self.chg[mask] = 0

    def createEther(self):
        self._refresh_counts()
        m = (self.occ != 1) & (self.lat.is_OCsite() | self.lat.is_BAsite()) \
            & (self.nLi <= 2)
        if self.cath_bc_band is not None:
            m &= ~self.cath_bc_band       # cathode BCO sites stay EMPTY
        self.occ[m] = 1; self.spc[m] = self.cETH; self.chg[m] = 0

    def updateEther(self):
        m = (self.occ == 1) & ((self.spc == self.cSOL) | (self.spc == self.cFSI)
                               | (self.spc == self.cETH))
        self._empty(m)
        self.createEther()

    def etherVol(self):
        rTE = 1.0
        nETH_sites = int(((self.occ == 1) & (self.spc == self.cETH)).sum())
        vol = nETH_sites * 4.0 * np.pi * rTE**3 / 3.0
        self.electrolyteVol = vol / 0.74

    def addSolvent(self):
        self.etherVol()
        to_add = int(round(self.p.solvMolar * self.electrolyteVol * MOLAR_TO_ATOMS_PER_A3))
        eth = np.where((self.occ == 1) & (self.spc == self.cETH))[0]
        if to_add > 0 and eth.size:
            pick = self.rng.choice(eth, size=min(to_add, eth.size), replace=False)
            self.spc[pick] = self.cSOL

    def addLithiumSalt(self):
        self.etherVol()
        to_add = int(round(self.p.saltMolar * self.electrolyteVol * MOLAR_TO_ATOMS_PER_A3))
        eth = np.where((self.occ == 1) & (self.spc == self.cETH))[0]
        if to_add > 0 and eth.size:
            pick = self.rng.choice(eth, size=min(to_add, eth.size), replace=False)
            self.spc[pick] = self.cFSI

    # ---- anode bulk reservoir helpers (shared mode) ----------------------
    def _bulk_take(self, n):
        """Framework wants to create n Li sites: draw from the foil.
        Unlimited foil (li_bulk_init < 0): always granted, li_bulk goes
        negative (= net Li drawn). Finite foil: capped at what is left."""
        if n <= 0:
            return 0
        if self.li_bulk_init >= 0:
            n = min(n, self.li_bulk)
        self.li_bulk -= n
        self.li_bulk_drawn += n
        return n

    def _bulk_return(self, n):
        if n > 0:
            self.li_bulk += n
            self.li_bulk_returned += n

    def wrap(self):
        if self.pool_shared:
            return self._wrap_accounted()
        return self._wrap_legacy()

    def _wrap_accounted(self):
        """Same geometry as _wrap_legacy, but every Li site the framework
        creates is drawn from the anode bulk reservoir (and capped by it)."""
        seiM = self._sei_mask()
        lat = self.lat

        def _fill(idx):
            f = idx[self._nonsei_mask()[idx]]
            new = f[~((self.occ[f] == 1) & (self.spc[f] == self.cLi))]
            k = self._bulk_take(new.size)
            new = new[:k]
            self.occ[new] = 1; self.spc[new] = self.cLi
            return f.size

        oc_sei = np.where(seiM & lat.is_OCsite())[0]
        for i in oc_sei:
            ba = lat.nbr_indices[BA][lat.nbr_indptr[BA][i]:lat.nbr_indptr[BA][i + 1]]
            if _fill(ba) < 6:
                bcb = lat.nbr_indices[BCB][lat.nbr_indptr[BCB][i]:lat.nbr_indptr[BCB][i + 1]]
                _fill(bcb)
        ba_sei = np.where(seiM & lat.is_BAsite())[0]
        for i in ba_sei:
            te = lat.nbr_indices[TE][lat.nbr_indptr[TE][i]:lat.nbr_indptr[TE][i + 1]]
            _fill(te)

    def _wrap_legacy(self):
        seiM = self._sei_mask()
        lat = self.lat
        nonsei = self._nonsei_mask()
        # SEI OC-sites -> fill nonSEI BA neighbours with Li; if <6 also BCB
        oc_sei = np.where(seiM & lat.is_OCsite())[0]
        for i in oc_sei:
            ba = lat.nbr_indices[BA][lat.nbr_indptr[BA][i]:lat.nbr_indptr[BA][i + 1]]
            fill = ba[self._nonsei_mask()[ba]]
            self.occ[fill] = 1; self.spc[fill] = self.cLi
            if fill.size < 6:
                bcb = lat.nbr_indices[BCB][lat.nbr_indptr[BCB][i]:lat.nbr_indptr[BCB][i + 1]]
                f2 = bcb[self._nonsei_mask()[bcb]]
                self.occ[f2] = 1; self.spc[f2] = self.cLi
        ba_sei = np.where(seiM & lat.is_BAsite())[0]
        for i in ba_sei:
            te = lat.nbr_indices[TE][lat.nbr_indptr[TE][i]:lat.nbr_indptr[TE][i + 1]]
            f = te[self._nonsei_mask()[te]]
            self.occ[f] = 1; self.spc[f] = self.cLi

    def updateLiMetal(self):
        lat = self.lat
        self._refresh_counts()
        liM = (self.occ == 1) & (self.spc == self.cLi)
        is_ba = lat.name_code == BA
        is_te = lat.name_code == TE
        is_bco = lat.name_code == BCO
        # remove floating / over-coordinated Li on BA/TE
        rm = liM & (is_ba | is_te) & (
            ((self.nBA_SEI + self.nOC_SEI) == 0) |
            ((self.nOC_SEI > 0) & (self.nBA_SEI > 0)))
        if self.pool_shared:
            self._bulk_return(int(rm.sum()))
        self._empty(rm)
        # remove BCO Li adjacent to BA-SEI
        rm2 = liM & is_bco & (self.nBA_SEI > 0)
        if self.pool_shared:
            self._bulk_return(int(rm2.sum()))
        self._empty(rm2)
        # fill vacuum below SEI
        self._refresh_counts()
        fill = (self.occ != 1) & lat.is_BC() & \
            ((self.nOC_SEI + self.nBA_SEI) == 0) & (self.nETH == 0)
        if self.cath_bc_band is not None:
            fill &= ~self.cath_bc_band    # no foil Li inside the cathode
        if self.pool_shared:
            idx = np.where(fill)[0]
            idx = idx[:self._bulk_take(idx.size)]
            self.occ[idx] = 1; self.spc[idx] = self.cLi; self.chg[idx] = 0
        else:
            self.occ[fill] = 1; self.spc[fill] = self.cLi; self.chg[fill] = 0
        # recount metrics
        liM = (self.occ == 1) & (self.spc == self.cLi)
        self.LiMetal = int((liM & (self.chg == 0)).sum())
        self.LiIon = int((liM & (self.chg != 0)).sum())
        self._refresh_counts()
        surf = liM & (self.nETH > 0)
        self.LiMetalS = int((surf & (self.chg == 0)).sum())
        self.LiIonS = int((surf & (self.chg != 0)).sum())

    def updateCharges(self):
        liM = (self.occ == 1) & (self.spc == self.cLi)
        self.chg[liM] = 0
        seiM = self._sei_mask()
        # edges SEI(src) -> Li(dst) within minR  => dst charge = 1
        sei_src = seiM[self.src]
        li_dst = (self.occ[self.dst] == 1) & (self.spc[self.dst] == self.cLi)
        minR = np.where(seiM, self.minR_geom, 1000.0)
        within = self.edist <= minR[self.src]
        hit_dst = self.dst[sei_src & li_dst & within]
        self.chg[np.unique(hit_dst)] = 1

    def countOandF(self):
        self.nO = int(((self.occ == 1) & (self.spc == self.cO)).sum())
        self.nF = int(((self.occ == 1) & (self.spc == self.cF)).sum())

    # ------------------------------------------------------------- electrodes
    def create_cell(self):
        self.createAnode()
        if self.g.cathode_enabled:
            self.createCathode()
        self.createEther()
        self.addSolvent()
        self.addLithiumSalt()

    def createAnode(self):
        z = self.lat.frac[:, 2]
        if self.g.anode_thickness_is_fraction:
            half = self.g.anode_thickness / 2.0
        else:
            half = (self.g.anode_thickness_A / self.lat.box[2, 2]) / 2.0
        c = self.g.anode_center
        m = self.lat.is_BC() & (z >= c - half) & (z <= c + half)
        Lz = self.lat.box[2, 2]
        self._anode_band = ((c - half) * Lz, (c + half) * Lz)
        base_code = self.spt.add(self.g.anode_species) if self.g.anode_species else self.cLi
        self.occ[m] = 1
        self.spc[m] = base_code
        self.chg[m] = 0
        # carve vacancies: knock out a random fraction of the anode BC sites
        vf = float(self.g.anode_vacancy_fraction)
        if vf > 0.0:
            idx = np.where(m)[0]
            n_hole = int(round(vf * idx.size))
            if n_hole > 0:
                holes = self.rng.choice(idx, size=n_hole, replace=False)
                if self.g.anode_vacancy_species:
                    vcode = self.spt.add(self.g.anode_vacancy_species)
                    self.spc[holes] = vcode          # filled hole (e.g. Zn)
                    self.chg[holes] = 0
                else:
                    self.occ[holes] = 0              # empty hole
                    self.spc[holes] = 0
                    self.chg[holes] = 0

    @staticmethod
    def _zband(z, c, half):
        """Sites with fractional z inside [c-half, c+half]. If the band crosses
        the periodic boundary (e.g. cathode_center 0.0 -> slab wrapping both
        box faces) use the minimum-image distance; otherwise keep the exact
        legacy inequality so validated cases stay bit-identical."""
        if c - half >= 0.0 and c + half <= 1.0:
            return (z >= c - half) & (z <= c + half)
        d = np.abs(z - c)
        d = np.minimum(d, 1.0 - d)
        return d <= half

    def createCathode(self):
        z = self.lat.frac[:, 2]
        half = self.g.cathode_thickness / 2.0
        c = self.g.cathode_center
        code = self.spt.add(self.g.cathode_species)
        m = self.lat.is_BC() & self._zband(z, c, half) & (self.occ != 1)
        self.occ[m] = 1; self.spc[m] = code; self.chg[m] = self.g.cathode_charge
        # static electronic/ionic accessibility mask (point 4, static version):
        # conversion only fires on accessible cathode sites.
        fa = float(getattr(self.g, "cathode_access_fraction", 1.0))
        if fa < 1.0:
            sites = np.where(m)[0]
            keep = self.rng.random(sites.size) < fa
            self.cath_access = np.zeros(self.N, dtype=bool)
            self.cath_access[sites[keep]] = True

    # ------------------------------------------------------------- the loops
    def _voltage_update(self):
        if self.p.potentialType == "step" and self.p.scanIntervalType == "time":
            if (self.timeCurrentV >= self.p.scanInterval
                    or self.stepCurrentV >= self.p.maxInterval
                    or (self.idle_limit > 0 and self.idle_events >= self.idle_limit)):
                if self.idle_limit > 0 and self.idle_events >= self.idle_limit:
                    self.out.log(f"half-cycle {self.cycleNumber} ended idle after "
                                 f"{self.idle_events} events without Li transfer "
                                 f"(step {self.currentStep})")
                self.idle_events = 0
                self.countOandF()
                self._write_xyz()
                self.ledger.close_half(self)
                self.presentV = self.p.EndV if self.presentV == self.p.BeginV else self.p.BeginV
                self.stepCurrentV = 0
                self.timeCurrentV = 0.0
                self.cycleNumber += 1
                self.currentStep += 1
                nsc = getattr(self.p, "snapshotEveryCycles", 0)
                if nsc and self.cycleNumber % nsc == 0:
                    self._write_snapshot()
        self.last["type"] = "NoRxn"

    def _draw_dt(self, W):
        u = self.rng.random()
        return -np.log(max(u, 1e-300)) / W

    def _diffusion_candidates(self):
        """Return (sites, rate_per_site, dest_pools) for SEI diffusion."""
        self._refresh_counts()
        seiM = self.seiM
        lat = self.lat
        eact = self._eact_all()
        # per-site diffusion rate (species-dependent k0/sigma; Ea=eact)
        sigma = np.zeros(self.N); k0 = np.zeros(self.N)
        alpha = np.zeros(self.N); E0 = np.zeros(self.N)
        for sym, pr in self.mob.params.get("Diffusion", {}).items():
            if sym in self.spt.sym2code:
                c = self.spt.sym2code[sym]
                sel = self.spc == c
                sigma[sel] = pr["sigma"]; k0[sel] = pr["k0"]
                alpha[sel] = pr["alpha"]; E0[sel] = pr["E0"]
        rate, _ = self._arr_rate(sigma, k0, eact, alpha, E0)

        sites, weights, pools = [], [], []
        # OC-type SEI diffuse into empty {BCO,OC} neighbours with nBA_SEI==0 & nLi>=2
        oc_sites = np.where(seiM & lat.is_OCsite() & (sigma > 0))[0]
        for i in oc_sites:
            dest = self._empty_neighbours(i, (BCO, OC),
                                          cond=lambda k: (self.nBA_SEI[k] == 0)
                                          & (self.nLi[k] >= 2))
            if dest.size and rate[i] > 0:
                sites.append(i); weights.append(rate[i] * dest.size); pools.append(dest)
        ba_sites = np.where(seiM & lat.is_BAsite() & (sigma > 0))[0]
        for i in ba_sites:
            dest = self._empty_neighbours(i, (BCB, BA),
                                          cond=lambda k: (self.nOC_SEI[k] == 0)
                                          & (self.nLi[k] >= 2))
            if dest.size and rate[i] > 0:
                sites.append(i); weights.append(rate[i] * dest.size); pools.append(dest)
        return sites, np.array(weights), pools, eact, rate

    def _empty_neighbours(self, i, rels, cond):
        lat = self.lat
        out = []
        for r in rels:
            ks = lat.nbr_indices[r][lat.nbr_indptr[r][i]:lat.nbr_indptr[r][i + 1]]
            ks = ks[self.occ[ks] != 1]
            if ks.size:
                ks = ks[cond(ks)]
            out.append(ks)
        return np.concatenate(out) if out else np.array([], dtype=int)

    def diffusion_step(self):
        sites, weights, pools, eact, rate = self._diffusion_candidates()
        if not sites or weights.sum() <= 0:
            return False
        W = weights.sum()
        dt = self._draw_dt(W)
        j = self.rng.choice(len(sites), p=weights / W)
        i = sites[j]
        k = pools[j][self.rng.integers(pools[j].size)]
        self._swap(i, k)
        self.last.update(type="Diffusion", Ospcs=self.spt.code2sym[self.spc[k]],
                         Ea=float(eact[k]), reactRate=float(rate[k] if rate[k] else 0),
                         time=float(dt), expValue=0.0)
        return True

    def _swap(self, a, b):
        self.spc[a], self.spc[b] = self.spc[b], self.spc[a]
        self.occ[a], self.occ[b] = self.occ[b], self.occ[a]
        self.chg[a], self.chg[b] = self.chg[b], self.chg[a]

    def _passivation_mask(self):
        """Boolean mask of electronically BLOCKED cathode sites (point 3).

        A site is blocked when >= pass_nmin of its neighbours (BCB, BCO, BA,
        OC shells, the same edge set as the nLi count) are occupied by a
        passivating species (e.g. Li2S, Li2S2). Evaluated on the current
        lattice state, so it is dynamic and needs no checkpoint entry.
        Returns None when the feature is off.
        """
        if self.pass_nmin <= 0:
            return None
        if self.pass_code_arr is None:
            codes = sorted({self.spt.add(sym) for sym in self.pass_species})
            self.pass_code_arr = np.array(codes, dtype=self.spc.dtype)
        passM = (self.occ == 1) & np.isin(self.spc, self.pass_code_arr)
        n_pass = self._count(self.eg4, passM)
        blocked = n_pass >= self.pass_nmin
        return blocked

    def _region_mask(self, region):
        """Boolean mask of sites inside a named z-band (anode / cathode / any)."""
        if region in ("any", "", None):
            return np.ones(self.N, dtype=bool)
        z = self.lat.frac[:, 2]
        if region == "cathode" and self.g.cathode_enabled:
            half = self.g.cathode_thickness / 2.0
            c = self.g.cathode_center
            return self._zband(z, c, half)
        if region == "anode":
            if self.g.anode_thickness_is_fraction:
                half = self.g.anode_thickness / 2.0
            else:
                half = (self.g.anode_thickness_A / self.lat.box[2, 2]) / 2.0
            c = self.g.anode_center
            return self._zband(z, c, half)
        if region == "anode_surface":
            # the live Li surface wherever it is (the anode grows past its
            # initial band): same criterion as LiStripping candidates
            return self._surface_li0_mask()
        if region == "cathode_surface":
            return self._cathode_surface_mask()
        return np.zeros(self.N, dtype=bool)

    def _cathode_surface_mask(self):
        """Electrolyte sites (ETH / SOL / FSI on OC or BA sites) with at least
        one neighbour occupied by cathode material (S8, Li2Sx on the lattice).
        This is where the cathode electrolyte interphase (CEI) forms. Dynamic:
        follows the cathode as it converts, dissolves and re-precipitates.
        CEI film species are NOT cathode material, so a film site does not by
        itself extend the surface (the film grows only next to active
        cathode sites)."""
        if not self.cath_code_arr.size:
            return np.zeros(self.N, dtype=bool)
        occ1 = self.occ == 1
        cathM = occ1 & np.isin(self.spc, self.cath_code_arr)
        n_cath = self._count(self.eg4, cathM)
        elec = occ1 & ((self.spc == self.cETH) | (self.spc == self.cSOL)
                       | (self.spc == self.cFSI))
        return elec & (n_cath >= 1)

    def sei_step(self):
        """One reaction event chosen by BKL across decomposition/plating/
        stripping/surface candidates."""
        self._refresh_counts()
        lat = self.lat
        occ1 = self.occ == 1
        liM = occ1 & (self.spc == self.cLi)
        li0 = liM & (self.chg == 0)
        liN = liM & (self.chg != 0)

        cand_types, cand_sites, cand_w, cand_dest = [], [], [], []

        atBegin = self.presentV == self.p.BeginV
        atEnd = self.presentV == self.p.EndV

        # --- decomposition / plating from FSI(salt) sites (BeginV) ----------
        if atBegin:
            sr_plate = self._decomp_sumrate("Plating")
            sr_fsi = self._decomp_sumrate("FSI")
            sr_psei = self._decomp_sumrate("PlatingSEI")
            if self.bv:
                # Butler-Volmer plating: one exchange rate for both Li placing
                # channels; DECOMPOSITION.in Plating/PlatingSEI entries ignored
                k_bv = self.bv_k0 * np.exp(-self.bv_alpha * self.eta_c / self.kT_V)
                sr_plate = sr_psei = k_bv * self.p.globalA
            salt = np.where(occ1 & (self.spc == self.cFSI))[0]
            # shared pool: every Li-placing event needs one Li+ and scales
            # with the remaining fraction of the pool (mean-field depletion)
            if self.pool_shared:
                fpool = (self.li_pool / self.li_pool0) if self.li_pool0 > 0 else 0.0
                if self.li_pool < 1 or fpool <= 0:
                    sr_plate = sr_fsi = sr_psei = 0.0
                else:
                    sr_plate *= fpool; sr_fsi *= fpool; sr_psei *= fpool
            if self.bv and sr_plate > 0:
                # BV plating: Li+ is available everywhere in a 1.2 M electrolyte,
                # so the plating weight is k_bv times the number of growth sites
                # (empty BC next to electrolyte and to 2 to 5 Li), not the number
                # of salt sites touching Li0. Removes the legacy artifact where a
                # slow side reaction fires whenever no salt site is adjacent.
                # LiPlating() then picks one of those growth sites at random.
                elig = self._eligible_BC() & (self.nETH > 0) & (self.nLi >= 2) & (self.nLi <= 5)
                n_el = int(elig.sum())
                if n_el:
                    cand_types.append("Plating"); cand_sites.append(-1)
                    cand_w.append(n_el * sr_plate); cand_dest.append(-1)
                sr_plate = sr_psei = 0.0     # no salt-based plating candidates
            if salt.size:
                mult_bc0 = self._neighbour_mult(salt, (BCB, BCO), li0)
                mult_allN = self._neighbour_mult(salt, (BCB, BCO, BA, OC, TE), liN)
                for s, mbc, mall in zip(salt, mult_bc0, mult_allN):
                    if mbc > 0 and sr_plate > 0:
                        cand_types.append("Plating"); cand_sites.append(s)
                        cand_w.append(mbc * sr_plate); cand_dest.append(-1)
                    if mbc > 0 and sr_fsi > 0:
                        cand_types.append("FSI"); cand_sites.append(s)
                        cand_w.append(mbc * sr_fsi); cand_dest.append(-1)
                    if mall > 0 and sr_psei > 0:
                        cand_types.append("PlatingSEI"); cand_sites.append(s)
                        cand_w.append(mall * sr_psei); cand_dest.append(-1)

            # --- solvent decomposition (BeginV) -----------------------------
            sr_sol = self._decomp_sumrate("SOL")
            sr_sol2 = self._decomp_sumrate("SOL2")
            sol = np.where(occ1 & (self.spc == self.cSOL))[0]
            if sol.size:
                mult_bc0 = self._neighbour_mult(sol, (BCB, BCO), li0)
                oxM = occ1 & (self.spc == self.cO)
                mult_o = self._neighbour_mult(sol, (BA,), oxM)
                for s, mbc, mo in zip(sol, mult_bc0, mult_o):
                    if mbc > 0 and sr_sol > 0:
                        cand_types.append("SOL"); cand_sites.append(s)
                        cand_w.append(mbc * sr_sol); cand_dest.append(-1)
                    if mo > 0 and sr_sol2 > 0:
                        cand_types.append("SOL2"); cand_sites.append(s)
                        cand_w.append(mo * sr_sol2); cand_dest.append(-1)

            # --- SFO / F5D further decomposition (BeginV) -------------------
            for sym, key in (("SFO", "SFO"), ("F5D", "F5D")):
                code = self.spt.sym2code.get(sym)
                sr = self._decomp_sumrate(key)
                if code is not None and sr > 0:
                    sites = np.where(occ1 & (self.spc == code))[0]
                    for s in sites:
                        cand_types.append(key); cand_sites.append(s)
                        cand_w.append(sr); cand_dest.append(-1)

        # --- Li stripping / surface diffusion (EndV) ------------------------
        if atEnd:
            eact = self._eact_all()
            strip = self.mob.params.get("LiStripping", {}).get("Li")
            surf = self.mob.params.get("LiSurface", {}).get("Li")
            # li_pool_max: a full pool cannot accept more Li+ -> no stripping
            # (surface diffusion is unaffected). Legacy path when cap is off.
            if self._pool_full() and strip:
                strip = None
                self.pool_full_blocks += 1
            cand_li = np.where(li0 & (self.nLi <= 4) & (self.nETH >= 1))[0]
            if self.bv and strip and cand_li.size:
                # Butler-Volmer stripping: absolute scale from j0, site
                # selection from the coordination-dependent Eact (relative)
                e_min = float(eact[cand_li].min())
                k_bv_d = self.bv_k0 * self.p.globalA * np.exp(self.bv_alpha * self.eta_d / self.kT_V)
            for i in cand_li:
                if strip and strip["sigma"] != 0:
                    if self.bv:
                        r = k_bv_d * np.exp(-(eact[i] - e_min) / self.kT)
                    else:
                        r, _ = self._arr_rate(strip["sigma"], strip["k0"], eact[i],
                                              strip["alpha"], strip["E0"])
                    if r > 0:
                        cand_types.append("LiStripping"); cand_sites.append(i)
                        cand_w.append(r); cand_dest.append(-1)
                if surf and surf["sigma"] != 0:
                    # eligible empty BCB (nLi<=5) and BCO (nLi>=4) neighbours
                    d1 = self._empty_neighbours(i, (BCB,),
                                                cond=lambda k: (self.nLi[k] <= 5)
                                                & (self.nETH[k] >= 1))
                    d2 = self._empty_neighbours(i, (BCO,),
                                                cond=lambda k: (self.nLi[k] >= 4)
                                                & (self.nETH[k] >= 1))
                    dest = np.concatenate([d1, d2])
                    if dest.size:
                        r, _ = self._arr_rate(surf["sigma"], surf["k0"], eact[i],
                                              surf["alpha"], surf["E0"])
                        if r > 0:
                            cand_types.append("LiSurface"); cand_sites.append(i)
                            cand_w.append(r * dest.size)
                            cand_dest.append(dest)

        # --- cathode / conversion / reservoir reactions (region + voltage) ---
        blockedM = None   # passivation mask, evaluated at most once per step
        n_dep_room = None # eligible deposit sites, evaluated at most once per step
        for r in self.conv_rxns:
            if r.voltage == "begin" and not atBegin:
                continue
            if r.voltage == "end" and not atEnd:
                continue
            # reservoir prerequisites (REQUIRE <spc_d> <n>)
            if any(self.n_dis.get(spc, 0) < n for (spc, n) in r.requires):
                continue
            # shared pool: a reaction needing n Li+ is ineligible below n
            if self.pool_shared and self.li_pool < self._li_need.get(r.name, 0):
                continue
            # li_pool_max: a reaction releasing Li+ is ineligible at the cap
            if self._li_release.get(r.name, 0) > 0 and self._pool_full():
                continue
            # deposit room: a reaction placing n DEPOSIT sites is ineligible
            # when fewer than n eligible BA surface sites exist (S conservation)
            dneed = self._dep_need.get(r.name, 0)
            if dneed > 0:
                if n_dep_room is None:
                    n_dep_room = self._deposit_eligible_count()
                if n_dep_room < dneed:
                    continue
            sr = self._decomp_sumrate(
                r.name, self._cathode_V() if (self.cath_bv and r.region == "cathode") else None)
            if sr <= 0:
                continue
            # mean-field scaling by reservoir concentration (RATE_SCALE)
            if r.rate_scale:
                if self.reservoir_sites <= 0:
                    continue
                sr *= self.n_dis.get(r.rate_scale, 0) / float(self.reservoir_sites)
                if sr <= 0:
                    continue
            reg = self._region_mask(r.region)
            if r.trigger == "EMPTY":
                sel = (self.occ != 1) & lat.is_BC() & reg
            else:
                tcode = self.spt.sym2code.get(r.trigger)
                if tcode is None:
                    continue
                sel = occ1 & (self.spc == tcode) & reg
            if r.region == "cathode" and self.cath_access is not None:
                sel &= self.cath_access
            # dynamic Li2S passivation: only Li-transfer reactions are gated
            if (r.region == "cathode" and self.pass_nmin > 0
                    and self._li_transfer.get(r.name, False)):
                if blockedM is None:
                    blockedM = self._passivation_mask()
                    self.n_cathode_blocked = int(
                        (blockedM & self._region_mask("cathode")
                         & lat.is_BC()).sum())
                sel &= ~blockedM
            sites = np.where(sel)[0]
            for s in sites:
                cand_types.append(("CONV", r.name)); cand_sites.append(s)
                cand_w.append(sr); cand_dest.append(-1)

        if not cand_w:
            self._last_W = 0.0
            return False
        w = np.array(cand_w, dtype=float)
        W = w.sum()
        self._last_W = float(W)
        if W <= 0:
            return False
        dt = self._draw_dt(W)
        j = int(self.rng.choice(len(w), p=w / W))
        rtype = cand_types[j]; site = cand_sites[j]
        is_conv = isinstance(rtype, tuple) and rtype[0] == "CONV"
        rname = rtype[1] if is_conv else rtype
        # C++ behaviour: a selected NON-diffusion event whose waiting time
        # exceeds scanInterval/5 is rejected and the loop is re-run (without
        # advancing the step). This caps the time a single slow event can skip.
        if dt > self.p.scanInterval / 5.0:
            return "REJECT"
        self.last.update(type=rname, Ospcs=(self.spt.code2sym[self.spc[site]] if site >= 0 else "Li"),
                         Ea=0.0, expValue=0.0, reactRate=float(w[j]), time=float(dt))

        # execute (and track whether this event transferred Li for end_half_when_idle)
        li_transfer = (rtype in ("Plating", "PlatingSEI", "LiStripping")
                       or (is_conv and self._li_transfer.get(rname, False)))
        self.idle_events = 0 if li_transfer else self.idle_events + 1
        if is_conv:
            self._apply_conversion(self.mech.get(rname), site)
            self.R[rname] += 1
        elif rtype in ("Plating", "PlatingSEI"):
            self.LiPlating(1); self.R[rtype] += 1
        elif rtype in ("FSI", "SFO", "SOL", "F5D", "SOL2"):
            self.electrolyteDecom(rtype)
        elif rtype == "LiStripping":
            self._empty_one(site); self.R["LiStripping"] += 1
            if self.pool_shared:
                self.li_pool += 1
        elif rtype == "LiSurface":
            dest = cand_dest[j]
            k = dest[self.rng.integers(dest.size)]
            self._swap(site, k); self.R["LiSurface"] += 1
        return True

    def _empty_one(self, i):
        self.occ[i] = 0; self.spc[i] = 0; self.chg[i] = 0

    def _neighbour_mult(self, sites, rels, value_mask):
        """For each site in `sites`, count neighbours (over rels) where value_mask True."""
        out = np.zeros(len(sites), dtype=int)
        lat = self.lat
        for n_idx, i in enumerate(sites):
            c = 0
            for r in rels:
                ks = lat.nbr_indices[r][lat.nbr_indptr[r][i]:lat.nbr_indptr[r][i + 1]]
                if ks.size:
                    c += int(value_mask[ks].sum())
            out[n_idx] = c
        return out

    # --------------------------------------------------------------- step
    def step(self):
        """One call of calculateAllReactions: either a diffusion event or a
        SEI/reaction event, with electrolyte/Li-framework bookkeeping."""
        self._voltage_update()

        fired = False
        attempts = 0
        reruns = 0
        max_attempts = int(getattr(self.p, "stall_attempts", 200))
        max_reruns = int(getattr(self.p, "stall_reruns", 10000))
        p_min = float(getattr(self.p, "stall_p_accept_min", 1e-4))
        while not fired and attempts < max_attempts and reruns < max_reruns:
            if self.stop < self.MaxToStop:
                self.stop += 1
                fired = self.diffusion_step()
                if not fired:
                    # refresh electrolyte and retry (mirrors NoRxn path)
                    self.flushElectrolyte(); self.updateEther()
                    self.addSolvent(); self.addLithiumSalt(); self.updateCharges()
                    attempts += 1
            else:
                self.stop = 0
                self.MaxToStop = int(self.rng.random() * 15) + 1
                fired = self.sei_step()
                if fired == "REJECT":
                    # slow non-diffusion event rejected: re-run the loop
                    # (stop is now 0, so diffusion is attempted next) without
                    # advancing the step or counting a failed attempt.
                    fired = False
                    reruns += 1
                    # hopeless: the expected waiting time 1/W is so far beyond
                    # scanInterval/5 that (almost) every draw will be rejected
                    p_acc = 1.0 - np.exp(-self._last_W * self.p.scanInterval / 5.0)
                    if p_acc < p_min:
                        reruns = max_reruns
                    continue
                if not fired:
                    self.flushElectrolyte(); self.updateEther()
                    self.addSolvent(); self.addLithiumSalt(); self.updateCharges()
                    attempts += 1
        if attempts >= max_attempts or reruns >= max_reruns:
            self._log_stall(attempts, reruns)
            self.stalled = True
            return

        # post-event framework maintenance
        self.flushElectrolyte()
        self.wrap()
        self.updateLiMetal()
        self.updateEther()
        self.addSolvent()
        self.addLithiumSalt()
        self.updateCharges()

        if self.last["type"] != "Diffusion":
            self.timeCurrentV += self.last["time"]
            self.currentTime += self.last["time"]
            self.countOandF()
            self._write_step()
            self.stepCurrentV += 1
            self.currentStep += 1

        if self.p.XYZprintFreq > 0 and self.currentStep % self.p.XYZprintFreq == 0:
            self._write_xyz(); self.countOandF()

    def _log_stall(self, attempts, reruns):
        """Diagnostic line for a run that cannot fire any event."""
        self._refresh_counts()
        elig = self._eligible_BC() & (self.nETH > 0) & (self.nLi >= 2) & (self.nLi <= 5)
        n_li = int(((self.occ == 1) & (self.spc == self.cLi)).sum())
        W = float(getattr(self, "_last_W", 0.0))
        wait = (1.0 / W) if W > 0 else float("inf")
        self.out.log(
            "Simulation stalled: no event can fire "
            f"(attempts={attempts}, rejected draws={reruns}) at step {self.currentStep}, "
            f"half-cycle {self.cycleNumber}, V={self.presentV}. Total non-diffusion rate "
            f"W={W:.3e} s^-1 (expected waiting time {wait:.3e} s, limit "
            f"{self.p.scanInterval / 5.0:.0f} s). Li+ pool={getattr(self, 'li_pool', 'n/a')}, "
            f"Li metal sites={n_li}, plating-eligible sites={int(elig.sum())}, "
            f"electrolyte sites={self._electrolyte_count()}. Stopping.")

    # --------------------------------------------------------------- logging
    def _row(self):
        return dict(
            stop=self.stop, MaxToStop=self.MaxToStop, currentStep=self.currentStep,
            cycleNumber=self.cycleNumber, stepCurrentV=self.stepCurrentV,
            timeCurrentV=self.timeCurrentV, currentTime=self.currentTime,
            presentV=self.presentV, type=self.last["type"], Ospcs=self.last["Ospcs"],
            Ea=self.last["Ea"], expValue=self.last["expValue"],
            reactRate=self.last["reactRate"], time=self.last["time"],
            RxnPlating=self.R["Plating"], RxnStripping=self.R["LiStripping"],
            RxnLiSurface=self.R["LiSurface"], RxnFSI=self.R["FSI"],
            RxnSFO=self.R["SFO"], RxnSOL=self.R["SOL"], RxnF5D=self.R["F5D"],
            RxnSOL2=self.R["SOL2"], RxnPlatingSEI=self.R["PlatingSEI"],
            LiMetal=self.LiMetal, LiIon=self.LiIon, LiMetalS=self.LiMetalS,
            LiIonS=self.LiIonS, nO=self.nO, nF=self.nF)

    def _write_step(self):
        self.out.write_step(self._row())

    def _write_xyz(self):
        self.out.write_xyz(self.currentStep, self.lat, self.occ, self.spc,
                           self.chg, self.spt.code2sym)

    def _write_snapshot(self):
        """Dump a POSCAR of the current frozen structure for Zeo++/RASPA.

        Only occupied, non-electrolyte sites (the SEI + metal skeleton) are
        written by default, since that is what pore/adsorption analysis cares
        about.  Files land in <out>/snapshots/snap_cycleNNN.vasp .
        """
        import os
        skel = {"ETH", "SOL", "FSI"}            # mobile electrolyte: excluded
        sym = self.spt.code2sym
        keep = (self.occ == 1) & np.array(
            [sym.get(int(c), "") not in skel for c in self.spc])
        d = os.path.join(self.out.root, "snapshots")
        os.makedirs(d, exist_ok=True)
        path = os.path.join(d, f"snap_cycle{self.cycleNumber:04d}.vasp")
        cart = self.lat.frac @ self.lat.box
        species = np.array([sym.get(int(c), "X") for c in self.spc])
        order = np.argsort(species[keep], kind="stable")
        ksel = np.where(keep)[0][order]
        uniq, counts = np.unique(species[ksel], return_counts=True)
        with open(path, "w") as fh:
            fh.write(f"kmc3d snapshot cycle {self.cycleNumber}\n1.0\n")
            for r in self.lat.box:
                fh.write(f"  {r[0]:.6f} {r[1]:.6f} {r[2]:.6f}\n")
            fh.write("  ".join(uniq) + "\n")
            fh.write("  ".join(str(int(c)) for c in counts) + "\n")
            fh.write("Cartesian\n")
            for i in ksel:
                p = cart[i]
                fh.write(f"  {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")
        self.out.log(f"snapshot -> {path} ({ksel.size} atoms)")

    # --------------------------------------------------------------- run
    def run(self, restart: str = ""):
        if restart:
            self.load_checkpoint(restart)
            self.out.log(f"Restarted from {restart}: step={self.currentStep} "
                         f"cycle={self.cycleNumber}")
        else:
            self.out.log(f"Starting kMC run with seed {self.p.seed}")
            self.create_cell()
            self.updateCharges()
            self._init_pool_and_reservoir()    # no RNG use
            self._write_xyz()                  # frame 0
        ckpt_every = int(getattr(self.p, "checkpointEveryCycles", 0))
        last_ckpt_cycle = self.cycleNumber
        last_stop_check = self.cycleNumber
        _anode_events = lambda: self.R["Plating"] + self.R["PlatingSEI"] + self.R["LiStripping"]
        anode_prev = _anode_events()
        anode_idle = 0
        try:
            while (self.currentStep <= self.p.totalSteps
                   and self.cycleNumber <= self.p.maxCycles):
                self.step()
                if self.stalled:
                    break
                if self.stop_anode_halves > 0 and self.cycleNumber != last_stop_check:
                    cur = _anode_events()
                    anode_idle = anode_idle + 1 if cur == anode_prev else 0
                    anode_prev = cur
                    if self.stop_elec_frac <= 0:
                        last_stop_check = self.cycleNumber
                    if anode_idle >= self.stop_anode_halves:
                        self.anode_dead = True
                        self.out.log(f"Anode inactive for {anode_idle} half-cycles (no plating, "
                                     f"no stripping) at half-cycle {self.cycleNumber}. Stopping.")
                        last_stop_check = self.cycleNumber
                        break
                if self.stop_elec_frac > 0 and self.cycleNumber != last_stop_check:
                    # evaluated once per half-cycle flip
                    last_stop_check = self.cycleNumber
                    n_el = self._electrolyte_count()
                    if n_el < self.stop_elec_frac * self.n_elec0:
                        self.dried_out = True
                        self.out.log(f"Cell dried out: {n_el} electrolyte sites "
                                     f"(< {self.stop_elec_frac:.2f} of {self.n_elec0}) "
                                     f"at half-cycle {self.cycleNumber}. Stopping.")
                        break
                if (ckpt_every and self.cycleNumber != last_ckpt_cycle
                        and self.cycleNumber % ckpt_every == 0):
                    self.save_checkpoint()
                    # flush the ledger with every checkpoint: a SLURM
                    # time-limit kill must not lose cycle_stats.csv. The file
                    # is rewritten in full each time (idempotent), so the final
                    # write at the end of the run is byte-identical to before.
                    self.ledger.finalize(self.out.root)
                    last_ckpt_cycle = self.cycleNumber
        finally:
            # always flush the ledger + a final checkpoint, even on interrupt
            path = self.ledger.finalize(self.out.root)
            if path:
                self.out.log(f"cycle stats -> {path}")
            if ckpt_every:
                self.save_checkpoint()
        self.out.log(f"Simulation terminated: steps={self.currentStep} "
                     f"cycles={self.cycleNumber}")

    def _init_pool_and_reservoir(self):
        """Initialise the shared Li+ pool and the reservoir normalisation.
        Deterministic (no RNG) so legacy runs are unaffected."""
        if self.reservoir_sites <= 0:
            self.reservoir_sites = int(((self.occ == 1) & (self.spc == self.cETH)).sum())
        if self.pool_shared:
            n0 = int(getattr(self.p, "li_pool_init", -1))
            if n0 < 0:
                self.etherVol()
                n0 = int(round(self.p.saltMolar * self.electrolyteVol
                               * MOLAR_TO_ATOMS_PER_A3))
            self.li_pool = self.li_pool0 = n0
            self.out.log(f"li_pool_mode=shared  li_pool0={n0}  "
                         f"reservoir_sites={self.reservoir_sites}  "
                         f"li_bulk_init={self.li_bulk_init}"
                         f"{' (unlimited foil)' if self.li_bulk_init < 0 else ''}")
        if self.pass_nmin > 0:
            self.out.log(f"cathode_passivation: species={','.join(self.pass_species)}  "
                         f"nmin={self.pass_nmin} (Li-transfer reactions gated)")
        if self.li_pool_max > 0:
            self.out.log(f"li_pool_max={self.li_pool_max}: stripping and RELEASE_LI "
                         f"gated when the pool is full")
        if self.bv:
            self.out.log(f"anode_kinetics=bv: k0_site={self.bv_k0} s^-1, alpha={self.bv_alpha}, "
                         f"eta_charge={self.eta_c} V, eta_discharge={self.eta_d} V "
                         f"(plating factor {np.exp(-self.bv_alpha * self.eta_c / self.kT_V):.2f}, "
                         f"stripping factor {np.exp(self.bv_alpha * self.eta_d / self.kT_V):.2f})")
        if self.presentV == self.p.EndV:
            self.out.log("first_half=end: the run starts with a discharge")
        if self.idle_limit > 0:
            self.out.log(f"end_half_when_idle={self.idle_limit}: half-cycles end after that "
                         f"many consecutive events without Li transfer")
        if self.cath_bv:
            self.out.log(f"cathode_kinetics=bv: V_charge={self.cath_V_c} V, "
                         f"V_discharge={self.cath_V_d} V (REGION cathode rates use "
                         f"alpha (V_cat - E0); alpha < 0 reduction, > 0 oxidation)")
        if self.stop_anode_halves > 0:
            self.out.log(f"stop_anode_inactive_halves={self.stop_anode_halves}: run stops "
                         f"after that many half-cycles without plating or stripping")
        if self.stop_elec_frac > 0:
            self.n_elec0 = self._electrolyte_count()
            self.out.log(f"stop_electrolyte_fraction={self.stop_elec_frac}: run stops "
                         f"below {int(self.stop_elec_frac * self.n_elec0)} of "
                         f"{self.n_elec0} electrolyte sites")

    # ------------------------------------------------------ checkpointing
    _CKPT_SCALARS = ("currentStep", "cycleNumber", "currentTime", "presentV",
                     "timeCurrentV", "stepCurrentV", "stop", "MaxToStop",
                     "LiMetal", "LiIon", "LiMetalS", "LiIonS", "nO", "nF",
                     "li_consumed", "li_released",
                     "li_pool", "li_pool0", "li_shuttled", "li_deposit",
                     "li_plated_pool", "reservoir_sites",
                     "li_bulk", "li_bulk_drawn", "li_bulk_returned", "s_cei", "li_cei",
                     "deposit_unplaced", "n_elec0", "pool_full_blocks")

    def save_checkpoint(self, path: str = ""):
        import json
        import os
        path = path or os.path.join(self.out.root, "checkpoint.npz")
        extra = {}
        if self.cath_access is not None:
            extra["cath_access"] = self.cath_access
        np.savez_compressed(path, occ=self.occ, spc=self.spc, chg=self.chg, **extra)
        meta = {k: getattr(self, k) for k in self._CKPT_SCALARS
                if hasattr(self, k)}
        meta["R"] = self.R
        meta["n_dis"] = self.n_dis
        meta["rng_state"] = self.rng.bit_generator.state
        meta["species_table"] = {s: int(c) for s, c in self.spt.sym2code.items()}
        meta["ledger_rows"] = self.ledger.rows
        meta["ledger_prevR"] = self.ledger._prevR
        meta["ledger_half_index"] = self.ledger._half_index
        with open(path + ".json", "w") as fh:
            json.dump(meta, fh)
        self.out.log(f"checkpoint -> {path}")

    def load_checkpoint(self, path: str):
        import json
        arr = np.load(path)
        self.occ[:] = arr["occ"]; self.spc[:] = arr["spc"]; self.chg[:] = arr["chg"]
        with open(path + ".json") as fh:
            meta = json.load(fh)
        # species codes must match the saved run for spc to be meaningful
        for s, c in meta["species_table"].items():
            self.spt.add(s)
            if self.spt.sym2code[s] != c:
                raise RuntimeError(
                    f"Species table mismatch on restart ('{s}'): the same "
                    f"MECHANISM/GEOMETRY inputs must be used.")
        for k in self._CKPT_SCALARS:
            if k in meta:
                setattr(self, k, meta[k])
        self.R.update(meta["R"])
        self.n_dis.update(meta.get("n_dis", {}))
        if "cath_access" in arr.files:
            self.cath_access = arr["cath_access"].astype(bool)
        self.rng.bit_generator.state = meta["rng_state"]
        self.ledger.rows = meta.get("ledger_rows", [])
        self.ledger._prevR = meta.get("ledger_prevR")
        self.ledger._half_index = meta.get("ledger_half_index", 0)
