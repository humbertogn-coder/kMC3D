"""
stats.py - per-half-cycle simulation ledger for ML-ready targets.

The engine calls CycleLedger.close_half() at every voltage flip. Each row of the
resulting cycle_stats.csv summarises the half-cycle that just ENDED:

  identity     : half_index, voltage of the half (V), charge/discharge label,
                 currentStep, sim time
  reactions    : per-reaction event counts DURING this half (columns rxn_*)
  CE           : coulombic-efficiency building blocks + a default CE per full
                 cycle (see below)
  Li inventory : n_Li, n_Li_ion, n_Li_dead (topologically disconnected),
                 n_Li_interface (Li touching SEI)
  SEI          : n_SEI total and per-species counts (sei_*), z-extent
                 (thickness), mean z
  electrolyte  : n_ETH, n_SOL, n_FSI
  morphology   : li_front_z, surface_roughness (std of per-column top Li z)
  cathode      : li_consumed / li_released deltas (full-cell runs)
  Li-S extras  : (only when active) li_pool, li_bulk, li_total (conservation
                 check, must stay constant), d_li_plated_pool, d_li_bulk_drawn,
                 d_li_bulk_returned, d_li_framework, n_dis_<spc>,
                 d_li_shuttled, d_li_deposit, n_cathode_blocked (Li2S passivation)
  CE columns   : CE_cycle (legacy), CE_mod (decomposition electrons in the
                 denominator), CE_shuttle (cathode-side, shuttle-dominated)

Default CE definition (documented, override downstream if you prefer another):
    CE_cycle = stripped_during_discharge_half / plated_during_charge_half
computed on consecutive (charge, discharge) half pairs; NaN when undefined.
Raw per-half reaction counts are always written so any convention can be
recomputed from the CSV.

Dead-Li metric: Li sites NOT connected (through the full neighbour graph,
Li-to-Li) to the anode anchor slab (the initial anode z-band). This is the
topological "stock" version of the buried/dead lithium descriptor.
"""

from __future__ import annotations

import os
import re
from typing import Dict, List, Optional

import numpy as np


class CycleLedger:
    def __init__(self, grid_xy: int = 32):
        self.rows: List[Dict] = []
        self._prevR: Optional[Dict[str, int]] = None
        self._prev_extra = {"li_consumed": 0, "li_released": 0,
                            "li_shuttled": 0, "li_deposit": 0,
                            "li_plated_pool": 0, "n_Li": None,
                            "li_bulk_drawn": 0, "li_bulk_returned": 0}
        self._half_index = 0
        self.grid_xy = grid_xy

    # ------------------------------------------------------------ engine API
    def close_half(self, eng) -> None:
        """Record the half-cycle that is ending now (presentV not yet flipped)."""
        R = dict(eng.R)
        if self._prevR is None:
            dR = dict(R)
        else:
            dR = {k: R.get(k, 0) - self._prevR.get(k, 0) for k in R}
        self._prevR = R

        row: Dict = {
            "half_index": self._half_index,
            "cycle": eng.cycleNumber,
            "voltage": float(eng.presentV),
            "phase": "charge" if eng.presentV == eng.p.BeginV else "discharge",
            "step": int(eng.currentStep),
            "time": float(getattr(eng, "currentTime", 0.0)),
        }
        for k, v in sorted(dR.items()):
            row[f"rxn_{k}"] = int(v)

        row.update(self._inventory(eng))

        # cathode bookkeeping deltas (0 for anode-only runs)
        for key in ("li_consumed", "li_released"):
            cur = int(getattr(eng, key, 0))
            row[f"d_{key}"] = cur - self._prev_extra[key]
            self._prev_extra[key] = cur

        # ---- shared Li+ pool / well-mixed polysulfide reservoir (opt-in) ----
        # Columns only appear when the feature is active, so legacy CSVs keep
        # exactly their old column set.
        n_dis = getattr(eng, "n_dis", {})
        if n_dis:
            for spc, n in sorted(n_dis.items()):
                row[f"n_dis_{spc}"] = int(n)
            for key in ("li_shuttled", "li_deposit"):
                cur = int(getattr(eng, key, 0))
                row[f"d_{key}"] = cur - self._prev_extra[key]
                self._prev_extra[key] = cur
        # ---- sulfur balance (only when the polysulfide reservoir is active) --
        # s_total = S on lattice cathode species + S in anode deposits
        #         + S dissolved in the reservoir + S sequestered in the CEI.
        # Must be CONSTANT over the run (= 8 * initial S8 sites).
        if n_dis:
            row["s_total"] = self._sulfur_total(eng, n_dis)
            row["deposit_unplaced"] = int(getattr(eng, "deposit_unplaced", 0))
        # ---- cathode electrolyte interphase (opt-in: CEI keyword) -----------
        if getattr(eng, "cei_active", False):
            row["s_cei"] = int(getattr(eng, "s_cei", 0))
            row["li_cei"] = int(getattr(eng, "li_cei", 0))
            occ1 = eng.occ == 1
            ceiM = occ1 & np.isin(eng.spc, eng.cei_code_arr)
            row["n_CEI"] = int(ceiM.sum())
            if ceiM.any():
                codes, counts = np.unique(eng.spc[ceiM], return_counts=True)
                for c, n in zip(codes, counts):
                    row[f"cei_{eng.spt.code2sym[int(c)]}"] = int(n)
        # ---- dynamic Li2S passivation of the cathode (opt-in, point 3) ----
        if getattr(eng, "pass_nmin", 0) > 0:
            bm = eng._passivation_mask()
            row["n_cathode_blocked"] = int(
                (bm & eng._region_mask("cathode") & eng.lat.is_BC()).sum())
        if getattr(eng, "pool_shared", False):
            row["li_pool"] = int(eng.li_pool)
            if getattr(eng, "li_pool_max", -1) > 0:
                row["pool_full_blocks"] = int(eng.pool_full_blocks)
            cur = int(getattr(eng, "li_plated_pool", 0))
            row["d_li_plated_pool"] = cur - self._prev_extra["li_plated_pool"]
            self._prev_extra["li_plated_pool"] = cur
            # Anode framework Li (wrap / updateLiMetal) now comes from the
            # bulk foil reservoir: d_li_framework = drawn - returned per half.
            for key in ("li_bulk_drawn", "li_bulk_returned"):
                cur = int(getattr(eng, key, 0))
                row[f"d_{key}"] = cur - self._prev_extra[key]
                self._prev_extra[key] = cur
            row["d_li_framework"] = row["d_li_bulk_drawn"] - row["d_li_bulk_returned"]
            row["li_bulk"] = int(getattr(eng, "li_bulk", 0))
            # Li conservation check: pool + foil + Li0 sites + Li bound in
            # polysulfides (cathode-consumed net + shuttle-corroded). Must be
            # CONSTANT across the run; any drift is a leak.
            row["li_total"] = (row["li_pool"] + row["li_bulk"] + row["n_Li"]
                               + int(getattr(eng, "li_consumed", 0))
                               - int(getattr(eng, "li_released", 0))
                               + int(getattr(eng, "li_shuttled", 0)))
            self._prev_extra["n_Li"] = row["n_Li"]

        self.rows.append(row)
        self._half_index += 1

    # ------------------------------------------------------------ metrics
    @staticmethod
    def _s_atoms(sym: str) -> int:
        """Sulfur atoms in a lattice/reservoir species label: S8 -> 8,
        Li2S8 -> 8, Li2S -> 1, Li2S2_an -> 2, Li2S4_d -> 4. Labels whose S
        is followed by another element letter (SOL, SFO) count as 0; CEI
        film species are tracked through the S_LOSS counter instead."""
        m = re.search(r"S(\d*)(?=_|$)", sym)
        if not m:
            return 0
        return int(m.group(1)) if m.group(1) else 1

    def _sulfur_total(self, eng, n_dis) -> int:
        occ1 = eng.occ == 1
        tot = 0
        codes = set(int(c) for c in getattr(eng, "cath_code_arr", []))
        codes |= set(int(c) for c in getattr(eng, "dep_code_arr", []))
        codes -= set(int(c) for c in getattr(eng, "cei_code_arr", []))
        for c in codes:
            ns = self._s_atoms(eng.spt.code2sym[c])
            if ns:
                tot += ns * int((occ1 & (eng.spc == c)).sum())
        for spc, n in n_dis.items():
            tot += self._s_atoms(spc) * int(n)
        tot += int(getattr(eng, "s_cei", 0))
        return tot

    def _inventory(self, eng) -> Dict:
        occ1 = eng.occ == 1
        spc = eng.spc
        liM = occ1 & (spc == eng.cLi)
        ethM = occ1 & (spc == eng.cETH)
        solM = occ1 & (spc == eng.cSOL)
        fsiM = occ1 & (spc == eng.cFSI)
        cath_codes = getattr(eng, "cath_code_arr", np.array([]))
        cathM = (occ1 & np.isin(spc, cath_codes)) if cath_codes.size \
            else np.zeros(eng.N, dtype=bool)
        seiM = occ1 & ~liM & ~ethM & ~solM & ~fsiM & ~cathM

        cart_z = (eng.lat.frac @ eng.lat.box)[:, 2]
        out: Dict = {
            "n_Li": int(liM.sum()),
            "n_Li_ion": int((liM & (eng.chg != 0)).sum()),
            "n_ETH": int(ethM.sum()),
            "n_SOL": int(solM.sum()),
            "n_FSI": int(fsiM.sum()),
            "n_SEI": int(seiM.sum()),
            "n_cathode": int(cathM.sum()),
        }
        # per-species cathode composition (cat_S8, cat_Li2S8, ..., cat_Li2S)
        if cathM.any():
            codes, counts = np.unique(spc[cathM], return_counts=True)
            for c, n in zip(codes, counts):
                out[f"cat_{eng.spt.code2sym[int(c)]}"] = int(n)

        # per-species SEI composition
        if seiM.any():
            codes, counts = np.unique(spc[seiM], return_counts=True)
            for c, n in zip(codes, counts):
                out[f"sei_{eng.spt.code2sym[int(c)]}"] = int(n)
            zs = cart_z[seiM]
            out["sei_thickness"] = float(zs.max() - zs.min())
            out["sei_mean_z"] = float(zs.mean())
        else:
            out["sei_thickness"] = 0.0
            out["sei_mean_z"] = 0.0

        # Li interface (touching SEI) and topologically dead Li
        out["n_Li_interface"] = self._li_touching_sei(eng, liM, seiM)
        out["n_Li_dead"] = self._li_dead(eng, liM)

        # Li front + surface roughness on an xy grid
        if liM.any():
            zli = cart_z[liM]
            out["li_front_z"] = float(zli.max())
            out["surface_roughness"] = self._roughness(eng, liM, cart_z)
        else:
            out["li_front_z"] = 0.0
            out["surface_roughness"] = 0.0
        return out

    def _li_touching_sei(self, eng, liM, seiM) -> int:
        src, dst = eng.src, eng.dst
        touch = np.zeros(eng.N, dtype=bool)
        sel = seiM[dst]
        if sel.any():
            touch[np.unique(src[sel])] = True
        return int((liM & touch).sum())

    def _li_dead(self, eng, liM) -> int:
        """Li not connected (Li-Li bonds over all relations) to the anode band."""
        if not liM.any():
            return 0
        band = getattr(eng, "_anode_band", None)
        if band is None:
            return 0
        lo, hi = band
        z = (eng.lat.frac @ eng.lat.box)[:, 2]
        anchor = liM & (z >= lo) & (z <= hi)
        if not anchor.any():
            return int(liM.sum())          # everything is detached from anchor
        sel = liM[eng.src] & liM[eng.dst]
        es, ed = eng.src[sel], eng.dst[sel]
        reach = anchor.copy()
        while True:
            new_dst = ed[reach[es]]
            frontier = np.zeros(eng.N, dtype=bool)
            frontier[new_dst] = True
            frontier &= ~reach
            if not frontier.any():
                break
            reach |= frontier
        return int((liM & ~reach).sum())

    def _roughness(self, eng, liM, cart_z) -> float:
        g = self.grid_xy
        xy = eng.lat.frac[liM][:, :2] % 1.0
        cols = (np.minimum((xy * g).astype(int), g - 1))
        key = cols[:, 0] * g + cols[:, 1]
        top = np.full(g * g, -np.inf)
        np.maximum.at(top, key, cart_z[liM])
        top = top[np.isfinite(top)]
        return float(top.std()) if top.size else 0.0

    # ------------------------------------------------------------ output
    def finalize(self, out_root: str) -> Optional[str]:
        if not self.rows:
            return None
        # union of keys, stable order: identity first, then sorted rest
        head = ["half_index", "cycle", "voltage", "phase", "step", "time"]
        rest = sorted({k for r in self.rows for k in r} - set(head))
        cols = head + rest

        # default CE per full cycle: pair (charge, discharge) halves
        ce = self._default_ce()
        ce_mod = self._ce_mod()
        ce_li = self._ce_shuttle()

        path = os.path.join(out_root, "cycle_stats.csv")
        with open(path, "w") as fh:
            fh.write(",".join(cols + ["CE_cycle", "CE_mod", "CE_shuttle"]) + "\n")
            for i, r in enumerate(self.rows):
                vals = [str(r.get(c, "")) for c in cols]
                for arr in (ce, ce_mod, ce_li):
                    vals.append("" if np.isnan(arr[i]) else f"{arr[i]:.6f}")
                fh.write(",".join(vals) + "\n")
        return path

    # one electron per decomposition event (FSI plates 1 Li into the SEI;
    # SFO/SOL/SOL2/F5D are electrolyte reductions)
    _DECOMP = ("rxn_FSI", "rxn_SFO", "rxn_SOL", "rxn_SOL2", "rxn_F5D")

    def _ce_mod(self) -> np.ndarray:
        """CE_mod = stripped(discharge) / [plated + decomposition](charge).
        Every decomposition event costs one electron that never comes back,
        so it is charge spent during the charge half without Li recovered."""
        n = len(self.rows)
        ce = np.full(n, np.nan)
        for i in range(1, n):
            a, b = self.rows[i - 1], self.rows[i]
            if a.get("phase") == "charge" and b.get("phase") == "discharge":
                q_ch = (a.get("rxn_Plating", 0) + a.get("rxn_PlatingSEI", 0)
                        + sum(a.get(k, 0) for k in self._DECOMP))
                q_dis = b.get("rxn_LiStripping", 0)
                if q_ch > 0:
                    ce[i] = q_dis / q_ch
        return ce

    def _ce_shuttle(self) -> np.ndarray:
        """CE_shuttle = Q_useful / (Q_useful + Q_shuttle) per full cycle, the
        Mikhaylik-Akridge picture: charge passed on charging is either cathode
        oxidation (Li released, 1 e- each) or wasted corroding Li0 by the
        polysulfide shuttle (1 e- per Li0). Shuttle counted over both halves
        (self-discharge included). NaN if the shuttle columns are absent or no
        cathode oxidation happened."""
        n = len(self.rows)
        ce = np.full(n, np.nan)
        for i in range(1, n):
            a, b = self.rows[i - 1], self.rows[i]
            if a.get("phase") == "charge" and b.get("phase") == "discharge":
                if "d_li_shuttled" not in a and "d_li_shuttled" not in b:
                    continue
                q_use = a.get("d_li_released", 0)
                q_sh = a.get("d_li_shuttled", 0) + b.get("d_li_shuttled", 0)
                if q_use + q_sh > 0:
                    ce[i] = q_use / (q_use + q_sh)
        return ce

    def _default_ce(self) -> np.ndarray:
        n = len(self.rows)
        ce = np.full(n, np.nan)
        for i in range(1, n):
            a, b = self.rows[i - 1], self.rows[i]
            if a.get("phase") == "charge" and b.get("phase") == "discharge":
                plated = a.get("rxn_Plating", 0) + a.get("rxn_PlatingSEI", 0)
                stripped = b.get("rxn_LiStripping", 0)
                if plated > 0:
                    ce[i] = stripped / plated
        return ce
