"""
mechanism.py
============
Data-driven reaction PRODUCTS.  This externalizes everything that used to be
hard-coded inside `electrolyteDecom()` in calculate_reactions.cpp.

A reaction TYPE (e.g. FSI, SFO, SOL, F5D, SOL2, Plating, PlatingSEI) owns:
  - a SELECT rule deciding which product channel fires
  - one or more CHANNELS, each a list of product ACTIONS

ACTIONS
-------
  PLATE_LI                 -> LiPlating(1)            (place one Li0 at a BC site)
  PACK_OC <spc> <q> <n>    -> packOC(spc, q, n)       (place n <spc> at OC sites)
  PACK_BA <spc> <q> <n>    -> packBA(spc, q, n)       (place n <spc> at BA sites)
  CONVERT <spc> [<q>]      -> trigger site becomes <spc> (cathode conversion)
  CONSUME_LI / RELEASE_LI n-> Li+ pool debit / credit (counters in fixed mode)

Well-mixed reservoir / shuttle extensions (Li-S full cell, see MECHANISM.in of
the fullcell_shuttle case):
  header keywords : REGION anode_surface | cathode ; TRIGGER EMPTY ;
                    REQUIRE <spc_d> <n> ; RATE_SCALE <spc_d>
  actions         : RES <spc_d> <delta> ; DISSOLVE <spc_d> ; PRECIP <spc> [<q>] ;
                    STRIP_LI0 <n> ; DEPOSIT <spc> <q> <n>

SELECT modes
------------
  SINGLE                   -> always channel 0
  RANDINT_INCLUSIVE <N>    -> draw k = int(rand*(N+1)) in 0..N and fire channel k.
                              Channels are addressed by an explicit INDEX; any
                              missing index (e.g. k==0) is a NULL event (the
                              reaction still counts but produces nothing) - this
                              faithfully reproduces the C++ `n==0` no-op branch.

FILE FORMAT (MECHANISM.in)
--------------------------
    REACTION FSI
      SELECT SINGLE
      CHANNEL 0
        PLATE_LI
        PACK_OC F   -1 1
        PACK_OC SFO -1 1
    END

    REACTION SFO
      SELECT RANDINT_INCLUSIVE 7
      CHANNEL 1
        PACK_OC SFO -1 2
      CHANNEL 2
        PACK_OC SFO -1 1
        PACK_OC F   -1 1
      ... (channels 3..7) ...
    END
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Tuple
import os


Action = Tuple[str, str, int, int]   # (op, species, charge, count)


@dataclass
class Reaction:
    name: str
    select: str = "SINGLE"            # SINGLE | RANDINT_INCLUSIVE
    select_n: int = 0                 # N for RANDINT_INCLUSIVE
    channels: Dict[int, List[Action]] = field(default_factory=dict)
    # --- cathode / region-aware extensions (default = legacy anode behaviour) --
    region: str = "any"               # any | anode | cathode | anode_surface |
                                      #   cathode_surface (electrolyte next to cathode)
    weights: Dict[int, float] = field(default_factory=dict)   # SELECT WEIGHTED
    voltage: str = "any"              # any | begin | end       (half-cycle restriction)
    trigger: str = ""                 # species that triggers a CONVERSION reaction;
                                      #   empty -> legacy product-only reaction
                                      #   "EMPTY" -> fires on an empty BC site of
                                      #   the region (re-precipitation)
    # --- well-mixed polysulfide reservoir extensions (Li-S shuttle) ---------
    requires: List[Tuple[str, int]] = field(default_factory=list)
                                      # REQUIRE <spc_d> <n>: reservoir must hold
                                      #   >= n of the dissolved species
    rate_scale: str = ""              # RATE_SCALE <spc_d>: multiply the rate by
                                      #   the reservoir concentration
                                      #   N_dis[spc_d] / reservoir_sites (mean field)

    @property
    def is_conversion(self) -> bool:
        return bool(self.trigger)

    def pick_channel(self, rng) -> List[Action]:
        if self.select == "SINGLE":
            return self.channels.get(0, [])
        if self.select == "RANDINT_INCLUSIVE":
            k = int(rng.random() * (self.select_n + 1))   # 0..N inclusive
            return self.channels.get(k, [])               # missing -> null event
        if self.select == "WEIGHTED":
            # channels carry explicit weights: CHANNEL <idx> <weight>
            keys = sorted(self.channels)
            w = [max(self.weights.get(k, 0.0), 0.0) for k in keys]
            tot = sum(w)
            if tot <= 0:
                return []
            r = rng.random() * tot
            acc = 0.0
            for k, wk in zip(keys, w):
                acc += wk
                if r <= acc:
                    return self.channels[k]
            return self.channels[keys[-1]]
        raise ValueError(f"Unknown SELECT mode {self.select} in reaction {self.name}")


@dataclass
class Mechanism:
    reactions: Dict[str, Reaction] = field(default_factory=dict)

    def __contains__(self, name: str) -> bool:
        return name in self.reactions

    def get(self, name: str) -> Reaction:
        return self.reactions[name]

    @classmethod
    def read(cls, path: str) -> "Mechanism":
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mechanism file not found: {path}")
        mech = cls()
        cur: Reaction | None = None
        cur_ch: int | None = None
        with open(path) as fh:
            for raw in fh:
                line = raw.split("#")[0].strip()
                if not line:
                    continue
                parts = line.split()
                kw = parts[0].upper()

                if kw == "REACTION":
                    cur = Reaction(name=parts[1]); cur_ch = None
                    mech.reactions[cur.name] = cur
                elif kw == "REGION":
                    cur.region = parts[1].lower()
                elif kw == "VOLTAGE":
                    cur.voltage = parts[1].lower()
                elif kw == "TRIGGER":
                    cur.trigger = parts[1]
                elif kw == "SELECT":
                    cur.select = parts[1].upper()
                    cur.select_n = int(parts[2]) if len(parts) > 2 else 0
                elif kw == "CHANNEL":
                    cur_ch = int(parts[1]) if len(parts) > 1 else 0
                    cur.channels.setdefault(cur_ch, [])
                    if len(parts) > 2:                 # optional channel weight
                        cur.weights[cur_ch] = float(parts[2])
                elif kw == "END":
                    cur, cur_ch = None, None
                elif kw == "PLATE_LI":
                    cur.channels[cur_ch].append(("PLATE_LI", "Li", 0, 1))
                elif kw == "CONVERT":
                    # transform the trigger site to a new species: CONVERT <spc> [<q>]
                    spc = parts[1]; q = int(parts[2]) if len(parts) > 2 else 0
                    cur.channels[cur_ch].append(("CONVERT", spc, q, 1))
                elif kw == "CONSUME_LI":
                    cur.channels[cur_ch].append(("CONSUME_LI", "Li", 0, int(parts[1])))
                elif kw == "RELEASE_LI":
                    cur.channels[cur_ch].append(("RELEASE_LI", "Li", 0, int(parts[1])))
                # ---- reservoir / shuttle keywords (Li-S full cell) ----------
                elif kw == "REQUIRE":
                    cur.requires.append((parts[1], int(parts[2])))
                elif kw == "RATE_SCALE":
                    cur.rate_scale = parts[1]
                elif kw == "RES":
                    # RES <spc_d> <delta>: change the reservoir count
                    cur.channels[cur_ch].append(("RES", parts[1], 0, int(parts[2])))
                elif kw == "DISSOLVE":
                    # DISSOLVE <spc_d>: empty the trigger site, reservoir[spc_d] += 1
                    cur.channels[cur_ch].append(("DISSOLVE", parts[1], 0, 1))
                elif kw == "PRECIP":
                    # PRECIP <spc> [<q>]: fill the (empty) trigger site with <spc>
                    spc = parts[1]; q = int(parts[2]) if len(parts) > 2 else 0
                    cur.channels[cur_ch].append(("PRECIP", spc, q, 1))
                elif kw == "STRIP_LI0":
                    # STRIP_LI0 <n>: remove the trigger Li0 site plus n-1 more
                    # surface Li0 sites (anode corrosion by the shuttle)
                    cur.channels[cur_ch].append(("STRIP_LI0", "Li", 0, int(parts[1])))
                elif kw == "DEPOSIT":
                    # DEPOSIT <spc> <q> <n>: place n insoluble <spc> on BA sites
                    # at the Li surface (passivating deposit on the anode)
                    spc = parts[1]; q = int(parts[2]); n = int(parts[3])
                    cur.channels[cur_ch].append(("DEPOSIT", spc, q, n))
                # ---- cathode electrolyte interphase (CEI) keywords -----------
                elif kw == "CEI":
                    # CEI <spc> [<q>]: convert the trigger electrolyte site (an
                    # OC/BA site next to cathode material) into the inert film
                    # species <spc>. Film species are excluded from the anode
                    # SEI class (no Li framework growth around them) and can be
                    # listed in cathode_passivating_species.
                    spc = parts[1]; q = int(parts[2]) if len(parts) > 2 else 0
                    cur.channels[cur_ch].append(("CEI", spc, q, 1))
                elif kw == "S_LOSS":
                    # S_LOSS <n>: n sulfur atoms leave the active inventory
                    # (sequestered in the CEI). Keeps the S balance closed.
                    cur.channels[cur_ch].append(("S_LOSS", "S", 0, int(parts[1])))
                elif kw in ("PACK_OC", "PACK_BA"):
                    spc = parts[1]; q = int(parts[2]); n = int(parts[3])
                    cur.channels[cur_ch].append((kw, spc, q, n))
                else:
                    raise ValueError(f"Unknown MECHANISM keyword: {parts[0]}")
        return mech
