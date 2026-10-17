# Full cell: where each kinetic value lives and what can replace it

Status 2026-10-13. Case: cases/fullcell_physical (25 reactions). Purpose: when
new values arrive (AIMD of F5DEE + LiFSI + Li2Sx, DFT barriers, experiments),
each one replaces ONE line, keeps its tag in KINETICS_TABLE.md, and is tested
alone. MOBILITY.in pair interactions are [G] and are never changed.

## 1. Conversion rules (one formula per kind of input)

| Input you get | Converts to | Formula (298 K, kT = 0.02568 eV) |
|---|---|---|
| Elementary barrier Ea (NEB, metadynamics, AIMD) | per-site rate | k = (kT/h) exp(-Ea/kT), kT/h = 6.21e12 s^-1; +0.059 eV = 10x slower |
| Reaction free energy dG (forward vs reverse) | ratio of a reaction pair | k_f / k_b = exp(-dG/kT) (detailed balance) |
| Solvation free energy of Li2Sx (dG_solv) | dissolution / precipitation ratio | k_diss / k_prec = exp(-dG_solv/kT) (times the reference concentration) |
| Marcus reorganisation energy lambda | electron-transfer barrier | dG* = (lambda + dG)^2 / (4 lambda) |
| Diffusion coefficient D (AIMD MSD) | hop rate / shuttle flux | Gamma = 6 D / a^2 (3D lattice); flux to the anode J = D c / L |
| Exchange current density j0 or i0 | per-site rate | k = j0 A_site / (n e), A_site = 4e-16 cm^2; 1 mA/cm^2 = 2.5 s^-1 |
| Products seen in an AIMD trajectory | MECHANISM.in channel | new or edited REACTION block; rate from one of the rows above |

Rates from the first row are elementary (per attempt at one site). An
elementary rate used directly can make a site an inexhaustible sink (SOL2, run
1: kT/h with 0.364 eV gave 4e6 s^-1). Where an effective electrode rate exists
(CE anchor, j0), the barrier enters as a RELATIVE factor exp(-dEa/kT) on that
anchor, as done for SOL2.

## 2. Map of the current values

| Reaction (DECOMPOSITION.in) | Value now | Tag | Replace with |
|---|---|---|---|
| Plating / stripping (bv_k0_site) | 74 s^-1 | [V] Boyle 2020 | j0 for LiFSI/F5DEE if measured; AIMD desolvation barrier as a check (0.646 eV equivalent) |
| FSI, SFO, SOL | 1.118, 0.01118, 0.01471 | [V-derived, CE-anchored] | AIMD: sequence and products of LiFSI and F5DEE reduction on Li (MECHANISM); relative barriers between channels; keep the CE anchor for the absolute scale |
| SOL2 (F5DEE on Li2O) | 0.01471 with Ea 0.364 | [V-derived] | AIMD/DFT barrier on Li2O vs on Li (relative factor) |
| cat_*_red / cat_*_ox | 2.5 / 1.25, E0 2.39 to 2.01 V | [V] i0, [V-assumed] E0 | not from AIMD; experimental i0 / OCV of the reference cell |
| cat_diss_Li2S8 / Li4S8 | 0.05 / 0.02 s^-1 | [P] | AIMD: desorption barrier of Li2Sx from the cathode surface into F5DEE + LiFSI (row 1) |
| cat_prec_Li2S8 / S6 / S4 | 5.0 | [P] | AIMD: dG_solv of Li2Sx (row 3) fixes k_diss / k_prec; calibrate the ratio, not each value |
| sh_Li2S8 / S6 / S4 | 0.01 / 0.01 / 0.004 | [V-derived] Mikhaylik 2004 ks | AIMD: D of Li2Sx in F5DEE (row 5) and the reduction barrier of Li2Sx on Li (row 1; likely barrierless) |
| cei_FSI_Li2S8 / S6 | 1.0 | [P] | AIMD: does Li2Sx react with FSI-? products and barrier (rows 1, 7) |
| cei_F5DEE_Li2S8 | sigma 0 (off) | no evidence | AIMD evidence of F5DEE + Li2Sx reactivity: set sigma 1 and the barrier |
| reservoir_sites, li_pool_max | 2000, 148 | structural | not kinetic; keep |

## 3. Time scale (decide before plugging per-second rates)

The box cycles ~0.16 uAh/cm^2 per half-cycle, so at mA/cm^2 a half lasts
~0.03 s of kMC time, while a real C/10 half lasts ~3.6e4 s. Any rate in s^-1
for a chemical step (dissolution, precipitation, shuttle, CEI) therefore fires
~1e6 times less per cycle than in the real cell. Two consistent ways:

1. Per-cycle ratio (same logic as the original engine): set each chemical rate
   so that its events per half-cycle, relative to the electrochemical events,
   match the measured ratio (shuttle: f_c = k_s q_H [S] / I; dissolution: the
   dissolved fraction at the solubility limit). The AIMD values then fix the
   RATIOS between chemical channels.
2. Time-scale factor (opt-in keyword, not implemented): multiply the
   chemical rates by S = t_half_real / t_half_kMC. Exact per cycle, but the
   chemical events then dominate the event budget.

Until this is decided, bv_strip_symmetric stays off in the full cell.
