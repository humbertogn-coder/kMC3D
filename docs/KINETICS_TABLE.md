# Kinetics table for the physically based engine (anode_physical, cathode)

Status tags used below
    [V]  verified against a PDF in hand (page cited)
    [M]  from memory of the literature: value and citation MUST be verified
         against the PDF before it enters any .in file (none left, 2026-09-22)
    [V-assumed] verified in the PDF, but the authors themselves label the
         value as assumed (not measured): usable, reported as such
    [P]  no literature value available: plausible range, marked as pending
         (target of a future DFT / NEB / AIMD calculation by the group)
    [G]  group result (DFT/AIMD, Balbuena group), not to be changed here

Nothing tagged [M] or [P] is used in a production case until it is [V] or [G].
The MOBILITY.in pair interactions are [G] and are never modified.

## 1. How every rate is built (one formula, three traceable pieces)

    k_site = A * exp(-Ea / kT) * exp(alpha * F * eta / (R T))       [s^-1 per site]

    A     prefactor: kT/h = 6.2e12 s^-1 (transition-state theory) for chemical
          steps; for electrochemical steps A * exp(-Ea/kT) is replaced by the
          measured exchange current density converted to a per-site rate,
          k0_site = j0 * A_site / (n e), with A_site = (latticeConstant/2)^2
          = (2.0 A)^2 = 4.0e-16 cm^2 and e = 1.602e-19 C.
    Ea    barrier from DFT/NEB (group) or from the literature.
    eta   overpotential of the electrode during the half-cycle; alpha = 0.5
          unless a Tafel slope says otherwise.

Worked conversions (298 K, kT = 25.7 meV):
    j0 = 1 mA cm^-2, n = 1  ->  k0_site = 1e-3 * 4e-16 / 1.6e-19 = 2.5 s^-1
    j0 = 10 mA cm^-2, n = 1 ->  25 s^-1
    i0 = 1.97 A m^-2 (= 0.197 mA cm^-2), n = 2 -> 0.25 s^-1
    Ea = 0.364 eV with A = kT/h -> 6.2e12 * exp(-14.2) = 4.3e6 s^-1
    eta = 50 mV, alpha = 0.5 -> factor exp(0.97) = 2.6

The engine already evaluates rate = globalA * sigma * k0 * exp(-(Ea - alpha (V - E0))/kT)
per channel (DECOMPOSITION.in: sigma k0 Ea alpha E0). Mapping: k0 <- A or
k0_site, Ea <- Ea, alpha <- alpha, (V - E0) <- eta of the electrode. Needed
engine addition (opt-in): a per-electrode eta instead of the cell-voltage
label BeginV/EndV (docs/ANODE_KINETICS.md).

## 2. kMC time versus cycle time (a modelling decision, not a parameter)

With physical absolute rates (10 s^-1 per surface site, ~300 surface sites)
the engine fires ~3000 events per simulated second; a 1 h half-cycle would
need 1e7 events. That is out of reach, and it is why the original engine used
relative rates and an event budget (maxInterval). The physically based engine
keeps the SAME strategy but makes it explicit: rate RATIOS are physical
(from this table), the half-cycle is an EVENT budget chosen from the depth of
discharge (Li moved per half-cycle = maxInterval / 2 for a coupled cell), and
the simulated clock is reported but not identified with wall-clock cycling.
The C-rate enters through eta (larger |eta| = larger current), not through the
clock. This must be stated in every report.

## 3. Anode (Li metal in 1.2 M LiFSI / F5DEE)

Per-site conversion: k0_site = j0 * A_site / e, A_site = 4.0e-16 cm^2.

| Process | Parameter | Value @298 K | Per-site rate | Source (verified) | Tag |
|---|---|---|---|---|---|
| Li plating / stripping, LiFSI in DME (ether, weakly solvating) | j0 | 29.8 +- 0.7 mA cm^-2 | 74 s^-1 | Boyle et al., ACS Energy Lett. 2020, 5, 701, Table 1 (transient CV on ultramicroelectrodes, electron-transfer controlled) | [V] |
| Li plating / stripping, LiFSI in EC:DEC (carbonate) | j0 | 4.0 +- 1.1 mA cm^-2 | 10 s^-1 | Boyle 2020, Table 1 | [V] |
| Li plating / stripping in 1.2 M LiFSI / F5DEE | j0 | 30 mA cm^-2 adopted (F5DEE is an ether weaker-solvating than DME; Boyle 2020 attributes larger j0 to weaker Li-solvent binding, p. 706) | 74 s^-1 | interpolation from Boyle 2020; no direct measurement for F5DEE | [P, bounded 10 to 74 s^-1] |
| transfer coefficient | alpha | 0.5 at low overpotential; Boyle 2020 shows Butler-Volmer holds only for small eta and Marcus-Hush (lambda 0.30 to 0.34 eV, Table 1) fits the curved Tafel plots | | Boyle 2020, Fig. 2, Table 1 | [V] |
| SEI-limited effective j0 (EIS) | j0_EIS | about 2 orders of magnitude below the kinetic j0 | | Boyle 2020, p. 703 and Fig. S4: EIS is dominated by ion transport through the SEI | [V] |
| anode overpotential per half-cycle | eta | 20 to 50 mV at C/10 to C/5 (working choice; sets the current via BV) | factor 1.5 to 2.6 | to be taken from the group's Li||Cu polarization data | [P] |
| LiFSI reduction a-1..a-8, F5DEE reduction b-1..b-3 | Ea | ~0 (LICET, no outstanding barriers) | electron-transfer limited | Tan et al. JACS 2024, 146, 11711, p. 11714 | [V] |
| F5DEE on Li2O (defluorination, LiOH path), SOL2 | Ea | 0.364 eV | 0.25 exp(-0.364/kT) = 1.8e-7 s^-1 per adjacent site: the SOL prefactor (CE-anchored, effective electrode rate) with the DFT barrier as the RELATIVE attenuation of the Li2O route vs the barrierless Li route. The bare TST value kT/h exp(-Ea/kT) = 4.3e6 s^-1 is an elementary-step rate; used as a per-site kMC rate it made every O site an inexhaustible sink (run 1 of fullcell_physical). Upper bound for sensitivity: Ea 0 (1.8e-7 to 0.25 s^-1) | Tan 2024 Fig. 4; Perez-Beltran 2024; anchoring: Yu 2022 CE | [V-derived, CE-anchored] |
| electron tunnelling through the SEI | beta, d | exp(-beta d), beta ~1 A^-1 | | needed for barrierless reductions to slow down as the SEI thickens (Boyle's EIS versus CV gap is the experimental face of this) | [P] |
| Li surface diffusion | Eb pairs | MOBILITY.in | | Balbuena group DFT | [G] |
| CE (validation target) | CE | 99.5 to 99.9 % Li||Cu, 1.2 M LiFSI/F5DEE | | Yu et al. Nat. Energy 2022, 7, 94, Fig. 4 | [V] |

Consequences: the legacy stripping k0 (56 * 0.5 = 28 s^-1 per site) is inside the
physical window (10 to 74); the legacy plating k0 (1e-4) is five to six orders
too slow. Kumaresan 2008 assumes i0(Li) = 0.394 A m^-2 (0.04 mA cm^-2), two
to three orders below Boyle's measurement: Boyle's kinetic value is used.

## 4. Cathode (S8 -> Li2S, S8-unit lattice model)

Per-site conversion for i0 given per electron (n = 1): k0_site = i0 * A_site / e,
A_site = 4.0e-20 m^2.

| Process (0D/1D model step) | Parameter | Value @298 K | Per-site rate | Source (verified) | Tag |
|---|---|---|---|---|---|
| S8(l) + e -> 1/2 S8^2- (reaction 2) | i0 | 1.972 A m^-2, alpha 0.5, U 2.39 V | 0.49 s^-1 | Kumaresan, Mikhaylik, White, J. Electrochem. Soc. 2008, 155, A576, Table II ("assumed parameters", footnote a) | [V-assumed] |
| 3/2 S8^2- + e -> 2 S6^2- (reaction 3) | i0 | 0.019 A m^-2, U 2.37 V | 4.7e-3 s^-1 | Kumaresan 2008 Table II | [V-assumed] |
| S6^2- + e -> 3/2 S4^2- (reaction 4) | i0 | 0.019 A m^-2, U 2.24 V | 4.7e-3 s^-1 | Kumaresan 2008 Table II | [V-assumed] |
| 1/2 S4^2- + e -> S2^2- (reaction 5) | i0 | 1.97e-4 A m^-2, U 2.04 V | 4.9e-5 s^-1 | Kumaresan 2008 Table II | [V-assumed] |
| 1/2 S2^2- + e -> S^2- (reaction 6) | i0 | 1.97e-7 A m^-2, U 2.01 V | 4.9e-8 s^-1 | Kumaresan 2008 Table II | [V-assumed] |
| high plateau S8 -> S4^2- (4 e) | iH,0 | 10 A m^-2 | 2.5 s^-1 | Marinescu, Zhang, Offer, PCCP 2016, 18, 584, Table 1 | [V] |
| low plateau S4^2- -> S^2- (4 e) | iL,0 | 5 A m^-2 | 1.25 s^-1 | Marinescu 2016 Table 1 | [V] |
| Li2S precipitation (0D) | kp | 100 s^-1 (per unit supersaturation) | see note | Marinescu 2016 Table 1 | [V] |
| S8(s) dissolution | kS8 | 1.0 s^-1 (Kumaresan) ; 5.0 s^-1 (Zhang) | | Kumaresan 2008 Table V; Zhang, Marinescu, Walus, Offer, Electrochim. Acta 2016, 219, 502, Table B.4 | [V-assumed] |
| Li2S precipitation (1D) | k | 3.45e-5 m^6 mol^-2 s^-1, Ksp 1e2 mol^3 m^-9 (Zhang) ; 27.5 m^6 mol^-2 s^-1, Ksp 3.0e-5 (Kumaresan) | | Zhang 2016 Table B.4; Kumaresan 2008 Table V (the two differ by 6 orders; both assumed) | [V-assumed] |
| S8 solubility | Ksp | 19 mol m^-3 (DOL/DME-type electrolyte) | | Kumaresan 2008 Table V; Zhang 2016 Table B.4 | [V-assumed] |
| shuttle constant | ks | 0.53 h^-1 = 1.5e-4 s^-1 at 298 K (1.85 m LiTFSI); 0.45 / 0.14 / 0.095 h^-1 for 0.5 / 1.85 / 2.5 m salt; Ea(shuttle) 0.56 eV | k0 per Li0 surface site = ks * reservoir_sites / N_surf = 1.5e-4 * 18211 / 300 = 9e-3 s^-1 (10x10x30 box) | Mikhaylik and Akridge, J. Electrochem. Soc. 2004, 151, A1969, Figs. 9 and 11 and the thermal-model table (activation energy) | [V] |
| shuttle constant (0D model) | ks | 2e-4 s^-1 | 1.2e-2 s^-1 per site (same conversion) | Marinescu 2016 Table 1 | [V] |
| polysulfide solubility in fluorinated ethers | c_sat | qualitative only: suppression TTE > TFPG > TFEG, UV-vis | | Yue et al., J. Power Sources 2018, 401, 271, Sec. 3.3 | [V, qualitative] |
| FSI + long-chain polysulfide -> CEI | k | qualitative: pure-LiFSI cell collapses in ~8 cycles at C/10 | | Soria-Fernandez et al., Small Methods 2026, Fig. 3a | [V, qualitative] |
| Li2S / Li2S2 electronic passivation | nmin | geometric threshold | | no literature value for a lattice criterion | [P] |

Mapping to the S8-unit lattice steps (per-site k0, alpha 0.5, E0 = U):
    S8    -> Li2S8   : reaction 2, 0.49 s^-1, E0 2.39 V
    Li2S8 -> Li4S8   : reactions 3 + 4 in series, limiting 4.7e-3 s^-1, E0 2.24 V
    Li4S8 -> Li8S8   : reaction 5, 4.9e-5 s^-1, E0 2.04 V
    Li8S8 -> Li16S8  : reaction 6, 4.9e-8 s^-1, E0 2.01 V
The last two are so slow in Kumaresan's set that only the BV exponential moves
them (at eta = 0.1 V, factor exp(0.5 * 0.1 / 0.0257) = 7); Marinescu's lumped
iL,0 (1.25 s^-1 per site for the whole low plateau) is 4 orders faster. Both
are literature; the spread IS the uncertainty and must be reported. Working
choice: Marinescu's two-plateau values (iH,0 for the two high-plateau lattice
steps, iL,0 for the two low-plateau steps), because they were fitted to
cycling data, with Kumaresan's set as the sensitivity bound.

Shuttle: the adopted 0.05 s^-1 per site (calibration step 1) is 4 to 5x the
Mikhaylik / Marinescu conversion (0.009 to 0.012); the legacy 0.5 was 40 to
50x. Next value: 0.01 s^-1 [V-derived].

## 5. PDFs verified (all six in hand, 2026-09-22)

1. Boyle et al., ACS Energy Lett. 2020, 5, 701.  2. Kumaresan, Mikhaylik, White,
J. Electrochem. Soc. 2008, 155, A576.  3. Marinescu, Zhang, Offer, PCCP 2016,
18, 584.  4. Mikhaylik and Akridge, J. Electrochem. Soc. 2004, 151, A1969.
5. Yue et al., J. Power Sources 2018, 401, 271.  6. Zhang, Marinescu, Walus,
Offer, Electrochim. Acta 2016, 219, 502.
Still open: j0 measured in F5DEE itself; polysulfide solubility in F5DEE;
tunnelling attenuation through the SEI.

## 6. Second stage (group calculations)

Once the literature values fix the orders of magnitude, the group's own
DFT/NEB/AIMD can replace [M] and [P] entries with system-specific ones:
polysulfide reduction barriers on the carbon/S8 surface, FSI attack by
Li2S6/Li2S8 (barrier and products), F5DEE reactivity toward polysulfides,
Li2Sx solubility in F5DEE (free energy of solvation). Each replacement is one
line in this table and one line in DECOMPOSITION.in, with the same tag logic.

## 7. The SEI acceleration factor (decision to state in every report)

With physical ratios, side reactions are (1 - CE)/CE ~ 2e-3 of the plating
events. In a 40 x 40 A patch cycled with 50 events per half-cycle that is
about one decomposition event every 10 cycles: the SEI that forms is a few
sites, which is physically right (experimental SEI growth per cycle is
sub-nanometre after formation) but useless for morphology studies. The
original engine obtained a visible SEI by setting k_decomposition ~ k_plating
(1e-4 : 1e-4), i.e. an implicit acceleration of ~500 on the side reactions.
anode_physical makes this explicit: the CE-anchored rates are the reference,
and any morphology study multiplies FSI / SFO / SOL k0 by a stated
acceleration factor A_SEI (e.g. 100 or 500), reported alongside the results.
Plating, stripping and the cathode kinetics are never accelerated, so the
electrode coupling stays physical.

## 8. Salt-site artifact removed in bv mode (2026-09-22)

In the legacy scheme plating is a candidate only where a salt (FSI) site
touches Li0. Salt is 75 randomly redistributed sites, so a few percent of the
time no plating candidate exists and the BKL fires the only remaining channel,
a solvent/salt decomposition, with a 100 to 600 s time jump. The apparent CE
was therefore a geometric property of the salt representation, not of the
kinetics (anode_physical iteration 1: 3.8 % side events at k_side/k_plate =
2e-3). In anode_kinetics bv the plating weight is k_bv times the number of
growth sites (Li+ is available throughout a 1.2 M electrolyte) and plating
never runs out of candidates; decompositions still require an adjacent salt or
solvent site. With this, the side-event fraction follows the k0 ratio as
intended.

## 9. Cathode implementation (2026-09-23)

cathode_kinetics bv evaluates every REGION cathode channel as
rate = globalA sigma k0 exp(-(Ea - alpha (V_cat - E0)) / kT) with V_cat the
half-cycle cathode potential (cathode_V_charge 2.45 V, cathode_V_discharge
1.9 V, tag [P] until replaced by the reference cell's voltage profile), E0 the
Kumaresan 2008 standard potentials and alpha signed: -0.5 reduction, +0.5
oxidation. Because E0 decreases along the cascade (2.39, 2.24, 2.04, 2.01 V),
the plateaus emerge from the BV factors: at 1.9 V the S8 step is driven by
exp(0.5 * 0.49 / 0.0257) ~ 1.4e4, the last step by exp(0.5 * 0.11 / 0.0257)
~ 8.5; at 2.45 V the order reverses. The absolute factors are large, but in a
BKL only their ratios and the pool coupling (li_pool_max) decide the event
split; the anode side keeps its own physical scale (74 s^-1 per site).
Dissolution / precipitation keep the milestone-1 placeholders (0.05 / 0.02
and 5.0 s^-1; Kumaresan's kS8 = 1.0 s^-1 is the only literature anchor and
refers to S8(s), which the S8-unit lattice does not dissolve as such).
