# Physical anode kinetics for the full cell (proposal, 2026-09-22)

Status: DESIGN, no code yet. To be reviewed with the advisor before implementation.

## Why

The legacy anode kinetics (Perez-Beltran 2024) are relative rates tuned so that
SEI growth and Li morphology came out right over ~150 half-cycles of ~50 events:
plating k0 1e-4 per salt site, stripping k0 56 per surface Li0 site,
decomposition 1e-4 to 1e-6, all with alpha = 0 (no potential dependence). They
reproduce the original engine byte for byte and stay untouched in
cases/fullcell_shuttle and validation/.

In a full cell they break the electrode coupling (series 2 and the N/P test,
2026-09-22): stripping is fast, so the cathode discharges fully, but plating is
~5 events per half-cycle, so during charge the pool fills to li_pool_max and
the cathode cannot oxidize (RELEASE_LI gated). The cell discharges but does not
charge. The asymmetry is unphysical: plating and stripping are the same
reaction in opposite directions and must share one exchange current density.

## Formulation (Butler-Volmer per lattice site)

Rate per eligible site, s^-1:

    k_plate = k0_BV * exp(-alpha_c * F * eta / (R T))
    k_strip = k0_BV * exp(+alpha_a * F * eta / (R T))
    k0_BV   = j0 * A_site / (n F)          [s^-1 per site]

with j0 the exchange current density (A cm^-2), A_site = a0^2 the area one
lattice site exposes to the electrolyte (a0 = latticeConstant / 2 = 2.0 A ->
4e-16 cm^2), n = 1, alpha_a = alpha_c = 0.5 (to be justified), eta the anode
overpotential (V) set by the half-cycle: eta < 0 during charge (plating),
eta > 0 during discharge (stripping). A constant |eta| per half-cycle
(galvanostatic-like) is the simplest defensible choice; a current-following eta
(solve eta so that the anode current equals the cathode consumption) is the
next step and follows naturally from the shared pool.

Side reactions (FSI, F5DEE decomposition): keep the DFT barriers of the anode
network (Ea from the Tan 2024 / Perez-Beltran 2024 network; F5DEE on Li2O
0.364 eV) with an electrochemical prefactor kT/h and the same alpha * eta
term, so that CE ~ 1 - (side reaction rate / plating rate) is an OUTPUT, not
an input. Target: CE 0.95 to 0.99 for LiFSI/F5DEE (Yu 2022 reports 99.5 to
99.9 % in Li||Cu half cells).

## Engine support (mostly present)

The DECOMPOSITION.in channel (sigma, k0, Ea, alpha, E0) already evaluates
rate = globalA * sigma * k0 * exp(-(Ea - alpha (V - E0)) / kT). Two additions
are needed: (1) a per-electrode potential, so that V - E0 is the anode
overpotential rather than the cell voltage label (BeginV/EndV), behind an
opt-in keyword (e.g. anode_eta_charge / anode_eta_discharge in PARAMETERS.in);
(2) the stripping rate read from the same BV block as plating instead of
MOBILITY.in LiStripping. Everything else is input values.

## Inputs to collect (literature, with citation per value)

- j0 for Li plating/stripping in LiFSI / fluorinated ether electrolytes at
  298 K (candidates: Yu et al. Nat. Energy 2022 exchange-current or Tafel
  data; Boyle et al. ACS Energy Lett. 2020 on transient voltammetry of Li
  deposition, j0 of order 1 to 10 mA cm^-2 in ether electrolytes; group data).
- alpha (0.5 default, or fitted Tafel slopes from the same sources).
- eta per half-cycle from the C-rate and the cell's measured overpotential
  (tens of mV at C/10 for Li metal).
- Ea of the decomposition network: already in the group's DFT.

## Validation plan

1. anode_physical case (half cell, 10x10x30, 150 half-cycles): CE, SEI
   thickness and Li growth against the original run and against experiment
   (CE 0.995 to 0.999). Not byte-identity: a new physics, a new reference.
2. Full cell with the same BV block: charge and discharge must be symmetric in
   Li+ moved (pool returns to li_pool0 each cycle) and the cathode must
   recharge.
