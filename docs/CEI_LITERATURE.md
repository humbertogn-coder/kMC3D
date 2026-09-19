# CEI in Li-S with LiFSI / F5DEE: literature basis (work in progress)

Purpose: every CEI reaction and rate placed in cases/*/MECHANISM.in and
DECOMPOSITION.in must point to an entry here. Status tags: VERIFIED (abstract
or full text read), METADATA (title/DOI confirmed, content not yet read),
TO-OBTAIN (needed, not yet available offline).

## Electrolyte definition (from the group's own work)

- 1.2 M LiFSI in F5DEE, 1-(2,2-difluoroethoxy)-2-(2,2,2-trifluoroethoxy)ethane,
  fluorinated-DEE family (Bao group, Battery500). VERIFIED in Perez-Beltran,
  Kuai, Balbuena, ACS Energy Lett. 2024, 9, 5268 (DOI 10.1021/acsenergylett.4c02019)
  and in the group's 2026 draft. Anode reduction network: LiFSI a-1..a-8
  (S-N cleavage, [O=S=O]2-, [O=S(=O)(F)N]2-, then F-, O2-, S2-), F5DEE b-1..b-4
  (defluorination, vinyl and alkoxide anions, organic amorphous phase);
  F5DEE on Li2O: Ea 0.364 eV, LiOH path.
- Yu, Z. et al. Rational solvent molecule tuning for high-performance lithium
  metal battery electrolytes. Nat. Energy 2022, 7, 94-106. METADATA (cited as
  ref 13 in the 2024 paper). TO-OBTAIN full text: does it report Li-S cells
  or polysulfide solubility in F-DEE solvents?
- Tan, S.; Kuai, D.; Yu, Z.; Perez-Beltran, S.; et al. Evolution and Interplay
  of Lithium Metal Interphase Components Revealed by Experimental and
  Theoretical Studies. JACS 2024, 146, 11711. METADATA. TO-OBTAIN.

## FSI anion vs polysulfides (the CEI-relevant chemistry)

- [SF26] Soria-Fernandez, A.; Castillo, J.; Cid, R.; Song, Z.; Wu, H.;
  Carriazo, D.; Armand, M.; Zhang, H.; Santiago, A. Beyond the Hype: Decoding
  Bis(fluorosulfonyl)imide Chemistry in Advanced Lithium-Sulfur Batteries.
  Small Methods 2026, DOI 10.1002/smtd.202502084. VERIFIED (full text).
  System: LiFSI/LiTFSI localized high-concentration electrolytes in
  sulfolane/TTE, E/S 7 uL/mg, C/10. Findings used here:
  * "FSI- anions undergo reductive decomposition upon interacting with LiPS,
    consuming both salt and active material while producing electronically
    insulating compounds that passivate the cathode surface"; driver is "the
    lability of the S-F bonds". Chemical path, no electrode electron transfer
    needed (VOLTAGE any in MECHANISM.in).
  * Products on the cycled cathode (XPS S 2p, Fig. 3c-d): sulfite/thiosulfate
    (SO3, 167 eV) and sulfate, plus LiF (F 1s, Fig. S13). Atomic % on the
    pure-LiFSI cathode: SO4 ~15, SO3 ~44, LixSy ~13, Li2Sx ~16; on LiTFSI:
    ~10, ~5, ~29, ~37. Long-chain signal depleted where the film grows.
  * Chain-length selectivity: with SPAN cathodes (only Li2S2 to Li2S4) the
    LiFSI cell is stable, so short-chain polysulfides do not drive the
    reaction. Hence cei reactions only for Li2S8_d and Li2S6_d.
  * Solvation dependence: in sparingly solvating, anion-rich environments the
    anion LUMO is lowered and FSI- is "more susceptible to reduction by
    long-chain polysulfides".
  * Magnitude (Fig. 3a): pure LiFSI cell drops from ~1050 to ~150 mAh/g
    within ~8 cycles at C/10; 0.2 M LiFSI co-salt still beneficial.
    This is the calibration target for k0 of cei_FSI_Li2Sx (TODO).
  * Anode side: sulfite/thiosulfate also found on the Li metal (Fig. S15),
    more with LiFSI: a possible future extension of the shuttle branch.

## Solvation structure of 1.2 M LiFSI / F5DEE (why [SF26] applies here)

- [Yu22] Yu, Z. et al. Rational solvent molecule tuning for high-performance
  lithium metal battery electrolytes. Nat. Energy 2022, 7, 94-106. VERIFIED
  (full text). 1.2 M LiFSI/F5DEE Li+ solvates: 7.5 % SSL + 11.9 % LASP, the
  rest (~80 %) Li-anion clusters (LAC, >= 2 anions). Weakly solvating; Li+
  binds the -CHF2 F (1.96 A) more than -CF3 (2.04 A). No Li-S data in the
  paper. Ref. 38 therein (Yue et al., J. Power Sources 2018, 401, 271) covers
  partially fluorinated ethers for Li-S electrolytes: TO-OBTAIN.
- [Tan24] Tan, S.; Kuai, D.; Yu, Z.; Perez-Beltran, S.; et al. Evolution and
  Interplay of Lithium Metal Interphase Components. JACS 2024, 146, 11711.
  VERIFIED (full text). LiFSI/F5DEE on Li metal: solvation shell is mostly
  FSI-Li-F5DEE contact ion pairs; S-F cleavage has thermodynamic priority
  over N-S when the anion is not solvent-coordinated (relevant to the S-F
  lability invoked by [SF26]); every LICET step of LiFSI is more favourable
  than that of F5DEE, i.e. LiFSI has faster decomposition kinetics under
  competitive conditions; F5DEE spin density on the -CHF2 side; F5DEE on
  Li2O: 0.364 eV barrier. No polysulfide data. Used only to rank the solvent
  path below the anion path (cei_F5DEE_Li2S8 disabled, sigma 0).

## Fluorinated ethers vs polysulfides (still open)

- Zu, C.; Manthiram, A. J. Mater. Chem. A 2015, DOI 10.1039/c5ta03195h (TTE
  co-solvent). METADATA. TO-OBTAIN.
- Okuda, D. et al. Electrochim. Acta 2026, DOI 10.1016/j.electacta.2026.149154.
  METADATA. TO-OBTAIN.
- Yue, Z. et al. J. Power Sources 2018, 401, 271. METADATA. TO-OBTAIN.
- Open question: do polysulfide dianions attack the C-F / C-O bonds of F5DEE
  (dehydrofluorination of the -CH2CHF2 end)? Until evidence exists the solvent
  CEI path stays disabled.

## What is NOT in this electrolyte

- No DOL/DME: ring-opening polymerization of DOL by polysulfides and DME
  deprotonation do not apply.
- No LiNO3: the LiNxOy / thiosulfate CEI of Aurbach et al. (J. Electrochem.
  Soc. 2009, 156, A694) does not apply unless nitrate is added later.

## Implemented reaction set (cases/fullcell_cei, 2026-09-19)

    cei_FSI_Li2S8   REGION cathode_surface, VOLTAGE any, TRIGGER FSI,
                    REQUIRE Li2S8_d 1, RATE_SCALE Li2S8_d
                    channel 0 (w 1.0): CEI CEI_SOx ; RES Li2S8_d -1 ; S_LOSS 8 ; LI_LOSS 2
                    channel 1 (w 0.0): CEI CEI_SOx ; RES Li2S8_d -1 ; RES Li2S6_d +1 ; S_LOSS 2
    cei_FSI_Li2S6   same with Li2S6_d (-> Li2S4_d in channel 1)
    cei_F5DEE_Li2S8 TRIGGER SOL, product CEI_org, DISABLED (sigma 0)

    Rates (DECOMPOSITION.in): k0 1.0 per FSI surface site, mean-field in
    c_x. PLACEHOLDER: same footing as the shuttle (0.5 per Li0 surface site);
    with ~20 FSI vs ~300 Li0 surface sites the CEI takes ~1/10 of the initial
    long-chain consumption. Calibration target: [SF26] Fig. 3a (collapse in
    ~8 cycles at C/10 for a pure-LiFSI cell) once the S8-unit cascade exists.
    Passivation: cathode_passivating_species CEI_SOx,CEI_org, nmin 10 of the
    18 to 22 electrolyte-facing neighbours (geometric placeholder).

    Smoke test (6x6x24, k0 boosted to 50, 12 half-cycles): CEI events fire,
    FSI surface salt depletes 18 -> 9, film 26 sites, li_total constant 416,
    s_total constant apart from the known cascade leak, no crash on restart
    path. Not a physical run.
