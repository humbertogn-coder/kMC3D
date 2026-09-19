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

- Soria-Fernandez, A.; Castillo, J.; Cid, R.; Song, Z.; Wu, H.; Carriazo, D.;
  Armand, M.; Zhang, H.; Santiago, A. Beyond the Hype: Decoding
  Bis(fluorosulfonyl)imide Chemistry in Advanced Lithium-Sulfur Batteries.
  Small Methods 2026, DOI 10.1002/smtd.202502084. VERIFIED (abstract):
  raising the LiFSI fraction (LiFSI/LiTFSI localized high-concentration
  electrolytes) improves conductivity and Li-metal compatibility but "reduces
  sulfur utilization through side reactions with long-chain polysulfides";
  0.2 M LiFSI co-salt found optimal. This is the key qualitative statement:
  FSI- is consumed by long-chain Li2Sx (x = 6, 8), i.e. a chemical (not
  electrochemical) CEI-forming path that also removes active sulfur.
  TO-OBTAIN full text: product identification (LiF, Li2SO4/LiSO2F,
  thiosulfate/polythionate, S-N species), whether the loss scales with
  [Li2Sx] and [FSI], any rate or extent numbers.
- Reference list of that review (Crossref) includes, among others:
  10.1021/ja412807w, 10.1039/C4CC06666A, 10.1021/jp408037e,
  10.1002/adfm.201505074, 10.1021/acsenergylett.6b00194,
  10.1021/acsami.1c09492, 10.1021/acsami.3c10977, 10.1021/acsami.3c14048,
  10.1016/j.ensm.2024.103501, 10.1002/aenm.202302378. TO-RESOLVE titles.

## Fluorinated ethers vs polysulfides

- Zu, C.; Manthiram, A. Insight into lithium-metal anodes in lithium-sulfur
  batteries with a fluorinated ether electrolyte. J. Mater. Chem. A 2015,
  DOI 10.1039/c5ta03195h. METADATA (TTE co-solvent). TO-OBTAIN.
- Okuda, D. et al. Mechanism of polysulfide dissolution suppression in
  lithium-sulfur batteries using a novel fluorinated ether electrolyte.
  Electrochim. Acta 2026, DOI 10.1016/j.electacta.2026.149154. METADATA.
  TO-OBTAIN.
- Needed: evidence on whether polysulfide nucleophiles attack the C-F or
  C-O bonds of fluorinated DEE solvents (F5DEE has CF3CH2O- and HCF2CH2O-
  termini). Expectation from the anode network (F5DEE defluorination is
  electron-transfer driven, negligible barrier) is that F5DEE is a minor
  chemical CEI contributor compared with FSI at 2 to 3 V; to be confirmed.

## What is NOT in this electrolyte

- No DOL/DME: ring-opening polymerization of DOL by polysulfides and DME
  deprotonation do not apply.
- No LiNO3: the LiNxOy / thiosulfate / polythionate CEI of Aurbach et al.
  (J. Electrochem. Soc. 2009, 156, A694) does not apply unless nitrate is
  added later (keep as an optional future block).

## Proposed reaction set (to be finalized once the full texts are read)

    cei_FSI_Li2S8   REGION cathode_surface, TRIGGER FSI, REQUIRE Li2S8_d 1,
                    RATE_SCALE Li2S8_d: CEI CEI_F 0 ; RES Li2S8_d -1 ;
                    RES Li2S6_d +1 ; S_LOSS 2 (placeholder split: FSI + Li2S8 ->
                    LiF/sulfoxy film + Li2S6, 2 S into the film)
    cei_FSI_Li2S6   same pattern with Li2S6_d -> Li2S4_d
    cei_F5DEE_Li2Sx optional, TRIGGER SOL, much slower (weight to be justified)

Rates: none assigned yet. Constraint from the Small Methods review: the
FSI side reaction must be strong enough to lower S utilization at high LiFSI
fraction, weak enough that 0.2 M LiFSI co-salt is still beneficial.
