# ECD epitope-region mutation design report

## Interpretation boundary

This workflow designs antigen mutants for experimental testing. It does not infer an antibody-specific epitope and it does not prove that a mutant is correctly folded. A construct supports epitope attribution only when antibody binding decreases while expression, folding controls and appropriate native-function controls remain acceptable.

## Target and structure

- UniProt: USER_STRUCTURE
- Protein: User-provided structure
- ECD: 19-238 (220 aa)
- Topology: user-defined ECD boundary
- Structure: structure.pdb
- Chain: A
- Observed ECD coverage: 100.0%
- Maximum internal coordinate gap: 0

## Surface and folding-risk screen

- SASA method: Bio.PDB Shrake-Rupley fallback 1.4A/960 points
- Surface criterion: rSASA >= 0.25 and side-chain SASA >= 20.0 A2
- Protected features: disulfide cysteines, annotated glycosylation/PTMs, missing mappings and low-confidence AF3 residues
- Native interfaces participate in design but retain a function-risk flag

## Partition and design

- Partition method: domain- and membrane-orientation-constrained geodesic surface partitioning
- Regions: 6
- Eligible-surface coverage: 100.0%
- In-silico substitutions evaluated: 912
- Final constructs: 6
- Stability interpretation: EvoEF2 energy units when the external engine is configured; otherwise a dimensionless within-target proxy

## Required experimental gates

1. Confirm total and surface expression relative to WT.
2. Confirm a non-competing conformational antibody or another domain-specific folding control remains bound.
3. Measure monodispersity/thermal stability when using purified ECD.
4. Measure native ligand/function when functional-interface mutations are interpreted.
5. Treat loss of all antibodies or controls as a folding/expression failure, not an epitope.
