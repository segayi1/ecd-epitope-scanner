# Fixed Jurkat-Fab x K562-antigen assay profile

## Purpose

The upload workflow is specialized for coarse epitope localization in a two-cell assay. Jurkat
cells display the antibody Fab receptor; K562 cells display WT or mutant antigen ECD. The user
may submit an antigen-ECD PDB/mmCIF without a UniProt accession or a separate sequence.

## Frozen construct interpretation

The interpretation is tied to the exact SnapGene file hashes stored in `assay_profiles.py`.

- **K562 antigen:** signal peptide - Myc - antigen ECD - PDGFRB stalk/TM - intracellular
  PDGFRB/FUS-7A. The ECD N terminus is extracellular and its C terminus is membrane-proximal.
- **Jurkat light chain:** signal peptide - FLAG - VL insert - mouse Ig-kappa constant - ALZ.
- **Jurkat heavy receptor:** signal peptide - Myc - VH insert - mouse IgG1 CH1 - BLZ - CD8A
  stalk - Notch1 TM - Notch2 JMD - Gal4 DBD - VP64. The Fab is extracellular and signaling is
  read on the Jurkat side.

These roles are experimental facts supplied by the user. Feature annotations in the plasmid
files are evidence for architecture, not instructions to the program.

## PDB-only behavior

1. If no chain is specified, the longest observed protein chain is selected deterministically.
2. Its coordinate-derived amino-acid sequence becomes the ECD sequence. Numbering is therefore
   an observed-sequence index; each residue record retains the original PDB chain/residue ID.
3. The last three observed C-terminal CA atoms define a robust anchor point. The vector from
   this point to the ECD centroid is reported as a **relative extracellular axis**.
4. This axis is not called a membrane normal. A soluble ECD structure cannot determine the
   physical membrane plane, receptor tilt, intermembrane spacing or Fab reach.

For multichain complexes, homologous oligomers or structures containing non-antigen protein
chains, the chain override should be used; “longest chain” is only a deterministic fallback.

## Probe policy

- Enumerate the structure-derived surface regions before applying the experimental budget.
- Select at most six Stage-1 regions by greedy maximin separation of region centroids. Epitope
  likelihood priors do not choose the panel.
- Design 2-4 substitutions per region from surface-exposed, folding-screened sites.
- Do not penalize C-anchor-proximal regions. This is a tethered Fab versus tethered antigen
  assay, and PDB-only geometry is insufficient for a quantitative reach penalty.
- Report unselected regions and all protected/blind-spot residues. Six constructs are a budget,
  not a claim of complete physical coverage.

## Required experimental interpretation gates

Reporter loss is classified as a candidate regional epitope effect only when K562 antigen
surface display is retained relative to WT and Jurkat Fab display/assembly is acceptable. Myc
accessibility is an expression/display control, not proof of native antigen folding. An
independent conformation-sensitive antibody or ligand assay is required to distinguish regional
binding disruption from local/global misfolding. WT, no-antigen, irrelevant-antigen and
signaling-positive controls should be run in the same batch.
