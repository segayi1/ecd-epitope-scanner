# PD-L1 fixture

`structure.pdb` is the public AlphaFold Protein Structure Database model for
human PD-L1 (UniProt Q9NZQ7), chain A. AlphaFold DB structures are distributed
under CC BY 4.0:

<https://alphafold.ebi.ac.uk/entry/Q9NZQ7>

`canonical.fasta` contains the canonical 290-residue UniProt sequence. The
fixture uses ECD residues 19-238 and `allow_partial=True` because the
AlphaFold model contains the full-length protein.

The integration test runs in `preview` mode with external engines disabled:

```powershell
python -m pytest --import-mode=importlib tests/test_pipeline_fixture.py -q
```

This exercises structure loading, residue mapping, surface partitioning,
preview candidate generation, missing-evidence handling, exports and the
offline report. Preview uses the legacy execution path; it does not validate
the site-first V2 scheduler, SPURS, DiscoTope, ProCaliper or FoldX.
