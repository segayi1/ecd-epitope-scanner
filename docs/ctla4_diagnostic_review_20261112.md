# CTLA4 diagnostic review — 2026-11-12

## Result

Job `1d812aa2f9be` completed using UniProt P16410 and experimental structure 1AH1 chain A. The ECD mapping covers 125 residues; the pipeline reported 76 surface-accessible residues, eight display regions and 12 selected constructs. All 12 selected combinations have SPURS multi and FoldX results recorded as `CONSENSUS_PASS`. This is computational screening, not experimental validation of folding or antibody binding.

The saved panel reports 89.10% of eligible surface area within its 12 Å geodesic support. This is geometric support, not a measured perturbation or binding-loss rate. Its synthetic hypothesis signatures leave 861 near-indistinguishable pairs. All 12 selected constructs are triples. The original job used development behavior under the `ordinary_site_first_v2` name and is not a version-isolated acceptance run.

The trace reached export, the job store records COMPLETE, all 557 manifest hashes verified, and the saved card-evidence completeness check passed. The ProCaliper, SPURS, DiscoTope and FoldX workers ran without a recorded worker failure; SPURS multi and FoldX worker identity metadata remain incomplete. The input is an experimental NMR structure, so AF3 confidence evidence does not apply. Engineered positions 162–165 were omitted from scoring and design.

## Findings and repairs

- Conservation failed with `KeyError` when an alignment column contained no identical residues. The calculation now records a conservation score of zero for that column; the historical result remains unchanged.
- The original panel ranked triple-to-parent margin before the marginal panel field, selecting two additional triples whose physical and weighted field gains were both zero. The selector now ranks field gains before that margin after sector priority; a regression covers the observed conflict. Positive geometric coverage remains eligible when no field gain is available.
- The V4 partition summary was hard-coded to report that epitope prior was not used. The updated export derives this flag from selector provenance and scored evidence among eligible candidates. The original analysis and fixed assay profile are historical records and were not rewritten.
- The 196 MB panel export and 190 MB analysis file are principally caused by repeated provenance objects, particularly the adaptive ledger and model metadata. This remains unresolved; the existing job-detail endpoint was slow during inspection.
- The FoldX diagnostics incorrectly report zero runtime: the adaptive engine reconstructs its result without carrying over the worker timings. The saved ledger separately records 309.31 seconds of BuildModel and 33.65 seconds of repair; the diagnostics bug remains open.
- Four regions have induced-region path diameters greater than 30 Å, while the configured hard limit uses shortest paths in the full eligible surface graph. The result satisfies the declared full-graph diameter convention; it must not be represented as satisfying a 30 Å bound on paths restricted to each region.

Focused regression after these changes: 97 passed, 2 skipped, with one Starlette deprecation warning. No post-fix real-model job or historical shadow replay has been performed.

## Remaining acceptance work

The surface field is a candidate-level chemical strength propagated uniformly over a 12 Å geometric support, with mutation-site epitope prior propagated across that support. It is a coarse ranking proxy, not a spatially resolved physical measurement or residue-specific antibody effect. Triple utility sums nominal chemical changes and has not been calibrated against experimental binding changes. Version isolation for the new selector, actual target-point epitope weighting, provenance deduplication, full replay and validation of the revised panel remain open. The SNIPR data were not used as production scores. FoldX/SPURS gates are predictive and do not prove correct experimental folding.
