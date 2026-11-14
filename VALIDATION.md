# Computational validation

Validation date: 2026-11-09. The local application environment used Windows and Python 3.12.10. SPURS and DiscoTope ran in existing WSL/CUDA environments; paired energy calculations used licensed FoldX 5.1.

## External setup audit

The external setup guide pins the upstream revisions and direct dependency versions observed in the recorded SPURS and DiscoTope environments. The four SPURS model/configuration files match the recorded inputs; 28 DiscoTope source/model-archive files match the pinned upstream revision. Offline checks cover shell syntax, worker-command resolution, and inspection of the recorded strict output. The installation commands have not been executed on a fresh machine. See [the setup audit](acceptance/external_setup_audit.json).

## Version 1.0.1

The installed distribution passed 91 local tests. Export redaction covers supported Windows, home, WSL and POSIX temporary-directory paths. Regression cases exercise nested request metadata, compressed artifacts, public URLs and application routes. The service and exported software metadata identify the installed patch version. Current Python 3.11/3.12 results are available in the repository's GitHub Actions runs.

## Installed distribution and interface

The version 1.0.0 baseline wheel passed **84 tests in 110.12 seconds** from a directory outside the source checkout. A separate extraction of the source archive rebuilt a wheel with identical payload contents.

The suite covers structure/sequence mapping, confidence provenance, AF3 atom mapping, connected surface geometry, site-first scheduling, construct/sampling budgets, evidence pairing, missing evidence, export consistency, privacy, and packaged browser resources. Automated scheduling tests use controlled scoring fixtures.

An Edge browser preview completed to 100%. Reloading its saved-job URL restored the structure viewer, 220 mapped residues, and 11 regions. Formal recommendations were empty in the absence of external stability evidence. No application console errors were observed during this check.

## Strict PD-L1 analysis

The input was the public PD-L1 fixture supplied in this repository: Q9NZQ7, chain A, inclusive ECD residues 19–238. The PDB was processed as an ordinary structure upload with no inferred pLDDT assignment from its B-factors.

| Measure | Recorded result |
| --- | --- |
| Protocol | ordinary_site_first_v2 |
| Surface regions | 11 |
| Exported designs | 57 |
| Stage 1 probes | 7 |
| Stage 2 single / pairwise designs | 16 / 8 |
| Regional designs | 26 |
| Principle queues | 33 complete; no unfinished gaps |
| SPURS combination requests | 62; 0 failed |
| FoldX paired slots | 113 valid; 0 failed; maximum budget 300 |
| Released FoldX reservations | 217, accumulated across reservations |
| Portable FoldX evidence | 33 records; all 113 pairs passed validate_pair |
| Evidence cards | 66; PASS; no missing fields |
| Global probe panel | 7 |
| Result inventory | 736 listed files verified within a 737-file ZIP |

The full result archive contained 544,203,356 bytes, with SHA-256:
`dbff36f77f069ca80cba83f3c5c1293e0e61e2a700cf0282dc45ae3a7e6c31e3`.

The executed-source snapshot had SHA-256:
`952f69ba3f66e498f8488f16ba47ab54686f4947ef67cfb8ec942fa358a27498`.

Comparison of 124 executed-source files with the release implementation identified differences in public API restrictions, web entry points, and assay metadata. Pipeline, model-worker, scheduler, scientific-gate, and export implementations were unchanged. The required SPURS inference modifications are provided as a patch against the documented upstream revision and were verified against the tested local implementation.

## Evidence availability

[acceptance/README.md](acceptance/README.md) indexes the published summaries, protocol, test log, browser record, and evidence-card receipt. Additional sanitized queue, design, and inventory records accompany the GitHub release. Full historical run archives and executed-source snapshots are retained locally; their identifiers in the summaries support provenance matching.

Publication files were scanned for known private identifiers and user-directory paths, including nested archive contents. Public assay metadata contains display-topology assumptions without private construct filenames, digests, or vector coordinates.

## Scope

ProCaliper was not configured or invoked in the strict analysis; its tool-inventory entry indicated installation only. The recorded real-model case covers the public ordinary-PDB workflow. AF3 confidence handling is covered by regression tests, while real AF3 multi-model inference and fresh-machine external installation were not part of this acceptance run. Python 3.11/Linux execution is tracked separately by repository CI.

The quantitative results above describe one computational case. They do not constitute measurements of experimental binding, expression, or folding, and do not prescribe candidate counts or runtime for other proteins.

## Version 1.1.0 candidate validation

The optional `ordinary_site_first_perturbation_v1` protocol adds sector-first scheduling, perturbation-field ranking and per-construct audit output. The default remains `ordinary_site_first_v2`. A formal result requires the configured model evidence and scientific gates; field values and DiscoTope scores are computational proxies, not experimental binding measurements.

The version 1.1.0 candidate passed 189 installed-wheel tests on Windows/Python 3.12, including integration fixtures and mock model scheduling. Four additional tests passed against the final source after the epitope-prior provenance correction. The browser submission check confirmed the optional protocol and 12-construct budget. The source archive rebuilt a wheel with identical member contents; source and wheel passed the bundled recursive privacy checker with local marker checks. Details are in acceptance/version_1_1_0.md. Validation used existing dependencies and did not repeat external-model calculations or prove generalization. A development-version CTLA4 run using the new protocol selected 12 formally qualified constructs; this computational result does not establish binding effects or wet-lab performance.
