# Version 1.1.0 validation

Date: 2026-11-14. Environment: Windows, Python 3.12, existing project dependencies.

- The installed candidate wheel passed 189 tests in 141.21 seconds, including integration fixtures; one Starlette deprecation warning was emitted. This run preceded the epitope-prior provenance correction.
- The final source passed four focused tests of the provenance flag, including missing prior, selector identity, eligible candidates and partial target-surface evidence.
- Task-center JavaScript tests: 2 passed. Browser submission: the optional perturbation protocol and a 12-construct maximum were submitted successfully without page errors.
- A wheel rebuilt from the source distribution had identical archive members. The source distribution and wheel passed recursive privacy scans for supported local paths and configured private markers.

The tests use saved inputs and controlled model fixtures. No new external-model calculations, fresh-machine installation or wet-lab experiments were performed for this release. Field and epitope-prior values are computational proxies; scientific and experimental generalization remains unverified.
