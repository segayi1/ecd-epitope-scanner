# Recorded computational acceptance

These summaries document the public PD-L1 example, Q9NZQ7 chain A / ECD 19–238, using `ordinary_site_first_v2`. They are computational acceptance evidence, not experimental binding data.

- [Strict acceptance summary](strict_acceptance.json): design/region counts, queue completion, 113 validated FoldX pairs and the full result archive's digest.
- [Execution protocol](execution_protocol.json): actual scoring budgets and policies.
- [Evidence-card receipt](card_evidence_completeness.json): displayed fields and recorded input digests. The complete referenced result files are not included here.
- [Installed-package test log](pytest.txt): 84 passed on Windows / Python 3.12.
- [Browser check](browser_acceptance.json): completed preview and saved-result reload.
- [SPURS patch check](spurs_patch_verification.json): applicability to the pinned upstream code and agreement with the tested implementation.
- [Executed-source comparison](source_comparison.json): comparison against the locally executed source snapshot before publication documentation was added.

Additional sanitized queue/design/inventory evidence is available in the acceptance ZIP attached to [v1.0.0](https://github.com/segayi1/ecd-epitope-scanner/releases/tag/v1.0.0). Original run archives and their executed-source snapshots are not public attachments. A digest in a summary identifies the corresponding local archive; it does not make that archive available for independent inspection here.

ProCaliper was not configured/invoked in the strict run. The PDB fixture was treated as an ordinary upload, not a verified AF3 confidence input. See [VALIDATION.md](../VALIDATION.md) for the full scope and limits.
