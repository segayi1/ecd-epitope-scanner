# AF3 preprocessing final state — 2026-10-22

Status: implementation and computational acceptance COMPLETE for the approved scope. No active work is required for this task. No worker/scoring/deployment was launched; existing dirty tree and historical outputs preserved.

Authoritative delivery: test_output/af3_final_delivery_20261022/evaluation_zh.md and acceptance.json. Final fresh (NO geometry reuse) comparison: test_output/af3_preprocessing_comparison_20261022_32de0d8a/report_zh.md, comparison.json, comparison.csv; fingerprints bind current implementation and each original input.

Final identities: AIB10 model1(m001); AIB203/AIB26/AIB56/AIB60 model0(m000). Only1/5changes. Prior claim2/5 and AIB26m002 was a mistaken pTM/ranking tie-baseline reading; do not reuse earlier summaries. All25models and100%canonicalmapping/coverage validated, originalsourcefiles unchanged.

Implementation: af3_preprocessing.py discovery/freeze/identity-validation/source groups/fingerprints/quality/coremedoid/near-tie/fallback/diagnostics; PipelineRequest.af3_manifest_path; --af3-directory; create/resume /api/jobs/.../af3-directory; originalselectedcoordinates/confidence; online/offline metrics, archives and replay. Common flow uses selected model exactly once. Diagnostics do not add mutation exclusion rules; workers/formal gates/80-300contract unchanged. Final comparison script cache removed to avoid attaching current source fingerprints to old geometric outputs.

Verification: focused92PASS/1Windows symlink permission SKIP(test_output/af3_final_focused_20261022.xml), real complete folder-vs-selected-pair1PASS(af3_final_equivalence_v4_20261022.xml). Mapping,reliability,surface,displayregions,AAcandidates equal; exported bundles PASS. Existing stochastic CDSoptimizer frozen only in test. Core per-residue statistics.mean matches legacy adapter exactly. PD1core replay v2 maps100%ECD,288confidencevalues,PAEmapped. Capacity2048tokens/2048target,bothfull2048²matricesPASS; exactmemoryin capacity_receipt.json, not100max-sizedmodeltest.

Scientific result: AIB10coreRMSD.252Å/global.616Å, pTM/rankingunchanged, meanpLDDT+.0166, contactsJaccard.967, eligible67to66, SurfaceV2regions3to4, regionARI.204. Main benefit is traceability and uncertainty diagnostics; no experiment-based accuracy gain demonstrated. No new SPURS/FoldX or validated mutation recommendation comparison. True multimodelECD/externalexperimentalbenchmarks remain future work.

Browser file:// action rejected by URLpolicy, no workaround; UI staticJS and offlineexportchecks passed. Live service was not restarted. No automation or heartbeat. No need for further polling.
