# Project continuity

Before changing the pipeline, surface geometry, model workers, reporting or UI,
read `docs/project_memory.md` and the relevant repair audit linked there.

- Preserve existing user changes and historical run directories. Do not stage
  the entire working tree when unrelated edits are present.
- Reproduce defects using saved inputs with their structure type, residue
  mapping and confidence provenance. Replay AF3 as AF3, with its confidence file.
- Display regions must be connected and satisfy configured diameter limits.
  Requested region count is a soft target; construct budget is independent.
- Keep geometry verification separate from validated mutation recommendations.
  Never fabricate missing evidence or relax scientific gates to produce output.
- Run package tests from the repository's parent directory, as documented in
  `docs/project_memory.md`. Verify meaningful invariants and export consistency.
- When fixing a recurrence, update the project memory with cause, affected code,
  regression evidence and remaining limitations.
