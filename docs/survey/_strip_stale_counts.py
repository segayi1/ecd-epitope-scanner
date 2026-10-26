"""Strip superseded literal row counts out of the reconciliation ledger's prose.

The numbers themselves now live in `counts`, recomputed by
`_fix_reconciliation_counts.py`. Leaving a stale "23" / "23 + 3 = 26" inside a
prose field of c06/c07 is what made two evaluation rounds read the ledger as
disagreeing with the document, so the prose points at `counts` instead.

Idempotent: running it twice is a no-op. Run from the repository root:

    python -X utf8 docs/survey/_strip_stale_counts.py
"""

from __future__ import annotations

import json
import pathlib

SURVEY = pathlib.Path(__file__).resolve().parent
REPO = SURVEY.parent.parent
RECON = SURVEY / "extension_directions_reconciliation.json"


def main() -> None:
    recon = json.loads(RECON.read_text(encoding="utf-8"))
    by_id = {c["id"]: c for c in recon["corrections"]}

    by_id["c06"]["found"] = (
        "Table grew by three rows when the three [A2] entries were added; it stands at "
        "24 rows after PMID 27216816 was added in round 5. The real errors in that round "
        "were a journal name and a title (42241701) and a journal name (39915771); the "
        "year cells were already PubDate-derived. (The row count this field used to state "
        "is superseded; the current count is in counts.table_rows_after.)"
    )

    by_id["c07"]["predicted"] = (
        "The table rows would all be mechanically checkable, with the count equal to the "
        "row count."
    )
    by_id["c07"]["found"] = (
        "snapshot_accounting.json: snapshot_A = 19 docsums, snapshot_B = 4 docsums, "
        "disjoint. The earlier additive reading summed a subset with three PMIDs that were "
        "already inside it, and one of those three (42800986) had no esummary at all; it "
        "is therefore a double count, not a layer. See counts.counting_note."
    )

    RECON.write_text(json.dumps(recon, ensure_ascii=False, indent=2), encoding="utf-8")

    text = RECON.read_text(encoding="utf-8")
    print("counts.table_rows_after:", recon["counts"]["table_rows_after"])
    print("counts.mechanically_checkable_after:", recon["counts"]["mechanically_checkable_after"])
    print("stale literal present:", "23 + 3 = 26" in text or "3 = 26" in text)
    print("wrote", RECON.relative_to(REPO))


if __name__ == "__main__":
    main()
