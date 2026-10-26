"""Apply the post-audit corrections to extension_directions_reconciliation.json.

Re-runnable; writes the JSON back in place with the corrected c04 record and a
post_audit_notes section. Run from the repository root:

    python docs/survey/_fix_reconciliation.py
"""

from __future__ import annotations

import json
from pathlib import Path

PATH = Path(__file__).resolve().parent / "extension_directions_reconciliation.json"

CORRECTED_C04 = {
    "id": "c04",
    "target": "docs/extension_directions.md section direction 4",
    "predicted": "no content-adjacent literature exists (first-round conclusion)",
    "found": (
        "Second round with different roots returned PMID 27401684 "
        "(Tools and techniques for computational reproducibility, GigaScience, "
        "PubDate 2016 Jul 11, DOI 10.1186/s13742-016-0135-4). "
        "Note: 10.1093/gigascience/giz095 belongs to a different paper (31675414); "
        "an earlier draft of this record mis-paired the two. "
        "Auxiliary: 35277501 (FAIRly big), 42800986 (FairFlow, efetch-only)."
    ),
    "action": (
        "Rewrote the direction-4 literature block: direct-adjacent evidence 27401684 "
        "plus two auxiliary; kept 31675414 and 39915771 explicitly as non-supporting "
        "side notes, with the reason each fails isomorphism."
    ),
    "status": "corrected",
}

NOTE = {
    "id": "n01",
    "target": "c04 / docs/survey/citation_ledger.md R1",
    "issue": (
        "27401684 was cited as 2017 by reading efetch <DateCompleted> 2017-10-30 as the "
        "publication year; its PubDate and ArticleDate are both 2016-07-11."
    ),
    "status": "corrected",
    "evidence": "docs/survey/date_field_audit.json -> by_pmid.27401684.efetch",
}


def main() -> None:
    data = json.loads(PATH.read_text(encoding="utf-8"))
    kept = []
    replaced = False
    for item in data["corrections"]:
        if item.get("id") == "c04":
            kept.append(CORRECTED_C04)
            replaced = True
        else:
            kept.append(item)
    if not replaced:
        raise SystemExit("c04 not found; refusing to write")
    data["corrections"] = kept

    notes = data.setdefault("post_audit_notes", [])
    if not any(n.get("id") == "n01" for n in notes):
        notes.append(NOTE)

    PATH.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding="utf-8")
    reloaded = json.loads(PATH.read_text(encoding="utf-8"))
    print("corrections:", len(reloaded["corrections"]))
    print("post_audit_notes:", [n["id"] for n in reloaded["post_audit_notes"]])
    print("c04 found:", reloaded["corrections"][3]["found"][:160])


if __name__ == "__main__":
    main()
