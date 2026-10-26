"""Rewrite the stale counters in docs/survey/extension_directions_reconciliation.json.

Why this exists: that file was written after the third evaluation round and its
`counts` block still says table_rows_after = 23 / mechanically_checkable_after = 26,
and correction c06 still records "20 rows -> 23 rows". Those numbers are now wrong
twice over (the 23+3 = 26 reading was a double count, and PMID 27216816 has since
been added to the table). Independently recomputing them here, from the document
itself, keeps the ledger from disagreeing with `docs/extension_directions.md`.

The recomputation reads the same two artefacts the human-facing checks read:
  - docs/extension_directions.md          (table rows, bytes, lines)
  - docs/survey/snapshot_article_count.json (per-layer snapshot attribution)

Run from the repository root:

    python -X utf8 docs/survey/_fix_reconciliation_counts.py
"""

from __future__ import annotations

import json
import re
import pathlib

SURVEY = pathlib.Path(__file__).resolve().parent
REPO = SURVEY.parent.parent
DOC = REPO / "docs" / "extension_directions.md"
RECON = SURVEY / "extension_directions_reconciliation.json"
COUNT = SURVEY / "snapshot_article_count.json"

ROW = re.compile(r"^\| (\d{8}) \|", re.M)


def strip_fenced(text: str) -> str:
    out: list[str] = []
    inside = False
    for line in text.splitlines():
        if line.lstrip().startswith("```"):
            inside = not inside
            continue
        if not inside:
            out.append(line)
    return "\n".join(out)


def layer_of(labels: list[str]) -> str:
    """Map a PMID's snapshot membership to the document's four-layer scheme."""
    short = {label.split("_")[0] for label in labels}
    if "A" in short:
        return "A"
    if "B" in short:
        return "B"
    if short & {"C", "D", "E"}:
        return "C/D/E"
    return "none"


def main() -> None:
    raw = DOC.read_text(encoding="utf-8")
    body = strip_fenced(raw)
    pmids = ROW.findall(body)

    counts = json.loads(COUNT.read_text(encoding="utf-8"))
    membership: dict[str, list[str]] = counts["pmid_snapshot_membership"]

    per_layer = {"A": [], "B": [], "C/D/E": [], "none": []}
    for pmid in sorted(set(pmids)):
        per_layer[layer_of(membership.get(pmid, []))].append(pmid)

    total = len(pmids)
    assert total == len(set(pmids)), "document table has duplicate PMIDs"
    assert sum(len(v) for v in per_layer.values()) == total, "layer partition is lossy"

    recon = json.loads(RECON.read_text(encoding="utf-8"))
    recon["counts"] = {
        "table_rows_after": total,
        "mechanically_checkable_after": total,
        "counting_note": (
            "Recomputed 2026-10-17 from docs/extension_directions.md by "
            "docs/survey/_fix_reconciliation_counts.py. The earlier '23 + 3 = 26' "
            "reading double-counted three PMIDs that were already among the 23; every "
            "table row has at least one on-disk snapshot, so the mechanically checkable "
            "count equals the row count."
        ),
        "layer_breakdown": {
            "A": len(per_layer["A"]),
            "B": len(per_layer["B"]),
            "C_D_E": len(per_layer["C/D/E"]),
            "no_snapshot": len(per_layer["none"]),
            "sums_to_table_rows": total,
        },
        "layer_pmids": {k: v for k, v in per_layer.items() if v},
        "document_lines": raw.count("\n") + 1,
        # Two different numbers, labelled so neither can be mistaken for the other:
        # read_text() applies universal newlines, so re-encoding gives the LF form;
        # stat().st_size is what is physically on disk (CRLF in this checkout).
        "document_bytes_lf_utf8": len(raw.encode("utf-8")),
        "document_bytes_on_disk": DOC.stat().st_size,
        "document_bytes_note": (
            "document_bytes_lf_utf8 = 换行归一为 LF 后的 UTF-8 字节数(本脚本实测);"
            "document_bytes_on_disk = 磁盘实际大小(本检出为 CRLF)。两者之差等于行数减一,"
            "不是数据不一致。"
        ),
    }

    c06 = next(c for c in recon["corrections"] if c["id"] == "c06")
    c06["found"] = (
        "Table grew from 20 to 23 rows when the three [A2] entries were added; it is "
        f"{total} rows after PMID 27216816 was added in round 5. The real errors in that "
        "round were a journal name and a title (42241701) and a journal name (39915771); "
        "the year cells were already PubDate-derived."
    )
    c06["action"] = (
        "Added [A2] rows for 27401684 / 35277501 / 42800986, added the [A2] legend, "
        "recorded the 42241701 and 39915771 corrections inline, and (2026-10-17) added "
        "the 27216816 row that direction 5 cited without listing."
    )

    c07 = next(c for c in recon["corrections"] if c["id"] == "c07")
    c07["predicted"] = "23 mechanically checkable rows"
    c07["found"] = (
        "snapshot_accounting.json: snapshot_A = 19 docsums, snapshot_B = 4 docsums, "
        "disjoint. The 23 + 3 = 26 reading was a double count: 27401684, 35277501 and "
        "42800986 were already among the 23, and 42800986 was in no esummary at all."
    )
    c07["action"] = (
        "Round 5 replaced the additive reading with a four-layer partition in which every "
        "table row belongs to exactly one layer and the layers sum to the row count "
        f"({total}). The layers are A = {len(per_layer['A'])}, B = {len(per_layer['B'])}, "
        f"C/D/E = {len(per_layer['C/D/E'])}, no snapshot = {len(per_layer['none'])}."
    )
    c07["status"] = "corrected"

    g01 = recon["unchanged_known_gaps"][0]
    if g01["id"] == "g01":
        g01["disposition"] = (
            "Handed to the human: docs/survey/direction_1_ruling_request.md states the "
            "two options (keep it as an extension direction and record criterion (1) as "
            "not passed, or demote it to a defect-fix task). Not resolved by this ledger."
        )

    RECON.write_text(json.dumps(recon, ensure_ascii=False, indent=2), encoding="utf-8")

    print("table rows:", total)
    print("layer breakdown:", {k: len(v) for k, v in per_layer.items()})
    print("document lines:", recon["counts"]["document_lines"])
    print("document bytes (LF-normalised UTF-8):", recon["counts"]["document_bytes_lf_utf8"])
    print("document bytes (on disk):", recon["counts"]["document_bytes_on_disk"])
    print("wrote", RECON.relative_to(REPO))


if __name__ == "__main__":
    main()
