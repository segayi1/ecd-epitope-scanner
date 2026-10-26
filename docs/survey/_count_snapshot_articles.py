"""Structured, line-count-free audit of the four on-disk PubMed snapshots.

Why this exists: a reviewer counted `<PubmedArticle>` occurrences by physical
line in `_tmp_retry_efetch.xml` and got 1, because PubMed returns each record
unwrapped onto very few physical lines. Counting XML elements is the only
honest count. This script prints and persists that count so the question never
has to be settled by recollection.

Run from the repository root:

    python -X utf8 docs/survey/_count_snapshot_articles.py
"""

from __future__ import annotations

import json
import pathlib
import xml.etree.ElementTree as ET

SURVEY = pathlib.Path(__file__).resolve().parent
REPO = SURVEY.parent.parent

# label -> file, plus how to read it
SOURCES = {
    "A_docsum": SURVEY / "lit_resnapshot_20261016.xml",
    "B_docsum_gap": SURVEY / "lit_resnapshot_gap_20261016.xml",
    "C_docsum_retry": SURVEY / "_tmp_retry_esummary2.xml",
    "D_efetch_retry": SURVEY / "_tmp_retry_efetch.xml",
    # E was pulled on 2026-10-17 to give PMID 42800986 a second channel; it is an
    # esummary view of one record that otherwise exists only in D.
    "E_esummary_42800986": SURVEY / "_fetch_esummary_42800986.xml",
}

DOC = REPO / "docs" / "extension_directions.md"


def _docsum_pmids(root: ET.Element) -> list[str]:
    """DocSum records carry the PMID in a bare <Id> child (not an <Item>)."""
    out: list[str] = []
    for docsum in root.findall("DocSum"):
        pmid = docsum.findtext("Id")
        if pmid:
            out.append(pmid.strip())
            continue
        # some esummary variants expose it as an Item instead
        for item in docsum.findall("Item"):
            if item.get("Name") == "PMID" and item.text:
                out.append(item.text.strip())
                break
    return out


def _pubmedarticle_pmids(root: ET.Element) -> list[str]:
    out: list[str] = []
    for art in root.findall("PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID")
        if pmid:
            out.append(pmid.strip())
    return out


def read_one(path: pathlib.Path) -> dict:
    with path.open(encoding="utf-8", newline="") as fh:
        raw = fh.read()
    root = ET.fromstring(raw)
    if root.tag == "PubmedArticleSet":
        kind = "PubmedArticleSet"
        ids = _pubmedarticle_pmids(root)
        # count both ways so a line-based reading cannot be mistaken for this one
        opens = raw.count("<PubmedArticle>")
        closes = raw.count("</PubmedArticle>")
    else:
        kind = root.tag
        ids = _docsum_pmids(root)
        opens = raw.count("<DocSum>")
        closes = raw.count("</DocSum>")
    return {
        "file": path.name,
        "root_tag": kind,
        "child_elements": len(list(root)),
        "pmids": ids,
        "pmid_count": len(ids),
        "open_tag_occurrences": opens,
        "close_tag_occurrences": closes,
        "chars": len(raw),
        "physical_line_count": len(raw.splitlines()),
    }


def main() -> None:
    per_file = {}
    for label, path in SOURCES.items():
        if not path.exists():
            per_file[label] = {"file": path.name, "error": "missing"}
            continue
        per_file[label] = read_one(path)

    where: dict[str, list[str]] = {}
    for label, info in per_file.items():
        for pmid in info.get("pmids", []):
            where.setdefault(pmid, []).append(label)

    cited = sorted(set(__import__("re").findall(r"\b(\d{8})\b", DOC.read_text(encoding="utf-8"))))
    cited_in_snapshots = {p: where.get(p, []) for p in cited if p in where}
    cited_without_snapshot = sorted(p for p in cited if p not in where)

    efetch = per_file.get("D_efetch_retry", {})
    cited_not_in_snapshots = sorted(p for p in where if p not in cited)
    cited_in_snapshots_only = sorted(p for p in where if p in cited)
    result = {
        "generated_by": "docs/survey/_count_snapshot_articles.py",
        "method": "xml.etree.ElementTree parsing; never physical-line counting",
        "per_file": per_file,
        "cross_file_counts": {label.split("_")[0]: info.get("pmid_count", 0)
                              for label, info in per_file.items()},
        "distinct_pmids": len(where),
        "pmid_snapshot_membership": where,
        "distinct_pmids_breakdown": {
            "note": (
                "distinct_pmids above is every PMID that appears in ANY snapshot; it is NOT "
                "the number of cited references. The two buckets below partition it exactly "
                "and must sum to distinct_pmids."
            ),
            "cited_in_document": len(cited_in_snapshots_only),
            "in_snapshot_but_not_cited": len(cited_not_in_snapshots),
            "cited_pmid_list": cited_in_snapshots_only,
            "not_cited_pmid_list": cited_not_in_snapshots,
            "sums_to_distinct_pmids": (
                len(cited_in_snapshots_only) + len(cited_not_in_snapshots) == len(where)
            ),
        },
        "efetch_pubmedarticle_count": efetch.get("pmid_count"),
        "efetch_pmids": efetch.get("pmids"),
        "efetch_physical_line_count": efetch.get("physical_line_count"),
        "line_vs_element_note": (
            "The efetch return has very few physical lines because PubMed does not wrap "
            "records; scanning line by line reports far fewer records than parsing XML "
            "elements. Element count is authoritative."
        ),
        "cited_pmids_total": len(cited),
        "cited_with_snapshot": cited_in_snapshots,
        "cited_without_snapshot": cited_without_snapshot,
    }

    out = SURVEY / "snapshot_article_count.json"
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")

    print("per-file element counts:")
    for label, info in per_file.items():
        print(
            f"  {label}: {info.get('pmid_count')} elements "
            f"({info.get('child_elements')} child nodes), "
            f"{info.get('chars')} chars, {info.get('physical_line_count')} physical lines"
        )
    print("efetch elements:", result["efetch_pubmedarticle_count"], result["efetch_pmids"])
    print("distinct PMIDs across snapshots:", result["distinct_pmids"])
    print("cited PMIDs in doc:", result["cited_pmids_total"],
          "| without any snapshot:", result["cited_without_snapshot"])
    print("wrote", out.relative_to(REPO))


if __name__ == "__main__":
    main()
