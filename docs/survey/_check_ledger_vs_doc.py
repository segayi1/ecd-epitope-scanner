"""Cross-check the citation ledger against docs/extension_directions.md section 5.

Run from the repository root:
    python -X utf8 docs/survey/_check_ledger_vs_doc.py

Counts are taken from the FENCED-BLOCK-FREE text only: the CLI examples embedded in
the markdown contain the same `| NNNNNNNN |` pattern and would otherwise be counted as
table rows (this produced a spurious "ledger 25 rows" reading on 2026-10-17).

Checks performed, all mechanical:
  1. set equality -- every PMID cited in the document has a ledger row, and reports
     the ledger-only leftovers (snapshot hits the document never cited);
  2. year -- the document's year cell vs the ledger's CURRENT-VALUE column;
  3. journal -- the document's journal cell vs the journal in the落盘 snapshot,
     matched after abbreviation-tolerant normalisation (the document abbreviates:
     "Nucleic Acids Res" vs snapshot "Nucleic acids research");
  4. title fragments -- every quoted English phrase (two or more words) in the
     document's usage cell must appear verbatim in the snapshot title. This is the
     check that catches "the document says X but the title says Y".

Item 3 and 4 were added on 2026-10-17 because the previous version compared years
only and left title/journal consistency to manual reading.
"""
from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LEDGER = ROOT / "docs" / "survey" / "citation_ledger.md"
DOC = ROOT / "docs" / "extension_directions.md"
OUT = ROOT / "docs" / "survey" / "ledger_vs_doc_check.json"

# Snapshot files that together cover every PMID cited in the document.
SNAPSHOTS = [
    ROOT / "docs" / "survey" / "lit_resnapshot_20261016.xml",
    ROOT / "docs" / "survey" / "lit_resnapshot_gap_20261016.xml",
    ROOT / "docs" / "survey" / "_tmp_retry_esummary2.xml",
    ROOT / "docs" / "survey" / "_tmp_retry_efetch.xml",
    ROOT / "docs" / "survey" / "_fetch_esummary_42800986.xml",
]

ROW = re.compile(r"^\| (\d{8}) \|", re.M)
QUOTED_ASCII = re.compile(r"`([A-Za-z][A-Za-z0-9 ,:;()\-]{7,})`|\"([A-Za-z][A-Za-z0-9 ,:;()\-]{7,})\"")

# 19055730 is the single documented calendar-year exception (online-first 2008-12-02,
# journal issue year 2009). 25048123's snapshot journal is the book-series long form,
# which the document shortens to "Methods Mol Biol".
ACCEPTED_YEAR_EXCEPTIONS = {"19055730": "快照 PubDate 2008-12-02,本文按期刊卷期年取 2009"}
ACCEPTED_JOURNAL_NOTES = {
    "25048123": "快照 FullJournalName 为 Methods in molecular biology (Clifton, N.J.),本文写作 Methods Mol Biol",
}


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


def section(text: str, start: str, end: str) -> str:
    i = text.index(start)
    j = text.index(end, i + len(start))
    return text[i:j]


def row_cells(block: str) -> dict[str, list[str]]:
    """Map PMID -> the pipe-separated cells of its ledger row (header/sep rows skipped)."""
    rows: dict[str, list[str]] = {}
    for line in block.splitlines():
        if not line.startswith("| "):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 4 or not re.fullmatch(r"\d{8}", cells[0]):
            continue
        rows[cells[0]] = cells
    return rows


def doc_rows(block: str) -> dict[str, list[str]]:
    """Map PMID -> cells of the document's section 5 row."""
    rows: dict[str, list[str]] = {}
    for line in block.splitlines():
        if not line.startswith("| "):
            continue
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 3 or not re.fullmatch(r"\d{8}", cells[0]):
            continue
        rows[cells[0]] = cells
    return rows


def lead_year(cell: str) -> str | None:
    """First 4-digit year in a cell, ignoring a bold marker. None when absent."""
    m = re.search(r"(\d{4})", cell.replace("*", ""))
    return m.group(1) if m else None


def norm_journal(value: str) -> str:
    """Lower-case, drop punctuation and 'the', collapse spaces."""
    text = value.lower().replace("&", "and")
    text = re.sub(r"\(.*?\)", " ", text)
    text = re.sub(r"[^a-z0-9 ]", " ", text)
    words = [w for w in text.split() if w and w != "the"]
    return " ".join(words)


def journal_compatible(doc_journal: str, snapshot_journal: str) -> bool:
    """Accept NLM-style abbreviation: each doc word must prefix a snapshot word, in order.

    "Brief Bioinform" vs "Briefings in bioinformatics" is the case that a positional
    word-by-word comparison gets wrong ("bioinform" would be aligned against "in").
    Walking the snapshot words in order and allowing skips matches real abbreviations
    while still failing on a genuinely different journal.
    """
    doc = norm_journal(doc_journal).split()
    snap = norm_journal(snapshot_journal).split()
    if not doc or not snap:
        return False
    if " ".join(doc) == " ".join(snap):
        return True
    cursor = 0
    for word in doc:
        while cursor < len(snap):
            candidate = snap[cursor]
            cursor += 1
            if candidate.startswith(word) or word.startswith(candidate):
                break
        else:
            return False
    return True


def load_snapshots() -> dict[str, dict[str, str]]:
    """PMID -> {title, journal, year, doi, snapshot_file} from every落盘快照."""
    found: dict[str, dict[str, str]] = {}
    for path in SNAPSHOTS:
        if not path.exists():
            continue
        root = ET.parse(path).getroot()
        for docsum in root.findall(".//DocSum"):
            pmid = (docsum.findtext("Id") or "").strip()
            if not re.fullmatch(r"\d{8}", pmid):
                continue
            items = {item.get("Name"): (item.text or "") for item in docsum.findall("Item")}
            found[pmid] = {
                "title": items.get("Title", ""),
                "journal": items.get("FullJournalName", "") or items.get("Source", ""),
                "year": (items.get("PubDate", "") or "")[:4],
                "doi": items.get("DOI", ""),
                "snapshot_file": path.name,
            }
        for art in root.findall(".//PubmedArticle"):
            pmid = (art.findtext(".//MedlineCitation/PMID") or "").strip()
            if not re.fullmatch(r"\d{8}", pmid):
                continue
            journal = art.findtext(".//Journal/Title") or art.findtext(".//Journal/ISOAbbreviation") or ""
            pubdate = art.find(".//Journal/JournalIssue/PubDate")
            year = ""
            if pubdate is not None:
                year = (pubdate.findtext("Year") or "")[:4] or (pubdate.findtext("MedlineDate") or "")[:4]
            doi = ""
            idlist = art.find("PubmedData/ArticleIdList")
            if idlist is not None:
                for node in idlist.findall("ArticleId"):
                    if node.get("IdType") == "doi":
                        doi = (node.text or "").strip()
            found.setdefault(
                pmid,
                {
                    "title": art.findtext(".//ArticleTitle") or "",
                    "journal": journal,
                    "year": year,
                    "doi": doi,
                    "snapshot_file": path.name,
                },
            )
    return found


def load_text_snapshot() -> dict[str, str]:
    """PMID -> Title from docs/survey/lit_resnapshot_20261016.txt (readable form).

    The text snapshot is the artifact a human reads; comparing it against the XML it
    was rendered from catches a re-typed or hand-edited title line.
    """
    path = ROOT / "docs" / "survey" / "lit_resnapshot_20261016.txt"
    titles: dict[str, str] = {}
    if not path.exists():
        return titles
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = [p.strip() for p in line.split("|")]
        if len(parts) >= 5 and re.fullmatch(r"\d{8}", parts[0]):
            titles[parts[0]] = parts[4]
    return titles


def main() -> None:
    ledger_raw = LEDGER.read_text(encoding="utf-8")
    doc_raw = DOC.read_text(encoding="utf-8")

    ledger_body = strip_fenced(ledger_raw)
    doc_body = strip_fenced(doc_raw)
    # Section 1 sits between the plain-text reference to "## 2." and the real heading.
    sec1 = section(ledger_body, "## 1.", "\n## 2. 本轮修正记录")

    ledger_rows = ROW.findall(sec1)
    doc_all = ROW.findall(doc_body)
    # Doc section 5 is the only table whose first column is an 8-digit PMID.
    ledger_set = set(ledger_rows)
    doc_set = set(doc_all)

    ledger_cells = row_cells(sec1)
    doc_cell_map = doc_rows(doc_body)
    snapshots = load_snapshots()

    year_table: dict[str, dict[str, str | None]] = {}
    year_mismatches: list[dict[str, str]] = []
    for pmid in sorted(doc_set & ledger_set):
        cells = ledger_cells.get(pmid)
        doc_cells = doc_cell_map.get(pmid)
        if cells is None or doc_cells is None:
            continue
        doc_year = lead_year(doc_cells[1])
        ledger_current = lead_year(cells[3]) if len(cells) > 3 else None
        year_table[pmid] = {
            "doc": doc_year,
            "ledger_current": ledger_current,
            "ledger_superseded": lead_year(cells[4]) if len(cells) > 4 else None,
        }
        if doc_year and ledger_current != doc_year:
            year_mismatches.append(
                {"pmid": pmid, "doc": doc_year, "ledger_current": ledger_current or ""}
            )

    journal_table: dict[str, dict[str, str]] = {}
    journal_mismatches: list[dict[str, str]] = []
    title_fragment_table: dict[str, dict[str, list[str]]] = {}
    title_fragment_missing: list[dict[str, str]] = []
    no_snapshot: list[str] = []

    for pmid in sorted(doc_set):
        cells = doc_cell_map[pmid]
        doc_journal = cells[2] if len(cells) > 2 else ""
        usage = cells[3] if len(cells) > 3 else ""
        snap = snapshots.get(pmid)
        if snap is None or not snap.get("title"):
            no_snapshot.append(pmid)
            continue

        journal_table[pmid] = {"doc": doc_journal, "snapshot": snap["journal"]}
        if not journal_compatible(doc_journal, snap["journal"]):
            journal_mismatches.append(
                {"pmid": pmid, "doc": doc_journal, "snapshot": snap["journal"]}
            )

        fragments = [a or b for a, b in QUOTED_ASCII.findall(usage)]
        # Keep only prose-like fragments (contain a space) - bare identifiers such as
        # `benchmark_landscape.json` are file references, not title quotes.
        prose = [f for f in fragments if " " in f.strip()]
        present = [f for f in prose if f.lower() in snap["title"].lower()]
        missing = [f for f in prose if f.lower() not in snap["title"].lower()]
        if prose:
            title_fragment_table[pmid] = {"quoted": prose, "found_in_title": present}
        for fragment in missing:
            title_fragment_missing.append(
                {"pmid": pmid, "fragment": fragment, "snapshot_title": snap["title"]}
            )

    txt_titles = load_text_snapshot()
    txt_vs_xml_mismatches: list[dict[str, str]] = []
    for pmid, txt_title in sorted(txt_titles.items()):
        snap = snapshots.get(pmid)
        if snap is None or not snap.get("title"):
            continue
        if txt_title.strip().lower() != snap["title"].strip().lower():
            txt_vs_xml_mismatches.append(
                {
                    "pmid": pmid,
                    "text_snapshot": txt_title,
                    "xml_snapshot": snap["title"],
                }
            )

    report = {
        "generated_by": "docs/survey/_check_ledger_vs_doc.py",
        "note": "fenced code blocks are stripped before counting",
        "ledger_section1_rows": len(ledger_rows),
        "ledger_section1_unique": len(ledger_set),
        "doc_pmid_rows": len(doc_all),
        "doc_unique": len(doc_set),
        "cited_but_not_in_ledger": sorted(doc_set - ledger_set),
        "in_ledger_but_not_cited": sorted(ledger_set - doc_set),
        "canonical_pmid_count": len(doc_set),
        "year_checked_count": len(year_table),
        "year_exceptions": ACCEPTED_YEAR_EXCEPTIONS,
        "year_mismatches": year_mismatches,
        "journal_checked_count": len(journal_table),
        "journal_notes": ACCEPTED_JOURNAL_NOTES,
        "journal_mismatches": journal_mismatches,
        "title_fragment_checked_count": len(title_fragment_table),
        "title_fragment_missing": title_fragment_missing,
        "text_vs_xml_title_checked_count": len(txt_titles),
        "text_vs_xml_title_mismatches": txt_vs_xml_mismatches,
        "cited_pmids_without_snapshot": no_snapshot,
        "year_table": year_table,
        "journal_table": journal_table,
        "title_fragment_table": title_fragment_table,
        "pass": not (doc_set - ledger_set)
        and not year_mismatches
        and not journal_mismatches
        and not title_fragment_missing
        and not txt_vs_xml_mismatches
        and not no_snapshot,
    }
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    for key, value in report.items():
        if key in {"year_table", "journal_table", "title_fragment_table"}:
            continue
        print(f"{key}: {value}")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
