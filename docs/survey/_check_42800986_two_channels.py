"""Cross-check PMID 42800986 between the two independent PubMed channels.

Run from the repository root:
    python -X utf8 docs/survey/_check_42800986_two_channels.py

Channel D: docs/survey/_tmp_retry_efetch.xml  (efetch, PubmedArticle record)
Channel E: docs/survey/_fetch_esummary_42800986.xml (esummary, DocSum record, fetched 2026-10-17)

The 42800986 entry was previously supported by the efetch channel only, which made
"verified by a second independent pull" untrue for that citation. This script is the
mechanical check behind whatever the ledger and the directions document now claim.
"""
from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SURVEY = ROOT / "docs" / "survey"
EFETCH = SURVEY / "_tmp_retry_efetch.xml"
ESUMMARY = SURVEY / "_fetch_esummary_42800986.xml"
OUT = SURVEY / "_esummary_42800986_check.json"

PMID = "42800986"


def efetch_record(pmid: str) -> dict | None:
    root = ET.parse(EFETCH).getroot()
    for art in root.findall(".//PubmedArticle"):
        got = art.findtext(".//PMID")
        if got != pmid:
            continue
        pd = art.find(".//Journal/JournalIssue/PubDate")
        year = pd.findtext("Year") if pd is not None else None
        medline_date = pd.findtext("MedlineDate") if pd is not None else None
        # Scope to this record's OWN identifier list. A bare ".//ArticleId" also matches
        # every ArticleId inside ReferenceList, whose last entry used to be picked up as
        # the article DOI (42800986 -> 10.1038/nature21406, a reference's DOI).
        doi = None
        own = art.find("PubmedData/ArticleIdList")
        if own is not None:
            for aid in own.findall("ArticleId"):
                if aid.get("IdType") == "doi":
                    doi = (aid.text or "").strip()
        return {
            "channel": "efetch",
            "file": EFETCH.name,
            "pmid": got,
            "title": (art.findtext(".//ArticleTitle") or "").strip().rstrip("."),
            "journal": art.findtext(".//Journal/Title"),
            "journal_iso": art.findtext(".//Journal/ISOAbbreviation"),
            "year": year,
            "medline_date": medline_date,
            "volume": art.findtext(".//JournalIssue/Volume"),
            "doi": doi,
        }
    return None


def esummary_record(pmid: str) -> dict | None:
    root = ET.parse(ESUMMARY).getroot()
    for doc in root.findall("DocSum"):
        if doc.findtext("Id") != pmid:
            continue
        items: dict[str, str] = {}
        for item in doc.findall("Item"):
            name = item.get("Name")
            if name:
                items[name] = (item.text or "").strip()
        return {
            "channel": "esummary",
            "file": ESUMMARY.name,
            "pmid": doc.findtext("Id"),
            "title": items.get("Title", "").rstrip("."),
            "journal": items.get("FullJournalName"),
            "journal_iso": items.get("Source"),
            "year": (items.get("PubDate") or "").split(" ")[0] or None,
            "pub_date_raw": items.get("PubDate"),
            "volume": items.get("Volume"),
            "doi": items.get("DOI"),
        }
    return None


def main() -> None:
    report: dict = {
        "generated_by": "docs/survey/_check_42800986_two_channels.py",
        "pmid": PMID,
        "why": (
            "42800986 was the only cited PMID supported by a single channel (efetch); "
            "this check adds an independent esummary pull and compares field by field."
        ),
    }
    try:
        summary_root = ET.parse(ESUMMARY).getroot()
        report["esummary_docsum_count"] = len(summary_root.findall("DocSum"))
    except Exception as exc:  # noqa: BLE001 - recorded, not raised
        report["esummary_docsum_count"] = 0
        report["parse_error"] = f"{type(exc).__name__}: {exc}"

    a = efetch_record(PMID)
    b = esummary_record(PMID)
    report["efetch"] = a
    report["esummary"] = b

    if a is None or b is None:
        report["verdict"] = "fetch_failed"
        report["missing"] = [k for k, v in (("efetch", a), ("esummary", b)) if v is None]
    else:
        matches = {}
        for field in ("title", "year", "doi", "volume"):
            va = (a.get(field) or "").strip().lower()
            vb = (b.get(field) or "").strip().lower()
            matches[field] = va == vb
            matches[f"{field}_efetch"] = a.get(field)
            matches[f"{field}_esummary"] = b.get(field)
        # Journal names differ by design: efetch gives "Journal/Title" verbatim,
        # esummary gives FullJournalName; compare only after normalising whitespace/case
        # and accept a match on either full name or ISO abbreviation.
        names = {
            (a.get("journal") or "").strip().lower(),
            (b.get("journal") or "").strip().lower(),
            (a.get("journal_iso") or "").strip().lower(),
            (b.get("journal_iso") or "").strip().lower(),
        } - {""}
        matches["journal_compatible"] = len(names) <= 2 and bool(
            names & {(a.get("journal") or "").strip().lower(), (b.get("journal") or "").strip().lower()}
        )
        matches["journal_efetch"] = a.get("journal")
        matches["journal_esummary"] = b.get("journal")
        report["field_comparison"] = matches
        core = [matches["title"], matches["year"], matches["doi"]]
        report["verdict"] = "consistent" if all(core) else "inconsistent"
        report["core_fields"] = ["title", "year", "doi"]

    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(
        {
            "verdict": report.get("verdict"),
            "esummary_docsum_count": report.get("esummary_docsum_count"),
            "efetch": {k: a[k] for k in ("title", "journal", "year", "doi")} if a else None,
            "esummary": {k: b[k] for k in ("title", "journal", "year", "doi")} if b else None,
        },
        ensure_ascii=False,
        indent=2,
    ))
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
