"""Parse publication-date fields straight out of the saved PubMed snapshots.

No manual transcription: every value below is read from an on-disk file under
docs/survey/ so that the audit can be re-run with:

    python docs/survey/_audit_dates.py
    python docs/survey/_audit_dates.py --out docs/survey/date_field_audit_<tag>.json

**Overwrite warning**: without --out, a re-run writes to
docs/survey/date_field_audit.json and REPLACES whatever is already there. Keep a
copy (or pass --out) if you need the previous version for comparison.

The point of this script is to settle one question mechanically: for each cited
PMID, which year does each snapshot say, and does the chosen publication year
come from PUBDATE rather than DATECOMPLETED.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime
from pathlib import Path

SURVEY = Path(__file__).resolve().parent

FILES = {
    "esummary_A_lit_resnapshot_20261016.xml": "lit_resnapshot_20261016.xml",
    "esummary_B_lit_resnapshot_gap_20261016.xml": "lit_resnapshot_gap_20261016.xml",
    "esummary_C_tmp_retry_esummary2.xml": "_tmp_retry_esummary2.xml",
    "efetch_D_tmp_retry_efetch.xml": "_tmp_retry_efetch.xml",
}


def _flat(fragment: str) -> str:
    return re.sub(r"\s+", " ", fragment).strip()


def parse_efetch(text: str) -> list[dict]:
    out = []
    for block in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", text, re.S):
        pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", block)
        title = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", block, re.S)
        journal = re.search(r"<Title>(.*?)</Title>", block, re.S)
        year = re.search(r"<JournalIssue>.*?<Year>(\d{4})</Year>", block, re.S)
        month = re.search(r"<JournalIssue>.*?<Month>([A-Za-z]+|\d+)</Month>", block, re.S)
        day = re.search(r"<JournalIssue>.*?<Day>(\d+)</Day>", block, re.S)
        doi = re.search(r'<ArticleId IdType="doi">(.*?)</ArticleId>', block)
        pubdate = re.search(r"<PubDate>(.*?)</PubDate>", block, re.S)
        artdate = re.search(r"<ArticleDate[^>]*>(.*?)</ArticleDate>", block, re.S)
        completed = re.search(r"<DateCompleted>(.*?)</DateCompleted>", block, re.S)
        vol = re.search(r"<Volume>(.*?)</Volume>", block)
        issue = re.search(r"<Issue>(.*?)</Issue>", block)
        pages = re.search(r"<MedlinePgn>(.*?)</MedlinePgn>", block)
        out.append(
            {
                "pmid": pmid.group(1) if pmid else None,
                "title": _flat(title.group(1)) if title else None,
                "journal": _flat(journal.group(1)) if journal else None,
                "volume": vol.group(1) if vol else None,
                "issue": issue.group(1) if issue else None,
                "pages": pages.group(1) if pages else None,
                "journal_issue_year": year.group(1) if year else None,
                "journal_issue_month": month.group(1) if month else None,
                "journal_issue_day": day.group(1) if day else None,
                "doi": doi.group(1) if doi else None,
                "pub_date_raw": _flat(pubdate.group(1)) if pubdate else None,
                "article_date_raw": _flat(artdate.group(1)) if artdate else None,
                "date_completed_raw": _flat(completed.group(1)) if completed else None,
            }
        )
    return out


def parse_esummary(text: str) -> list[dict]:
    out = []
    for docsum in re.findall(r"<DocSum>.*?</DocSum>", text, re.S):
        pmid = re.search(r'<Id>(\d+)</Id>', docsum)
        items = dict(re.findall(r'<Item Name="([^"]+)"[^>]*>(.*?)</Item>', docsum, re.S))
        out.append(
            {
                "pmid": pmid.group(1) if pmid else None,
                "title": _flat(items.get("Title", "")) or None,
                "journal": items.get("FullJournalName") or items.get("Source"),
                "pubdate_esummary": _flat(items.get("PubDate", "")) or None,
                "epubdate_esummary": _flat(items.get("EPubDate", "")) or None,
                "doi": None,
            }
        )
    return out


def year_from_pubdate(raw: str | None) -> str | None:
    """The publication year is the year inside PubDate -- never DateCompleted."""
    if not raw:
        return None
    m = re.search(r"(\d{4})", raw)
    return m.group(1) if m else None


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Re-derive publication years from the on-disk PubMed snapshots. "
            "Without --out the default report at docs/survey/date_field_audit.json "
            "is REPLACED."
        )
    )
    parser.add_argument(
        "--out",
        default=None,
        help=(
            "write the report here instead of replacing "
            "docs/survey/date_field_audit.json; an existing file at --out is "
            "refused unless --force is also given"
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="allow --out to replace an existing file",
    )
    args = parser.parse_args(argv)

    report: dict = {"date_policy": {}, "per_file_counts": {}, "by_pmid": {}}

    report["date_policy"] = {
        "publication_year_source": "PubDate (fall back to JournalIssue Year, then ArticleDate)",
        "never_use": "DateCompleted -- this is the MEDLINE processing date, not the publication date",
        "reason": (
            "PMID 27401684 was cited as 2017 because DateCompleted 2017-10-30 was read as the "
            "publication year; its PubDate is 2016 Jul 11."
        ),
    }

    efetch_records: list[dict] = []
    for label, filename in FILES.items():
        path = SURVEY / filename
        text = path.read_text(encoding="utf-8")
        if label.startswith("efetch"):
            records = parse_efetch(text)
            kind = "PubmedArticle"
        else:
            records = parse_esummary(text)
            kind = "DocSum"
        report["per_file_counts"][label] = {
            "file": filename,
            "record_kind": kind,
            "record_count": len(records),
            "pmids": [r["pmid"] for r in records],
        }
        if label.startswith("efetch"):
            efetch_records = records
        for rec in records:
            entry = report["by_pmid"].setdefault(
                rec["pmid"],
                {"pmid": rec["pmid"], "seen_in": [], "esummary": [], "efetch": None},
            )
            entry["seen_in"].append(label)
            if label.startswith("efetch"):
                entry["efetch"] = rec
            else:
                entry["esummary"].append(
                    {
                        "source": label,
                        "title": rec["title"],
                        "journal": rec["journal"],
                        "pubdate_esummary": rec["pubdate_esummary"],
                        "epubdate_esummary": rec["epubdate_esummary"],
                    }
                )

    for entry in report["by_pmid"].values():
        ef = entry["efetch"]
        if ef:
            entry["chosen_publication_year"] = year_from_pubdate(
                ef["pub_date_raw"]
            ) or ef["journal_issue_year"]
            entry["chosen_publication_year_source"] = (
                "efetch PubDate" if ef["pub_date_raw"] else "efetch JournalIssue Year"
            )
        else:
            years = []
            for es in entry["esummary"]:
                y = year_from_pubdate(es["pubdate_esummary"])
                if y:
                    years.append(y)
            entry["chosen_publication_year"] = years[0] if years else None
            entry["chosen_publication_year_source"] = (
                "esummary PubDate" if years else None
            )
            if len(set(years)) > 1:
                entry["pubdate_conflict"] = sorted(set(years))

    report["by_pmid"] = dict(sorted(report["by_pmid"].items()))
    report["summary"] = {
        "distinct_pmids": len(report["by_pmid"]),
        "efetch_only": sorted(
            p for p, e in report["by_pmid"].items() if not e["esummary"]
        ),
        "esummary_only": sorted(
            p for p, e in report["by_pmid"].items() if e["efetch"] is None
        ),
    }

    if args.out:
        out = Path(args.out)
        if not out.is_absolute():
            out = SURVEY / out if (SURVEY / out).parent.exists() else Path(args.out)
        if out.exists() and not args.force:
            raise SystemExit(
                f"refusing to overwrite {out}; pass --force or choose another --out"
            )
    else:
        out = SURVEY / "date_field_audit.json"
        print("NOTE: writing to the default path REPLACES the existing report.")

    report["generated_at"] = datetime.now().astimezone().isoformat(timespec="seconds")
    report["written_to"] = str(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")

    print("per_file_counts:")
    for label, info in report["per_file_counts"].items():
        print(f"  {label}: {info['record_count']} {info['record_kind']} {info['pmids']}")
    print("distinct PMIDs:", report["summary"]["distinct_pmids"])
    print("efetch-only:", report["summary"]["efetch_only"])
    probe = report["by_pmid"].get("27401684", {})
    print("27401684 efetch:", json.dumps(probe.get("efetch"), ensure_ascii=False))
    print("wrote", out)


if __name__ == "__main__":
    main()
