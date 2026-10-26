import re, json, xml.etree.ElementTree as ET, pathlib

survey = pathlib.Path("docs/survey")
doc = pathlib.Path("docs/extension_directions.md").read_text(encoding="utf-8")

# 1) which PMIDs are cited anywhere in the document?
cited = set(re.findall(r"\b(\d{7,8})\b", doc))
# restrict to ones that look like real PMIDs (present in any snapshot)
paths = {
 "A": survey/"lit_resnapshot_20261016.xml",
 "B": survey/"lit_resnapshot_gap_20261016.xml",
 "C": survey/"_tmp_retry_esummary2.xml",
 "D": survey/"_tmp_retry_efetch.xml",
}
pmids_in = {}
for k, p in paths.items():
    t = ET.parse(p)
    if k == "D":
        ids = [x.findtext(".//PMID") for x in t.getroot()]
    else:
        ids = [i.text if hasattr(i, "text") else i.get("Name")
               for i in t.getroot() for i in [next((c for c in i if c.get("Name") == "PMID"), None)] if i is not None]
        ids = [i for i in ids if i]
    pmids_in[k] = ids

# 2) which PMIDs appear as rows in the section 5 table?
rows = re.findall(r"^\| (\d{8}) \|", doc, re.M)
print("table rows:", len(rows), "unique:", len(set(rows)))

# 3) cited-but-not-in-table, among PMIDs that appear in any snapshot
allids = {i for v in pmids_in.values() for i in v}
missing = sorted(i for i in allids if i in cited and i not in set(rows))
print("cited in text but NOT a table row:", missing)
for m in missing:
    where = [k for k, v in pmids_in.items() if m in v]
    print(" ", m, "in snapshots:", where)

out = {"snapshot_counts": {k: len(v) for k, v in pmids_in.items()},
       "snapshot_ids": pmids_in,
       "table_row_count": len(rows),
       "table_unique": len(set(rows)),
       "cited_but_not_in_table": missing}
(survey/"_tmp_cite_audit.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
print("wrote docs/survey/_tmp_cite_audit.json")
