import json
import pathlib
import sys

NEW_BLOCK = {
    "snapshots": [
        "docs/survey/lit_resnapshot_20261016.xml(19 条 DocSum,原始返回)",
        "docs/survey/lit_resnapshot_gap_20261016.xml(4 条 DocSum)",
        "docs/survey/_tmp_retry_esummary2.xml(10 条 DocSum)",
        "docs/survey/_tmp_retry_efetch.xml(5 条 PubmedArticle)",
        "docs/survey/_fetch_esummary_42800986.xml(1 条 DocSum)",
    ],
    "distinct_pmids_in_snapshot_union": 34,
    "distinct_pmids_cited_in_artifact": 24,
    "cited_pmids_without_snapshot": [],
    "note": "rev3 的 23/21 是当时的量(只算了两份 esummary 快照),与 docs/survey/snapshot_article_count.json 的 34/24 不是同一个口径;2026-10-17 起以脚本重算的 34(distinct)/24(被引)为准,差值 10 条是快照命中了但正文未引用(清单见该 json 的 in_snapshot_but_not_cited)。",
}


def main() -> int:
    target = pathlib.Path("docs/survey/self_check.json")
    if not target.exists():
        print(f"missing {target}")
        return 1
    doc = json.loads(target.read_text(encoding="utf-8"))

    hit = None
    for item in doc.get("criteria_check", []):
        if item.get("criterion", "").startswith("(4)"):
            hit = item
            break
    if hit is None:
        print("criterion (4) not found")
        return 1

    before = hit.get("evidence_added_this_revision", {})
    before_union = before.get("distinct_pmids_in_snapshot_union")
    if before_union == NEW_BLOCK["distinct_pmids_in_snapshot_union"]:
        print("already updated; nothing to do")
        return 0

    hit["evidence_added_this_revision"] = NEW_BLOCK
    hit["note"] = (
        "rev2 把本条记为 pass 是**按声称而非按落盘凭据**判的,属自查缺陷,已在 revision 3 更正。"
        "2026-10-17 再更正一次:本条 rev3 记的 23(union)/21(被引)只算了两份 esummary 快照,"
        "未含 C/D/E 三份补拉快照;按 docs/survey/snapshot_article_count.json 脚本重算应为 34/24,"
        "差集仍为空(24 条被引**每一条**都有落盘快照)。verdict 与 scope_of_this_pass 未改动。"
    )
    hit["revision_of_this_update"] = 4
    doc["revision"] = max(int(doc.get("revision", 0)), 3) + 1

    target.write_text(
        json.dumps(doc, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"updated: revision={doc['revision']}")
    print(
        "evidence_now:",
        json.dumps(hit["evidence_added_this_revision"], ensure_ascii=False),
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
