"""Apply the round-7 factual corrections to docs/extension_directions.md.

Each replacement asserts it hits exactly once, so a silent miss fails the run
instead of leaving a half-corrected document.

Run from the repository root:
    python -X utf8 docs/survey/_fix_round7.py
"""

from __future__ import annotations

import pathlib
import sys

DOC = pathlib.Path("docs/extension_directions.md")

# (old, new, why)
EDITS = [
    (
        "方向可以分成三档:能用本机已有冻结数据立刻做的(方向 1–4)、",
        "方向可以分成三档:能用本机已有冻结数据立刻做的(方向 2–4;原方向 1 已移入附录 A)、",
        "方向 1 已移出清单",
    ),
    (
        "- `sectoring.target_sector_residues`(=12,被 `__post_init__` 限制为 10–15 的整数)、`preferred_minimum/maximum_sector_residues`(10/15)",
        "- `sectoring.target_sector_residues`(=12,被 `__post_init__` 限制为 **7–15** 的整数,见 `pipeline.py:119-124`,越界报错原文 `Display target must be an integer from 7 to 15 residues`)、`preferred_minimum/maximum_sector_residues`(10/15,这是**偏好目标**不是校验区间,12 落在 10–15 内)",
        "校验区间实测为 7–15(pipeline.py:122),原写 10–15 有误",
    ),
    (
        "- 该文第 9.1 节的 PD-L1 数字(215 规范化残基、476 组合、60×5、105 节点 10 区、10 个专属首选、8 个正式面板)与本文方向 1 引用的验收收据一致。",
        "- 该文第 9.1 节的 PD-L1 数字(215 规范化残基、476 组合、60×5、105 节点 10 区、10 个专属首选、8 个正式面板)与本文**附录 A**引用的验收收据一致。",
        "方向 1 已移入附录 A",
    ),
    (
        "**本文新增但该文未含的一项**:方向 1 记录的那个导出契约缺口(`panel_coverage_matrix.csv` 两列为空)。这是本次侦察的实测结果,已作为未解释项登记并定位到 `exporting.py:221-226`;是否修复需另立任务。",
        "**本文新增但该文未含的一项**:附录 A 记录的导出契约缺口(`panel_coverage_matrix.csv` 两列为空)。这是本次侦察的实测结果,已登记为缺陷修复任务并定位到 `exporting.py:332-336`(表头 `:335`)与静默取空的 `exporting.py:200-205`,完整登记见 `docs/survey/export_contract_gap.md`。(此前把根因写成 `exporting.py:221-226` 是错的,该处是 `assay_profile` 的写法。)",
        "根因行号更正 + 方向 1 移到附录 A",
    ),
    (
        "5. 方向 1–4 的判据是**本次新写的**,尚未执行过;执行它们需要另立任务并获得授权。",
        "5. 方向 2–4 的判据是**本次新写的**,尚未执行过;执行它们需要另立任务并获得授权。",
        "方向 1 已移出清单",
    ),
    (
        "要补结合侧,文献给出三条路线,全部需要抗体-抗原复合物结构:复合物下的图签名预测(27216816)、FEP+ 类自由能扰动(38844044,成本高且需离群校正)、共识分类程序(32761082)。",
        "要补结合侧,文献给出三条路线,全部需要抗体-抗原复合物结构:①**工具型**(27216816)—— mCSM-AB 是一个 **web server**,输入抗体-抗原复合物结构、用图签名预测突变引起的**亲和力变化**(快照标题即 `mCSM-AB: a web server for predicting antibody-antigen affinity changes upon mutation with graph-based signatures.`),属「拿来就能算」的实现,不是方法学新路线;②FEP+ 类自由能扰动(38844044,成本高且需离群校正);③共识分类程序(32761082)。三者的定位差异见 §5 表该三行。",
        "27216816 是工具/web server 论文,原写法把它当作方法路线",
    ),
    (
        "| 27216816 | 2016 | Nucleic Acids Res | mCSM-AB:用图签名在**抗体-抗原复合物**上预测突变引起的亲和力变化(工具论文);方向 5「要补结合侧需复合物结构」的三条路线之一 [A](在线先发于 2016-05-23) |",
        "| 27216816 | 2016 | Nucleic Acids Res | mCSM-AB:用图签名在**抗体-抗原复合物**上预测突变引起的亲和力变化;**工具/web server 论文**,不是方法学新路线;方向 7「要补结合侧需复合物结构」三条路线中的**工具型**一条 [A](在线先发于 2016-05-23) |",
        "该文属方向 7,原写「方向 5」是章节归属错误",
    ),
    (
        "| 42620036 | 2026 | bioRxiv(**预印本,未经同行评审**) | AF3 抗体基准:幻觉表位、偏好有序区(title 作 \"Bias for Orderness\")[L]† |",
        "| 42620036 | 2026 | bioRxiv(**预印本,未经同行评审**) | AF3 抗体基准:快照标题原文含 \"Hallucinated Epitopes\" 与 \"a Bias for Orderness\"(本文未通读全文,不据此引任何数值)[L]† |",
        "与 §6 限制 7 的「只按底账转述、未回原文核」保持一致",
    ),
]


def main() -> int:
    text = DOC.read_text(encoding="utf-8")
    failures: list[str] = []
    for old, new, why in EDITS:
        count = text.count(old)
        if count != 1:
            failures.append(f"{count} hits (want 1): {why} :: {old[:70]}")
            continue
        text = text.replace(old, new)

    if failures:
        print("ABORTED, nothing written:")
        for item in failures:
            print(" -", item)
        return 1

    DOC.write_text(text, encoding="utf-8")
    print(f"applied {len(EDITS)} edits")
    print("lines:", text.count("\n") + 1, "bytes:", DOC.stat().st_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
