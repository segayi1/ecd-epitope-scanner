"""Round-7 part 2: align the 42620036 wording and add per-direction verdict markers.

Two gaps from the sixth evaluation round:
  * u3 - direction 6 leaned on 42620036 for "hallucinated epitopes / Bias for Orderness"
    while §5 and §6 list that item as "only transcribed from the ledger, never
    re-checked against the source". Both statements must be the same statement.
  * criterion (2) - directions 8 and 9 never carry the explicit
    "可 / 需外部数据 / 需湿实验" marker that direction 7 has.

Every replacement asserts exactly one hit.

Run from the repository root:
    python -X utf8 docs/survey/_fix_round7b.py
"""

from __future__ import annotations

import pathlib
import sys

DOC = pathlib.Path("docs/extension_directions.md")

EDITS = [
    (
        "另有 AF3 在抗体复合物上\"幻觉表位、偏好有序区\"(原文作 \"a Bias for Orderness\")的预印本结论 [PMID 42620036, `method_landscape.json`] —— 这与本仓库\"AF3 必须带置信度文件、不得当实验结构\"的要求方向一致。",
        "另有 AF3 在抗体复合物上的预印本结论,其**快照标题原文**为 \"A Structural Antibody Benchmark of AlphaFold3 reveals **Hallucinated Epitopes** and **a Bias for Orderness**\" [PMID 42620036]。**说明**:本文只读到该条的 esummary 标题层,未取摘要或全文,因此这里引用的是**标题自身的措辞**,不引任何数值,也不把它当作已验证的实测结论(与 §5 表该行、§6 限制 7 的记载一致)。它与本仓库\"AF3 必须带置信度文件、不得当实验结构\"的要求方向一致。",
        "u3:方向 6 与 §5/§6 对 42620036 的说法必须一致",
    ),
    (
        "**为什么慢**:需要冻结输入、构建预算、模型资源与判定合同,并且必须区分\"有结构/候选准备/真实评分/物理解析/实验确认\"五个层次,避免把准备完成写成验证完成。",
        "**为什么慢**:需要冻结输入、构建预算、模型资源与判定合同,并且必须区分\"有结构/候选准备/真实评分/物理解析/实验确认\"五个层次,避免把准备完成写成验证完成。\n\n**本机可验证性标记:需外部数据 + 需算力(当前不可做)**。本机虽有 PD-L1 冻结运行,但只有一个蛋白、一套引擎参数,不构成\"跨蛋白、跨引擎\"比较所需的输入集;冻结新的蛋白输入与重跑评分都需要模型资源。**不得**把它写成\"现在可立即执行\"。",
        "判据(2):方向 8 缺可验证性标记",
    ),
    (
        "**为什么慢**:这是唯一能把\"几何代理\"变成\"真实结合后果\"的路径,但它需要湿实验,且登记状态 COMPLETE 只说明登记事务状态,不能称实验分析完成。",
        "**为什么慢**:这是唯一能把\"几何代理\"变成\"真实结合后果\"的路径,但它需要湿实验,且登记状态 COMPLETE 只说明登记事务状态,不能称实验分析完成。\n\n**本机可验证性标记:需湿实验(当前不可做)**。本机侧只做到\"真实结果表已登记\"(162 个样本、3 份 FlowJo 百分比导入),人工 QC 为 0、效应分析为 0;把几何代理校准成真实结合后果必须新做湿实验,本机任何脚本改动都不构成这条方向的验证。",
        "判据(2):方向 9 缺可验证性标记",
    ),
]


def main() -> int:
    text = DOC.read_text(encoding="utf-8")
    failures: list[str] = []
    for old, new, why in EDITS:
        count = text.count(old)
        if count != 1:
            failures.append(f"{count} hits (want 1): {why}")
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
    for marker in ("本机可验证性标记:需外部数据 + 需算力", "本机可验证性标记:需湿实验"):
        print(f"{marker} -> {text.count(marker)} hit(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
