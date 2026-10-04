"""Rebuild the standalone Joint V4-ML method schematic (PNG and SVG)."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

font_path = Path(r"C:\Windows\Fonts\msyh.ttc")
font_family = font_manager.FontProperties(fname=str(font_path)).get_name() if font_path.is_file() else "DejaVu Sans"
plt.rcParams.update({"font.family": font_family, "svg.fonttype": "none"})

NAVY = "#17324D"
INK = "#203247"
MUTED = "#5E7185"
TEAL = "#168C88"
BLUE = "#3978B5"
GOLD = "#D69A2D"
CORAL = "#C86A5B"
PALE = "#F4F7FA"
LINE = "#CFD9E2"
WHITE = "#FFFFFF"

fig, ax = plt.subplots(figsize=(19.2, 10.8), facecolor=WHITE)
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

ax.text(0.035, 0.955, "Joint V4-ML | ECD 突变探针面板设计", fontsize=23,
        fontweight="bold", color=NAVY, va="top")
ax.text(0.036, 0.915,
        "结构映射  →  表面/表位/稳定性证据  →  连通区域分区  →  组合复核  →  strict 门控",
        fontsize=12, color=MUTED, va="top")


def card(x, y, w, h, number, title, accent=NAVY):
    patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.014",
                           facecolor=WHITE, edgecolor=LINE, linewidth=1.25)
    ax.add_patch(patch)
    ax.add_patch(FancyBboxPatch((x, y + h - 0.075), w, 0.075,
                                boxstyle="round,pad=0.008,rounding_size=0.014",
                                facecolor=accent, edgecolor=accent, linewidth=0.5))
    ax.add_patch(plt.Circle((x + 0.027, y + h - 0.037), 0.015, color=WHITE, zorder=4))
    ax.text(x + 0.027, y + h - 0.037, str(number), fontsize=9, color=accent,
            fontweight="bold", ha="center", va="center", zorder=5)
    ax.text(x + 0.052, y + h - 0.037, title, fontsize=12.2, color=WHITE,
            fontweight="bold", va="center")


def block(x, y, w, h, title, body, color, body_size=9.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.007,rounding_size=0.010",
                                facecolor=PALE, edgecolor=LINE, linewidth=0.75))
    ax.add_patch(FancyBboxPatch((x, y), 0.008, h,
                                boxstyle="round,pad=0.002,rounding_size=0.006",
                                facecolor=color, edgecolor=color, linewidth=0))
    ax.text(x + 0.018, y + h - 0.022, title, fontsize=10.1, fontweight="bold",
            color=color, va="top")
    ax.text(x + 0.018, y + h - 0.055, body, fontsize=body_size, color=INK,
            va="top", linespacing=1.35)


cards = [
    (0.025, 0.245, 0.170, 0.585),
    (0.210, 0.245, 0.275, 0.585),
    (0.500, 0.245, 0.175, 0.585),
    (0.690, 0.245, 0.170, 0.585),
    (0.875, 0.245, 0.100, 0.585),
]
titles = ["输入与映射", "并行证据", "表面分区", "候选与复核", "正式筛选"]
accents = [NAVY, BLUE, TEAL, GOLD, CORAL]
for index, ((x, y, w, h), title, accent) in enumerate(zip(cards, titles, accents, strict=True), 1):
    card(x, y, w, h, index, title, accent)

# Step 1
x, y, w, h = cards[0]
block(x + 0.012, 0.650, w - 0.024, 0.145, "目标输入",
      "ECD 序列与坐标\nPDB/mmCIF + 目标链\n测定场景：cell_surface / purified", BLUE, 8.8)
block(x + 0.012, 0.465, w - 0.024, 0.155, "映射 QC",
      "校验链、WT 与残基数\n处理 author 编号/插入码\n拒绝歧义或 WT 不一致", TEAL, 8.7)
block(x + 0.012, 0.280, w - 0.024, 0.155, "规范化上下文",
      "单链连续编号 1…L\n生成 residue_mapping.csv\n原始结构保留为生物学上下文", NAVY, 8.4)

# Step 2: three parallel evidence channels
x, y, w, h = cards[1]
block(x + 0.012, 0.645, w - 0.024, 0.150, "表面可及性｜SASA / rSASA",
      "Native + ProCaliper（Shrake–Rupley）\n探针 1.4 Å；100 点；Tien MaxASA\n默认 rSASA ≥ 0.25 且侧链 SASA ≥ 20 Å²", TEAL, 8.4)
block(x + 0.012, 0.465, w - 0.024, 0.150, "表位倾向｜DiscoTope-3.0",
      "ESM-IF1 结构表示 + 阳性/未标记学习\n保留原始分数、标签与 provenance\n仅作先验排序，不作为硬门槛", BLUE, 8.4)
block(x + 0.012, 0.285, w - 0.024, 0.150, "单点稳定性｜SPURS single",
      "ESM2 + ProteinMPNN 特征重连\n逐残基输出 19 种非 WT 替换 ΔΔG\n默认可进入 strict：≤ 0.5 kcal/mol", GOLD, 8.4)

# Step 3
x, y, w, h = cards[2]
block(x + 0.012, 0.645, w - 0.024, 0.150, "表面残基图",
      "高暴露位点作节点\n表面边界位点可作路径中继\n边权表示表面测地距离", TEAL, 8.7)
block(x + 0.012, 0.465, w - 0.024, 0.150, "先验岛与连通组件",
      "高 DiscoTope 位点形成软先验岛\n保留低先验 COVERAGE_ONLY 区域\n不把岛屿等同真实抗体 footprint", BLUE, 8.3)
block(x + 0.012, 0.285, w - 0.024, 0.150, "确定性 sector",
      "图 k-center / minimax medoid\n测地 Voronoi；目标数默认 7\n不连通组件多时保留各区；直径 QC：35 / 30 Å", NAVY, 8.1)

# Step 4
x, y, w, h = cards[3]
block(x + 0.012, 0.645, w - 0.024, 0.150, "位点保护与菜单",
      "二硫键/Cys、主链风险、糖基化\n及用户标注位点作硬排除\n分层保留表位先验与回退选项", CORAL, 8.2)
block(x + 0.012, 0.465, w - 0.024, 0.150, "组合候选搜索",
      "同一连通 sector 内组合 2/3/4 位点\n最大测地跨度 24 Å\n直接 / 9 / 12 / 15 Å 覆盖情景", GOLD, 8.2)
block(x + 0.012, 0.285, w - 0.024, 0.150, "组合稳定性",
      "SPURS multi 预测组合 ΔΔG / 非加和\nFoldX RepairPDB + BuildModel ×5\n中位数汇总；短名单上限 42", BLUE, 8.0)

# Step 5 compact card
x, y, w, h = cards[4]
block(x + 0.010, 0.645, w - 0.020, 0.150, "完整性",
      "worker 必须\n成功且完整\n部分评分不放行", BLUE, 8.5)
block(x + 0.010, 0.465, w - 0.020, 0.150, "候选门控",
      "单点/组合通过\n无冲突、无 QC 失败\n模型分歧则阻断", CORAL, 8.2)
block(x + 0.010, 0.285, w - 0.020, 0.150, "面板输出",
      "确定性 sector 平衡\n最多 7 个主构建体\nstrict 或 BLOCKED", TEAL, 8.2)

# Between-stage arrows
for i in range(len(cards) - 1):
    x0 = cards[i][0] + cards[i][2] + 0.003
    x1 = cards[i + 1][0] - 0.003
    ax.add_patch(FancyArrowPatch((x0, 0.535), (x1, 0.535), arrowstyle="-|>",
                                 mutation_scale=14, linewidth=1.7, color=MUTED))

# Footer: reproducibility and interpretive boundary
ax.add_patch(FancyBboxPatch((0.025, 0.075), 0.950, 0.115,
                            boxstyle="round,pad=0.008,rounding_size=0.012",
                            facecolor="#EEF3F7", edgecolor=LINE, linewidth=0.9))
ax.text(0.045, 0.157, "可复现审计", fontsize=10.2, color=NAVY, fontweight="bold", va="center")
ax.text(0.155, 0.157,
        "结构/代码 SHA-256 · 参数快照 · worker 命令/状态/错误 · 模型/二进制摘要 · 候选排除原因",
        fontsize=9.2, color=INK, va="center")
ax.text(0.045, 0.112, "解释边界", fontsize=10.2, color=CORAL, fontweight="bold", va="center")
ax.text(0.155, 0.112,
        "表位倾向 ≠ 抗体特异性结合；预测稳定性 ≠ 实测表达/折叠；当前未建模 cell-surface Fab/IgG 膜可达性。",
        fontsize=9.2, color=INK, va="center")
ax.text(0.965, 0.035, "V4 thresholds are engineering heuristics, not calibrated probabilities",
        fontsize=8.0, color=MUTED, ha="right", va="center", style="italic")

fig.savefig(OUT / "joint_v4_method.png", dpi=240, bbox_inches="tight", facecolor=WHITE)
fig.savefig(OUT / "joint_v4_method.svg", bbox_inches="tight", facecolor=WHITE)
fig.savefig(OUT / "joint_v4_method.pdf", bbox_inches="tight", facecolor=WHITE)
plt.close(fig)
