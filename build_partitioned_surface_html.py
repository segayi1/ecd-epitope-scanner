from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "CTLA4_PDL1_partitioned_surface_models.html"
WEB_OUT = ROOT / "web" / "surface_partitions.html"
VENDOR = ROOT / "web" / "vendor" / "3Dmol-min.js"
WORKBENCH_CSS = ROOT / "web" / "workbench.css"

SPECS = (
    {
        "label": "PD-L1",
        "expected_accession": "Q9NZQ7",
    },
    {
        "label": "CTLA-4",
        "expected_accession": "P16410",
    },
)


def _created_at(value: str) -> datetime:
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return parsed.replace(tzinfo=timezone.utc) if parsed.tzinfo is None else parsed


def latest_run(spec: dict[str, str]) -> tuple[Path, dict, dict]:
    candidates: list[tuple[datetime, Path, dict]] = []
    for base in (ROOT / "workspace_data").iterdir():
        manifest_path = base / "run_manifest.json"
        analysis_path = base / "analysis.json"
        structure_path = base / "structure_for_energy.pdb"
        if not (manifest_path.is_file() and analysis_path.is_file() and structure_path.is_file()):
            continue
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            request = manifest.get("request", {})
            source = manifest.get("source", {})
            accessions = {
                source.get("normalized_accession"),
                request.get("accession"),
                request.get("identifier"),
            }
            if spec["expected_accession"] not in accessions:
                continue
            candidates.append((_created_at(manifest["created_at"]), base, manifest))
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue

    for _, base, manifest in sorted(candidates, reverse=True, key=lambda row: row[0]):
        try:
            analysis = json.loads((base / "analysis.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (
            analysis.get("protein", {}).get("accession") == spec["expected_accession"]
            and analysis.get("regions")
            and analysis.get("residues")
            and (base / "structure_for_energy.pdb").stat().st_size
        ):
            return base, manifest, analysis
    raise FileNotFoundError(f"No complete run found for {spec['expected_accession']}")


def payload(spec: dict[str, str]) -> dict:
    base, manifest, analysis = latest_run(spec)
    protein = analysis["protein"]
    if protein["accession"] != spec["expected_accession"]:
        raise ValueError(f"{base.name}: expected {spec['expected_accession']}, got {protein['accession']}")
    source = manifest.get("source", {})
    structure_id = str(source.get("identifier") or source.get("normalized_accession") or "structure")
    structure = f"{source.get('provider', 'Source')} {structure_id.upper()}"
    input_label = f"{manifest.get('request', {}).get('input_type', 'input')} / {source.get('category', source.get('kind', 'structure'))}"
    residues = [
        {
            "pos": item["uniprot_pos"],
            "chain": item["chain_id"],
            "resi": item["auth_seq_id"],
            "region": item.get("region_id"),
        }
        for item in analysis["residues"]
    ]
    # Do not use ResidueRecord.region_id here: the V4 display projection may
    # retain a different legacy assignment.  RegionRecord.residues is the
    # authoritative membership list for the exported sector.
    residue_by_position = {item["pos"]: item for item in residues}
    regions = [
        {
            "id": item["id"], "name": item["name"], "color": item["color"],
            "n": len(item["residues"]), "area": item["surface_area"],
            "geodesic": item["geodesic_diameter"], "prior": item["epitope_prior"],
            "tier": item["tier"], "auto": item["auto_design_allowed"],
            "flags": item.get("quality_flags", []),
            # Author residue IDs are carried into the HTML so 3Dmol can render
            # an independently colored molecular surface for every sector.
            "atoms": [
                {"chain": residue_by_position[position]["chain"],
                 "resi": residue_by_position[position]["resi"]}
                for position in item["residues"]
                if position in residue_by_position
            ],
        }
        for item in analysis["regions"]
    ]
    return {
        "id": base.name, "label": spec["label"], "accession": protein["accession"],
        "protein_name": protein["protein_name"], "structure": structure,
        "input": input_label, "created": manifest["created_at"],
        "engine": manifest.get("request", {}).get("design_engine", "unknown"),
        "mode": manifest.get("request", {}).get("ml_mode", "unknown"),
        "pdb": (base / "structure_for_energy.pdb").read_text(encoding="utf-8"),
        "residues": residues, "regions": regions,
    }


def main() -> None:
    models = [payload(spec) for spec in SPECS]
    vendor = VENDOR.read_text(encoding="utf-8")
    workbench_css = WORKBENCH_CSS.read_text(encoding="utf-8")
    models_json = json.dumps(models, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")
    template = r'''<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>PD-L1 与 CTLA-4 已分区 3D 表面模型</title><style>__WORKBENCH_CSS__
:root{--paper:var(--bg);--card:var(--surface);--accent:var(--teal-dark)}body{margin:0;color:var(--ink);font:14px/1.5 "Segoe UI","Microsoft YaHei",sans-serif;background:var(--bg)}h1{margin:0;font-size:20px}.site-page-intro{padding:14px 24px 0}.site-page-intro h1{margin:0;font-size:clamp(20px,2vw,27px)}.site-page-intro p{max-width:1100px;margin:5px 0 0;color:var(--muted)}#tabs{display:flex;gap:8px;padding:12px 24px 0;overflow-x:auto}.tab,.btn{border:1px solid var(--line);background:#fff;border-radius:8px;padding:7px 12px;cursor:pointer;color:var(--ink);font:inherit;white-space:nowrap}.tab[aria-selected="true"]{background:var(--teal-dark);color:#fff;border-color:var(--teal-dark)}main{display:grid;grid-template-columns:minmax(0,2.15fr) minmax(315px,.9fr);gap:14px;padding:12px 24px 20px;min-width:0}#stage{position:relative;height:clamp(430px,calc(100vh - 250px),720px);min-width:0;border:1px solid var(--line);border-radius:12px;overflow:hidden;background:#f6f9fa}#viewer{position:absolute;inset:0}.hud{position:absolute;z-index:2;top:12px;left:12px;display:flex;gap:6px;flex-wrap:wrap}.hud .btn.on{color:#fff;background:var(--teal-dark);border-color:var(--teal-dark)}.hint{position:absolute;z-index:2;bottom:10px;left:12px;background:#fffffff0;padding:5px 8px;border-radius:7px;color:var(--muted);font-size:12px}aside{min-width:0;max-height:clamp(430px,calc(100vh - 250px),720px);overflow:auto}.card{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:12px;margin-bottom:12px}.meta{color:var(--muted);font-size:12px}.meta b{color:var(--ink)}h2{font-size:16px;margin:0 0 7px}.region{display:grid;grid-template-columns:15px minmax(0,1fr) auto;gap:8px;align-items:center;padding:7px 5px;border-radius:7px;cursor:pointer}.region:hover,.region.active{background:#eaf4f3}.region.active{outline:1px solid var(--teal)}.swatch{width:14px;height:14px;border-radius:3px}.region small{display:block;color:var(--muted)}.badge{font-size:11px;padding:1px 5px;border-radius:5px;background:#edf0f2;color:#52616a}.warn{background:#fae5ce;color:#8b4d00}.legend{display:flex;gap:10px;flex-wrap:wrap;font-size:12px}.legend i{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px}.note{font-size:12px;color:var(--muted)}details summary{cursor:pointer;font-weight:600}@media(max-width:860px){.site-page-intro{padding:12px 12px 0}#tabs{padding-left:12px;padding-right:12px}main{grid-template-columns:1fr;padding:10px 12px 16px}#stage{height:58vh;min-height:380px}aside{max-height:none}}
</style><style>.site-page-intro{font-size:14px;line-height:1.5;font-weight:400;letter-spacing:normal}@media(min-width:861px){#stage,aside{max-height:clamp(430px,calc(100vh - 280px),720px)}#stage{height:clamp(430px,calc(100vh - 280px),720px)}}</style></head><body>
<header class="site-header"><a class="site-brand" href="/" aria-label="ECD Epitope Workbench 首页"><span class="site-mark" aria-hidden="true">E</span><span>ECD Epitope Workbench</span></a><nav class="site-links" aria-label="主导航"><a href="/">结构设计</a><a href="/experiments">实验项目</a><details class="site-nav-menu"><summary>结构图册</summary><div><a href="/surface_partitions">分区表面</a><a href="/ctla4_six">CTLA-4 结构</a><a href="/pd1-antibodies">PD-1 抗体</a></div></details></nav></header>
<section class="site-page-intro site-heading"><p class="site-kicker">结构图册 · 分区表面</p><h1>PD-L1 与 CTLA-4：已分区的交互式 3D 表面模型</h1><p>单文件离线报告：已内嵌 3Dmol.js、PDB 坐标与本次 Joint V4-ML 分区结果。整条目标链作为一张连续表面绘制；分区交界处由网格颜色插值柔化，视觉过渡不会改变 sector 残基归属。点击右侧 sector 聚焦视角。模型颜色表示显示 sector，不等同于实验确认的抗体表位。</p></section>
<div id="tabs"></div><main><section id="stage"><div id="viewer"></div><div class="hud"><button class="btn on" id="surface">分区表面</button><button class="btn" id="cartoon">卡通 + 分区</button><button class="btn" id="reset">复位视角</button></div><div class="hint">左键旋转 · 滚轮缩放 · 右键平移 · 点击右侧分区聚焦</div></section><aside><div class="card" id="summary"></div><div class="card"><h2>显示分区（sector）</h2><div id="regions"></div></div><div class="card"><details><summary>阅读与解释边界</summary><p class="note">sector 是显式溶剂表面图上的互斥显示归属，用于粗定位探针覆盖组织；V4 的候选 footprint hypothesis 可以彼此重叠。本文件显示输入 WT 结构坐标，不包含突变后的原子构象预测。SASA、先验分数和颜色均不能单独证明有限尺寸 Fab 可达性、结合变化或真实表位。</p></details></div></aside></main>
<script>__VENDOR__</script><script>
const MODELS=__MODELS__;let current=0,viewer=null,viewerVersion=0,representation='surface',selected=null;const $=id=>document.getElementById(id);
function fmt(n,d=3){return Number(n).toFixed(d).replace(/\.0+$/,'').replace(/(\.\d*?)0+$/,'$1')}function model(){return MODELS[current]}function selector(region){const byChain={};region.atoms.forEach(a=>(byChain[a.chain]??=[]).push(a.resi));return Object.entries(byChain).map(([chain,resi])=>({chain,resi}))}
async function renderViewer(){const draw=++viewerVersion,m=model();if(viewer)viewer.clear();$('viewer').innerHTML='';const activeViewer=viewer=$3Dmol.createViewer($('viewer'),{backgroundColor:'#f4f0e7',antialias:true});activeViewer.addModel(m.pdb,'pdb');if(representation==='surface'){const colors=['#c8cdd1',...m.regions.map(r=>r.color)],scheme={gradient:'linear_'+colors.join('_'),min:0,max:colors.length},tasks=[],chains=[...new Set(m.residues.map(r=>r.chain))];chains.forEach(chain=>activeViewer.mapAtomProperties(atom=>{atom.properties=atom.properties||{};atom.properties.surface_region_color=0},{chain}));m.regions.forEach((r,i)=>selector(r).forEach(sel=>activeViewer.mapAtomProperties(atom=>{atom.properties=atom.properties||{};atom.properties.surface_region_color=i+1},sel)));chains.forEach(chain=>tasks.push(activeViewer.addSurface($3Dmol.SurfaceType.MS,{map:{prop:'surface_region_color',scheme},opacity:1,smoothness:3},{chain},{chain})));await Promise.all(tasks);if(draw!==viewerVersion)return}else{activeViewer.setStyle({},{cartoon:{color:'#c8cdd1'}});m.regions.forEach(r=>selector(r).forEach(sel=>activeViewer.setStyle(sel,{cartoon:{color:r.color}})))}const focused=selected===null?null:m.regions.find(r=>r.id===selected);if(focused){const sel=selector(focused)[0];if(sel)activeViewer.zoomTo(sel);activeViewer.zoom(1.45)}else{activeViewer.zoomTo();activeViewer.zoom(1.28)}activeViewer.render()}function focusRegion(id){selected=selected===id?null:id;renderRegions();if(selected!==null){const r=model().regions.find(x=>x.id===selected),s=selector(r);if(s.length){viewer.zoomTo(s[0]);viewer.zoom(1.45);viewer.render()}}}function renderTabs(){$('tabs').innerHTML=MODELS.map((m,i)=>`<button class="tab" aria-selected="${i===current}" onclick="switchModel(${i})">${m.label} · ${m.accession}</button>`).join('')}function renderSummary(){const m=model();$('summary').innerHTML=`<h2>${m.label} <span class="meta">${m.accession}</span></h2><div class="meta"><b>${m.protein_name}</b><br>结构：${m.structure}<br>输入：${m.input}<br>运行：${m.engine} / ${m.mode}<br>运行 ID：${m.id}<br>生成：${new Date(m.created).toLocaleString('zh-CN',{timeZone:'Asia/Taipei'})} 台北时间<br>显示 sector：${m.regions.length}</div><div class="legend" style="margin-top:9px"><span><i style="background:#c8cdd1"></i>未归属/非选 sector</span><span><i style="background:#1565c0"></i>sector 色见下表</span></div>`}function renderRegions(){const m=model();$('regions').innerHTML=m.regions.map(r=>{const flags=r.flags||[],q=flags.length?`<span class="badge warn" title="${flags.join(', ')}">QC 标记</span>`:`<span class="badge">${r.tier}</span>`;return `<div class="region ${selected===r.id?'active':''}" onclick="focusRegion(${r.id})"><span class="swatch" style="background:${r.color}"></span><span><b>R${r.id} · ${r.name}</b><small>${r.n} 残基 · ${fmt(r.area)} Å² · 测地直径 ${fmt(r.geodesic)} Å · 先验 ${fmt(r.prior,4)}${r.auto?'':' · 不自动设计'}${flags.length?' · '+flags.join(', '):''}</small></span>${q}</div>`}).join('')}function renderAll(){renderTabs();renderSummary();renderRegions();renderViewer();$('surface').classList.toggle('on',representation==='surface');$('cartoon').classList.toggle('on',representation==='cartoon')}function switchModel(i){current=i;selected=null;renderAll()}$('surface').onclick=()=>{representation='surface';renderAll()};$('cartoon').onclick=()=>{representation='cartoon';renderAll()};$('reset').onclick=()=>{viewer.zoomTo();viewer.zoom(1.28);viewer.render()};renderAll();
</script></body></html>'''
    html = template.replace("__WORKBENCH_CSS__", workbench_css).replace("__VENDOR__", vendor).replace("__MODELS__", models_json)
    html = html.replace('href="/surface_partitions"', 'href="/assets/surface_partitions.html"')
    html = html.replace('href="/ctla4_six"', 'href="/assets/ctla4_six_mutants.html"')
    for output in (OUT, WEB_OUT):
        output.write_text(html, encoding="utf-8")
    print(f"Created {OUT.name} and {WEB_OUT.relative_to(ROOT)}: {OUT.stat().st_size} bytes each; models={len(models)}")


if __name__ == "__main__":
    main()
