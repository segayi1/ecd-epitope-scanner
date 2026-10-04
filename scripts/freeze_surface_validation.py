"""Overlay an explicit source snapshot onto an existing detached checkout.

Never copies runtime/history/model folders. The manifest makes the overlay
(including required pre-existing user dependencies) distinct from Git HEAD.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import shutil
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('destination',type=Path)
    args=parser.parse_args();source=Path(__file__).resolve().parents[1];target=args.destination.resolve()
    if target==source or source in target.parents or not (target/'.git').is_file():
        raise ValueError('requires a separate, existing detached Git worktree')
    paths=list(source.glob('*.py'))
    for directory in ('surface_perturbation','joint_design_v4','joint_design','integrations','surface_patch','scripts','schemas','web','tests'):
        for path in (source/directory).rglob('*'):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix in {'.py','.js','.mjs','.html','.css','.json','.csv','.yaml','.toml','.md','.svg','.png','.pdb','.cif'}:
                paths.append(path)
    paths.extend(source/name for name in ('AGENTS.md','pyproject.toml','requirements.txt',
        'docs/project_memory.md','docs/ctla4_coverage_v3_implementation_audit_20261004.md',
        'docs/ctla4_stability_and_physical_evidence_research_20261004.md'))
    manifest={}
    for path in sorted(set(paths)):
        if not path.is_file(): continue
        relative=path.relative_to(source);destination=target/relative
        destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,destination)
        manifest[relative.as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    (target/'source_overlay_manifest.json').write_text(json.dumps({'source':str(source),'files':manifest,
        'interpretation':'detached HEAD plus explicit current source overlay; no clean-HEAD claim'},indent=2),encoding='utf-8')
    print(json.dumps({'destination':str(target),'files':len(manifest)}))
if __name__=='__main__':main()
