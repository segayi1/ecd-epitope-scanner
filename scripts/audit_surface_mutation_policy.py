"""Read-only requalification of saved panel identities against current policy."""
import argparse
import json
from pathlib import Path
from ..stability import AUTO_EXCLUDED_MUTANTS
from ..util import sha256_file, write_json


def audit(folder):
    folder = Path(folder)
    summary_path = folder / 'surface_perturbation_preview.json'
    summary = json.loads(summary_path.read_text(encoding='utf-8')) if summary_path.exists() else {}
    ids = summary.get('nested_verified_candidate_ids', summary.get('selected_candidate_ids', []))
    if not ids:
        return {'status': 'INCOMPLETE_POLICY_EVIDENCE', 'panel_ids': [],
                'reason': 'No saved final panel identities; do not infer acceptance from an empty set.'}
    path = folder / 'surface_perturbation' / 'paired_manifest.json'
    if not path.exists():
        return {'status': 'INCOMPLETE_POLICY_EVIDENCE', 'panel_ids': ids}
    data = json.loads(path.read_text(encoding='utf-8'))
    violations = []
    unresolved = []
    for candidate_id in ids:
        pairs = data.get(candidate_id, [])
        if not pairs:
            unresolved.append(candidate_id)
            continue
        for mutation in pairs[0].get('mutations', []):
            if mutation[-1] in AUTO_EXCLUDED_MUTANTS or mutation[0] in 'CGP':
                violations.append({'candidate_id': candidate_id, 'mutation': mutation})
    return {'status': 'INCOMPLETE_MUTATION_POLICY' if violations else
            'INCOMPLETE_POLICY_EVIDENCE' if unresolved else 'CURRENT_POLICY_PASS',
            'excluded_mutants': sorted(AUTO_EXCLUDED_MUTANTS),
            'panel_ids': ids, 'violations': violations, 'unresolved': unresolved,
            'summary_sha256': sha256_file(summary_path),
            'paired_manifest_sha256': sha256_file(path),
            'interpretation': 'Policy pass alone does not establish physical or budget acceptance.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('folder', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.folder)
    write_json(args.output, result)
    print(json.dumps(result))


if __name__ == '__main__': main()
