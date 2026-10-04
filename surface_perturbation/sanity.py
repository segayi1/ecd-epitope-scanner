"""Post-design structural contact checks; never candidate or threshold inputs."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from ..util import sha256_file


def contact_sanity(label_path,expected_sha256,reference,hard,candidate_ids):
    path=Path(label_path)
    if sha256_file(path)!=expected_sha256:
        raise ValueError('frozen contact label fingerprint changed')
    labels=json.loads(path.read_text(encoding='utf-8'))
    if labels.get('numbering')!='UniProt precursor numbering':
        raise ValueError('contact labels require explicit UniProt precursor mapping')
    hard=np.asarray(hard,bool);area=np.asarray(reference['area']);owner=np.asarray(reference['owner'])
    if hard.shape!=(len(candidate_ids),len(area)):
        raise ValueError('contact sanity matrix/reference mismatch')
    report={'label_sha256':expected_sha256,'accession':labels['accession'],
            'source_structure_sha256':labels.get('source_sha256'),
            'excluded_declared_sequence_conflicts':labels.get('excluded_declared_sequence_conflicts',[]),
            'candidate_ids':list(candidate_ids),'contacts':{},
            'interpretation':'post-design structural center-overlap sanity check; not functional detection'}
    for name,item in labels['contacts'].items():
        positions=set(map(int,item['uniprot_positions']))
        contact=np.isin(owner,list(positions)) & (area>0)
        denominator=float(area[contact].sum());covered=np.zeros(len(area),bool);curve=[]
        for k in range(1,13):
            if k<=len(hard): covered|=hard[k-1]
            curve.append({'K':k,'actual_constructs':min(k,len(hard)),
                'contact_center_overlap':float(area[contact&covered].sum()/denominator) if denominator else None})
        report['contacts'][name]={'source':item.get('source',item.get('source_file')),
            'kind':item.get('kind','unspecified_geometric_contact'),
            'heavy_atom_cutoff_a':item.get('heavy_atom_cutoff_a'),
            'target_contact_area_a2':denominator,
            'positions_outside_reliable_target':sorted(positions-set(owner[area>0].tolist())),
            'curve':curve}
    return report
