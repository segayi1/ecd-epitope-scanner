"""Real-weight regression, run with the installed WSL SPURS environment."""
import argparse
import copy
import json
from pathlib import Path
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=Path)
    parser.add_argument('source_job')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    sys.path[:0] = [str(root/'vendor/spurs'), str(root/'integrations/workers')]
    import torch
    from spurs_gpu_worker import _load_multi
    from spurs.inference import parse_pdb, parse_pdb_for_mutation, _seed_everything
    from gpu_runtime import serialized_gpu

    source = root/'workspace_data'/args.source_job
    analysis = json.loads((source/'analysis.json').read_text())
    groups = {}
    for order in (2, 3):
        candidates = [c for c in analysis['joint_v4_ml']['candidates'] if c['order']==order][:3]
        assert len(candidates)==3
        groups[order] = [[f"{m['residue']['one_letter']}{m['residue']['uniprot_pos']}{m['mutant_aa']}" for m in c['mutations']] for c in candidates]
    with serialized_gpu(), torch.inference_mode():
        model, cfg = _load_multi(root/'models/spurs/spurs_multi', 'cuda')
        model.eval()
        seed = int(cfg['train']['seed'])
        pdb = parse_pdb(str(source/'structure_for_energy.pdb'), 'encoding_regression', 'A', cfg, device='cuda')
        scores={}; timings={}
        for mode in ('legacy_per_order', 'shared_wt'):
            _seed_everything(seed)
            started=time.perf_counter()
            work=copy.deepcopy(pdb)
            representation = model.encode_wt(work) if mode=='shared_wt' else None
            values=[]
            for order, mutations in groups.items():
                ids, appends=parse_pdb_for_mutation(mutations)
                work['mut_ids']=ids.to('cuda');work['append_tensors']=appends.to('cuda')
                pred=model(work) if mode=='legacy_per_order' else model.decode_mutations(work,representation)
                values.extend(pred.detach().float().cpu().reshape(-1).tolist())
            torch.cuda.synchronize()
            scores[mode]=values;timings[mode]=time.perf_counter()-started
        delta=max(abs(a-b) for a,b in zip(scores['legacy_per_order'],scores['shared_wt']))
        receipt={'status':'PASS' if delta<=1e-6 else 'FAIL','source_job':args.source_job,
                 'candidate_mutations':groups,'scores':scores,'maximum_absolute_difference':delta,
                 'tolerance':1e-6,'inference_seconds':timings,'device':torch.cuda.get_device_name(0),
                 'peak_reserved_bytes':torch.cuda.max_memory_reserved(),
                 'scope':'Same real weights/input, three doubles and three triples; no FoldX calls.'}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt))
    return 0 if receipt['status']=='PASS' else 1


if __name__=='__main__':
    raise SystemExit(main())
