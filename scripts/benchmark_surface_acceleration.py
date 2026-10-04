"""Engineering timing only: fixed PQR/grid, independent fresh serial/parallel caches."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from ..surface_perturbation.tools import PhysicalTools,read_pqr,make_grid
from ..surface_perturbation.execution import potential_pair,solvent_distances
from ..surface_perturbation.runtime import RunBudget
from ..util import write_json,sha256_file


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wt',type=Path,required=True);parser.add_argument('--mut',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    budget=RunBudget(args.output/'engineering_budget.json',1200)
    atoms=[read_pqr(p) for p in (args.wt,args.mut)]
    xyz=[v['xyz'] for a in atoms for v in a.values()];grid=make_grid(xyz,.5)
    points=np.asarray(xyz,float)[::10]+[4,0,0]
    report={'scope':'engineering benchmark; not a scientific CTLA4 rerun or panel validation',
        'pqr_sha256':[sha256_file(p) for p in (args.wt,args.mut)],'grid':grid}
    measured=[]
    try:
        for workers in (1,2):
            tools=PhysicalTools(args.output/f'fresh_{workers}',{'apbs_workers':workers},budget)
            started=time.perf_counter();fields=potential_pair(tools,args.wt,args.mut,grid)
            report[f'workers_{workers}']={'wall_seconds':time.perf_counter()-started,'activity':tools.activity,
                'calls':dict(tools.calls)}
            measured.append(np.stack([field(points) for field in fields]))
        error=float(np.max(np.abs(measured[0]-measured[1])))
        report['max_potential_difference_kt_e']=error
        report['numerically_equal']=bool(np.allclose(measured[0],measured[1],rtol=1e-12,atol=1e-12))
        report['apbs_pair_speedup']=report['workers_1']['wall_seconds']/report['workers_2']['wall_seconds']
        tools=PhysicalTools(args.output/'solvent_cache',{},budget)
        values=[]
        for label in ('fresh','cached'):
            started=time.perf_counter();values.append(solvent_distances(tools,atoms[0],points))
            report[f'solvent_{label}_seconds']=time.perf_counter()-started
        report['solvent_cache_exact']=all(np.array_equal(a,b,equal_nan=True) for a,b in zip(*values))
        report['status']='PASS' if report['numerically_equal'] and report['solvent_cache_exact'] else 'FAIL'
    finally:
        budget.save(False);report['elapsed_seconds']=budget.elapsed
        write_json(args.output/'acceleration_report.json',report)
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
