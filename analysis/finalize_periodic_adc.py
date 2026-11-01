#!/usr/bin/env python3
"""Freeze each completed two-start experiment without overwriting prior outputs."""
import json
import subprocess
import sys
from research_periodic_adc import ROOT, OUT, CONFIGS
from research_flight_conditional import write,digest
from research_conditional_mc import stamp


def main():
    rows=[]
    for dataset in ['flight','mc']:
        for family,_,_ in CONFIGS[dataset]:
            for variant in ['static','step','spline']:
                folder=OUT/dataset/(family+'_'+variant)
                starts=[folder/('start_%d.json'%i) for i in range(2)]
                present=[p.exists() for p in starts]
                row=dict(dataset=dataset,family=family,variant=variant,completed_starts=sum(present))
                if all(present):
                    records=[json.loads(p.read_text()) for p in starts]
                    best=min(records,key=lambda item:item['penalized_nll'])
                    row.update(converged=[r['numerically_converged'] for r in records],
                        best_nll=best['penalized_nll'],best_gradient=best['projected_gradient_max'])
                    if best['numerically_converged'] and not (folder/'freeze.json').exists():
                        log=folder/'freeze.log'
                        if log.exists():raise FileExistsError('An earlier freeze attempt needs inspection: '+str(folder))
                        with log.open('w') as stream:
                            result=subprocess.run([sys.executable,'analysis/research_periodic_adc.py','freeze','--dataset',dataset,
                                '--family',family,'--variant',variant],cwd=str(ROOT),stdout=stream,stderr=subprocess.STDOUT)
                        if result.returncode:raise RuntimeError('Period freeze failed: '+str(log))
                row['frozen']=(folder/'measurement_evaluation.json').exists()
                rows.append(row)
    print('PERIOD_STATUS',json.dumps(rows),flush=True)
    complete=all(row['frozen'] for row in rows)
    if complete and not (OUT/'fit_completion.json').exists():
        write(OUT/'fit_completion.json',dict(complete=True,created_utc=stamp(),models=rows,
            total_starts=sum(row['completed_starts'] for row in rows),
            freeze_sha256={str(p.relative_to(ROOT)):digest(p) for p in OUT.glob('*/*/freeze.json')},
            source_sha256=digest(__file__)))
        print('PERIOD_ALL_MODELS_FROZEN',len(rows),flush=True)


if __name__=='__main__':main()
