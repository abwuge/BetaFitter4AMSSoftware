#!/usr/bin/env python3
"""Continue ratio training through new-prediction center-plane figures."""
import json
from pathlib import Path
import subprocess
import sys
import time
import time_local_ratio as model


def run(script,*args):
    subprocess.run([sys.executable,str(model.ROOT/'analysis'/script)]+list(args),check=True)


def main():
    out=model.OUT;expected=len(json.loads((model.OLD/'plan.json').read_text())['folds'])
    last=-1
    while True:
        count=len(list(out.glob('folds/*/freeze.json')))
        if count!=last:print('FROZEN',count,'/',expected,flush=True);last=count
        failures=list(out.glob('folds/*/failure.json'))
        if failures:raise RuntimeError('Training failures: '+str(failures[:4]))
        if count==expected:break
        queue=json.loads(subprocess.check_output(['condor_q','-name','scheduler@schedd01.ams.sdiat.cn','399127','-json'],universal_newlines=True) or '[]')
        held=[j for j in queue if j['JobStatus']==5]
        if held:raise RuntimeError('Held jobs: '+str([(j['ProcId'],j.get('HoldReason')) for j in held]))
        if not queue:
            count=len(list(out.glob('folds/*/freeze.json')))
            if count==expected:break
            raise RuntimeError('Queue ended before all folds froze')
        time.sleep(15)
    run('time_local_ratio.py','evaluate')
    run('plot_time_local_ratio.py','--workers','8')
    run('plot_ratio_material_closure.py')
    final=out/'final_comparison.pdf'
    subprocess.run(['pdfunite',str(out/'figures/center_comparison.pdf'),str(out/'material_closure/material_closure.pdf'),str(final)],check=True)
    info=subprocess.check_output(['pdfinfo',str(final)],universal_newlines=True)
    if 'Pages:          16' not in info:raise RuntimeError('Unexpected final PDF page count')
    model.dump(out/'completion.json',dict(complete=True,figures_generated_after_new_model_freeze=True,
        final_pdf=str(final),final_pdf_sha256=model.local.w.ev.digest(final),pages=16,
        reference_direction='RICH to center (L5)',MC_RICH_used_for_training=False))
    print('RATIO_CENTER_FIGURES_COMPLETE',final,flush=True)


if __name__=='__main__':main()
