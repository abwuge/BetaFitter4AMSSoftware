#!/usr/bin/env python3
"""Check the full NEWTOF population, isolation, numerical bins and figures."""
import json
import re
from pathlib import Path

import numpy as np
import flight_beta_wide as w
from prepare_tofcalib_validation import OUT


def read(path):
    return json.loads(path.read_text())


def main():
    freeze=w.check_freeze()
    config=read(OUT/'selection.json')
    notes=[]
    binary=w.ev.digest(OUT/'bin/supplement_light')
    for source in config['file_indices']:
        part=OUT/'parts'/('%05d'%source)
        note=read(part/'supplement.json')
        assert note['complete'] and note['source']==source
        if note['requested']:
            assert (part/'.native_complete').read_text().strip()==binary
            assert note['binary_sha256']==binary
            assert not any(note['association_mismatches'].values())
        notes.append(note)
    progress=read(OUT/'supplement_progress.json')
    assert progress['complete'] and progress['completed_files']==len(notes)==4096
    assert sum(n['requested'] for n in notes)==read(OUT/'native_preparation.json')['reference_preselected']
    assert sum(n['selected'] for n in notes)==progress['selected']
    proof=read(OUT/'probe/batch_equivalence.json')
    assert proof['required_measurements_equal'] and proof['light_binary_sha256']==binary
    assert read(OUT/'probe/el9_equivalence.json')['passed']
    training=read(OUT/'training_runs.json')
    assert training['training_cache_sha256']==freeze['train_sha256']
    with np.load(str(OUT/'training_runs.npz'),allow_pickle=False) as f:train_runs=f['run']
    assert len(train_runs)==training['actual_training_run_count']==341
    assert w.ev.digest(OUT/'training_runs.npz')==training['run_set_sha256']
    cache=w.CACHE/'train.npz'
    if cache.exists():
        assert w.ev.digest(cache)==freeze['train_sha256']
        with np.load(str(cache),allow_pickle=False) as f:
            assert np.array_equal(np.unique(f['run']),train_runs)
    path=OUT/'predictions.npz';digest=w.ev.digest(path)
    audit=read(OUT/'prediction_audit.json')
    assert audit['complete'] and audit['prediction_sha256']==digest
    assert audit['model_freeze_sha256']==w.ev.digest(w.DEST/'freeze.json')
    with np.load(str(path),allow_pickle=False) as f:p={k:f[k] for k in f.files}
    ids=np.rec.fromarrays([p['run'],p['event']],names='run,event')
    assert len(np.unique(ids))==len(ids)==audit['selected']
    assert not np.intersect1d(p['run'],train_runs).size
    assert len(ids)+audit['identical_duplicates_removed']==progress['selected']
    result=read(OUT/'results.json');axis_result=read(OUT/'axis_results.json')
    hist=read(OUT/'distributions/histograms.json')
    assert result['complete'] and not result['retrained'] and not result['RICH_used_for_training']
    assert hist['prediction_sha256']==digest
    mask=p['reference_mask'];n=int(mask.sum())
    assert n==audit['reference_n']==result['reference_n']
    assert len(result['rows']['rigidity'])==48 and len(result['rows']['beta'])==29
    distributions={(x['axis'],x['low'],x['high']):x for x in hist['records']}
    for axis,rows in result['rows'].items():
        values=p['R'] if axis=='rigidity' else p['rich_beta']
        assert sum(x['n'] for x in rows)==n
        for row,plot_row in zip(rows,axis_result['rows'][axis]):
            select=mask&(values>=row['low'])&(values<row['high'])
            assert int(select.sum())==row['n']==plot_row['n']
            assert (row['low'],row['high'])==(plot_row['low'],plot_row['high'])
            if row['n']>=2:
                metrics=w.ev.compare(p['raw_beta_rich'][select],p['new_beta_rich'][select],p['rich_beta'][select])
                for method in ('baseline_beta','new_beta'):
                    for metric in ('w68','mean','rmse','p95_abs'):
                        assert np.isclose(metrics[method][metric],row[method][metric],rtol=0,atol=1e-14)
            if row['n']>=20:
                assert row['bootstrap']['repeats']==400
                expected=np.array(row['bootstrap']['beta_gain_pct_q025_q50_q975'])[[0,2]]
                assert np.allclose(expected,plot_row['errors']['gain_95ci'],rtol=0,atol=1e-9)
            entry=distributions[(axis,row['low'],row['high'])]
            assert entry['n']==row['n']
            for view in ('core','full'):
                edges=hist['core_edges'] if view=='core' else entry['full_edges']
                for counts,key in zip(entry[view],('raw_beta_rich','new_beta_rich')):
                    residual=(p[key][select]-p['rich_beta'][select])*1000
                    actual,_=np.histogram(residual,bins=edges)
                    assert np.array_equal(actual,counts['histogram'])
                    outside=int(((residual<edges[0])|(residual>edges[-1])).sum())
                    assert outside==counts['outside_view']
                    assert sum(counts['histogram'])+outside==row['n']
                    if view=='full':assert outside==0
        for suffix in ('','_full'):
            pdf=OUT/'distributions'/(axis+'_residuals'+suffix+'.pdf')
            pages=len(re.findall(rb'/Type\s*/Page\b',pdf.read_bytes()))
            assert pages==len(rows),(pdf,pages,len(rows))
        for metric in ('width','bias','gain'):
            for extension in ('pdf','png'):
                assert (OUT/'figures'/(axis+'_'+metric+'.'+extension)).stat().st_size>1000
    assert len(distributions)==77
    assert result['retained_raw_beta_above_one']==int((mask&(p['raw_beta_rich']>1)).sum())
    assert result['retained_new_beta_above_one']==int((mask&(p['new_beta_rich']>1)).sum())
    broad=[]
    for low,high in [(2.5,6),(6,20),(20,200)]:
        select=mask&(p['R']>=low)&(p['R']<high)
        broad.append(dict(low=low,high=high,n=int(select.sum()),**w.ev.compare(
            p['raw_beta_rich'][select],p['new_beta_rich'][select],p['rich_beta'][select])))
    assert w.ev.digest(path)==digest
    w.ev.dump(OUT/'verification.json',dict(passed=True,source_files=len(notes),reference_n=n,
        unique_events=True,actual_run_count=len(np.unique(p['run'])),training_run_overlap=0,
        frozen_models_unchanged=True,rigidity_bins=48,beta_bins=29,
        all_bin_distributions_verified=True,full_range_pdf_pages=77,core_pdf_pages=77,
        beta_above_one_retained=True,prediction_sha256=digest,broad_rigidity_comparison=broad,
        script_sha256=w.ev.digest(Path(__file__))))
    print('Verified:',n,'independent reference events; 77 bins and all distributions',flush=True)


if __name__=='__main__':main()
