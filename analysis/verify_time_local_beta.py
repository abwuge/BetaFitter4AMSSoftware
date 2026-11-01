#!/usr/bin/env python3
"""Independently check time exclusions, measurement lineage and paired metrics."""
import json
from pathlib import Path
import sys
import numpy as np
import time_local_beta as t


def main():
    out = t.OUT
    plan = json.loads((out/'plan.json').read_text())
    result = json.loads((out/'results.json').read_text())
    config = json.loads((out/'configuration.json').read_text())
    t.w.check_freeze()
    assert t.w.ev.digest(t.w.DEST/'freeze.json') == config['source_freeze_sha256']
    assert t.w.ev.digest(t.OLD/'predictions.npz') == config['frozen_prediction_sha256']
    assert t.w.ev.digest(out/'training.npz') == plan['training_sha256']
    assert t.w.ev.digest(out/'predictions.npz') == result['prediction_sha256']
    a = t.load(out/'training.npz')
    p = t.load(out/'predictions.npz')
    original = t.load(t.OLD/'predictions.npz')
    assert set(a) == t.tr.TRAIN_KEYS
    assert all(np.allclose(original[k],p[k],rtol=0,atol=0,equal_nan=True) for k in original)
    for data in (a,p):
        ids=np.rec.fromarrays([data['run'],data['event']],names='run,event')
        assert len(np.unique(ids)) == len(ids)
    mask=p['reference_mask']; raw=p['raw_beta_rich']; new=p['local_beta_rich']; ref=p['rich_beta']
    expected=(p['run'].astype(np.int64)-plan['anchor'])//t.PERIOD
    assert np.array_equal(expected,p['fold'])
    assert len(np.unique(p['fold'])) == len(plan['folds'])
    assert np.isfinite(new).all() and (new>0).all()
    fhash=json.loads((out/'freeze.json').read_text())['folds']
    boundary_count_differences=[]
    for f in plan['folds']:
        dest=out/'folds'/('%04d'%f['fold'])
        freeze=json.loads((dest/'freeze.json').read_text())
        assert t.w.ev.digest(dest/'freeze.json') == fhash[str(f['fold'])]
        target=p['run'][p['fold']==f['fold']]
        assert np.all((target>=f['start']) & (target<f['stop']))
        union=[]
        fitted=json.loads((dest/'training.json').read_text())['blocks']
        for b,edges in zip(f['blocks'],t.w.BLOCKS):
            time_selected=t.training_mask(a['run'],f['start'],f['stop'],b['training_days'])
            # The immutable plan used np.histogram, whose last bin includes the
            # upper edge. Fits use half-open blocks. Check both independently;
            # retain and disclose any exact-edge difference instead of changing
            # the frozen plan or weakening the actual fitted-sample check.
            planned_count=np.histogram(a['R'][time_selected],edges)[0]
            assert np.array_equal(planned_count,b['train_bins'])
            assert int(planned_count.sum())==b['train_n']
            selected=time_selected & (a['R']>=edges[0]) & (a['R']<edges[-1])
            count=np.histogram(a['R'][selected],edges)[0]
            actual=fitted[b['block']]
            assert np.array_equal(count,actual['counts']) and count.min()>=t.MIN_BIN
            assert int(selected.sum())==actual['n']
            right_edge_n=int((time_selected & (a['R']==edges[-1])).sum())
            expected_difference=np.zeros(len(count),dtype=int)
            expected_difference[-1]=right_edge_n
            assert np.array_equal(planned_count-count,expected_difference)
            if right_edge_n:
                boundary_count_differences.append(dict(fold=f['fold'],block=b['block'],
                    upper_edge=edges[-1],planned_n=b['train_n'],actual_n=actual['n'],
                    right_edge_n=right_edge_n))
            runs=np.unique(a['run'][selected]); assert np.intersect1d(runs,target).size==0
            union.extend(runs.tolist())
            for candidate in b['choices'][:-1]:
                assert min(candidate['bins'])<t.MIN_BIN
            model=dest/('block_%d_model.json'%b['block'])
            assert t.w.ev.digest(model)==freeze['model_sha256'][str(b['block'])]
        saved=t.load(dest/'training_runs.npz')['run']
        assert np.intersect1d(saved,target).size==0
        assert np.isin(np.unique(union),saved).all()
    def compare(select,record):
        actual=t.w.ev.compare(raw[select],new[select],ref[select])
        for k in ('gain_beta_pct','gain_u_pct'):
            assert np.isclose(actual[k],record[k],rtol=0,atol=1e-10),(k,actual[k],record[k])
        for k in ('baseline_beta','new_beta','baseline_u','new_u'):
            for stat,value in actual[k].items():
                if isinstance(value,(int,float)):
                    assert np.isclose(value,record[k][stat],rtol=0,atol=1e-10,equal_nan=True),(k,stat)
    compare(mask,result['overall'])
    for r in result['period_rows']:
        select=mask & (p['fold']==r['fold'])
        assert select.sum()==r['reference_n']
        if select.sum()>=2: compare(select,r['local'])
    years=p['run'].astype('datetime64[s]').astype('datetime64[Y]').astype(int)+1970
    for r in result['year_rows']: compare(mask & (years==r['year']),r['local'])
    for axis,values in [('rigidity',p['R']),('beta',ref)]:
        total=0
        for r in result['rows'][axis]:
            select=mask & (values>=r['low']) & (values<r['high']);total+=select.sum()
            assert select.sum()==r['n']
            if select.sum()>=2:compare(select,r)
        assert total==mask.sum()
    # Replay a deterministic cross-time sample using full-precision measured
    # joined features. Quantify the existing diagnostic-cache float32 effect.
    rows=np.unique(np.linspace(0,len(p['run'])-1,60).astype(int))
    max_difference=0.
    for source in np.unique(p['source'][rows]):
        selected=rows[p['source'][rows]==source]
        original_a=t.load(t.OLD/'parts'/('%05d'%source)/'joined.npz')
        target_ids=np.rec.fromarrays([p['run'][selected],p['event'][selected]],names='run,event')
        ids=np.rec.fromarrays([original_a['run'],original_a['event']],names='run,event')
        order=np.argsort(ids);index=order[np.searchsorted(ids[order],target_ids)]
        assert np.array_equal(ids[index],target_ids)
        measured=t.w.ev.sub(original_a,index)
        t.w.set_edges(t.w.R_EDGES);d=t.w.rt.center_features(measured)
        for j,row in enumerate(selected):
            fold=int(p['fold'][row]); block=next(i for i,e in enumerate(t.w.BLOCKS) if e[0]<=d['R'][j]<e[-1])
            edges=t.w.BLOCKS[block]; t.w.set_edges(edges)
            local={k:v[j:j+1] for k,v in d.items()};local['bin']=np.searchsorted(edges,local['R'],side='right')-1
            model=t.w.rt.arrays(json.loads((out/'folds'/('%04d'%fold)/('block_%d_model.json'%block)).read_text()))
            x=t.tr.feature_matrix(local,model)
            corrected=local['u'][0]-x.dot(model['parameters'][:x.shape[1]])[0]
            beta=1/(corrected+p['delay_u_C12'][row])
            max_difference=max(max_difference,abs(beta-new[row]))
    assert max_difference<1e-6, max_difference
    for name in ('time_gain','rigidity_gain','beta_gain'):
        assert (out/'figures'/(name+'.png')).exists()
    histogram_path=out/'distributions/histograms.json'
    distribution_verified=False
    if histogram_path.exists():
        hist=json.loads(histogram_path.read_text())
        for axis,values in [('rigidity',p['R']),('beta',ref)]:
            rows=hist['axes'][axis]
            assert len(rows)==len(result['rows'][axis])
            for row in rows:
                select=mask & (values>=row['low']) & (values<row['high'])
                assert int(select.sum())==row['n']
                assert hist['methods']==['raw','local']
                samples=[raw[select]-ref[select],new[select]-ref[select]]
                for view,detail in row['views'].items():
                    for i,sample in enumerate(samples):
                        h,_=np.histogram(sample,np.array(detail['edges']))
                        assert np.array_equal(h,np.array(detail['counts'][i]))
                        assert int(len(sample)-h.sum())==detail['outside'][i]
                    if view=='full': assert detail['outside']==[0,0]
            import re
            for view in ('core','full'):
                pdf=out/'distributions'/(axis+'_distributions_'+view+'.pdf')
                assert len(re.findall(rb'/Type\s*/Page\b',pdf.read_bytes()))==len(rows)
        distribution_verified=True
    t.w.ev.dump(out/'verification.json',dict(passed=True,folds=len(plan['folds']),
        reference_n=int(mask.sum()),train_n=len(a['run']),target_run_overlap=0,
        each_event_predicted_once=True,original_frozen_outputs_unchanged=True,
        RICH_fields_in_training=False,RICH_selection_for_training=False,
        full_precision_replay_n=60,max_float32_beta_difference=max_difference,
        distributions_verified=distribution_verified,
        planned_upper_edge_count_differences=boundary_count_differences,
        actual_fitted_counts_verified=True,
        prediction_sha256=t.w.ev.digest(out/'predictions.npz')))
    print('VERIFIED',len(plan['folds']),'folds;',int(mask.sum()),'reference events; float32 max beta difference',max_difference)


if __name__=='__main__':main()
