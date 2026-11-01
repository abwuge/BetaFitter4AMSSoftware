import json,sys
from pathlib import Path
import numpy as np
from analyze import metric,boot
H=Path(__file__).resolve().parent;sys.path.insert(0,str(H.parent))
from summarize_tkadc import load_tmva
a=np.load(str(H/'aligned.npz'));z=np.load(str(H/'decomposition.npz'));t=load_tmva('tkadc_base','TestTree');lookup={int(v):i for i,v in enumerate(a['ref_event_index'])};ix=np.array([lookup[int(v)] for v in t['event_index']]);keep=z['validtruth'][ix];ix=ix[keep]
pred=np.asarray(t['BDTG']).reshape(-1)[keep];err=np.log(a['ref_R'][ix]/a['ref_R_true'][ix]);eb=err-pred;src=a['ref_source'][ix]
r={'n':len(ix),'baseline':metric(eb),'raw':metric(err),'baseline_gain_CI':boot(err,eb,src),'transport_correlation_after_baseline':float(np.corrcoef(z['transport'][ix],eb)[0,1]),'caveat':'Truth removal and new hit correction use approximate native-field GLS response, not exact Kalman reruns.'}
for component in ['measurement','transport']:
 ec=eb-z[component][ix];r['remove_'+component]={'metric':metric(ec),'CI_percent':boot(eb,ec,src)}
for stage in ['base','adc']:
 ph=np.load(str(H/('hit_'+stage+'_prediction_cm.npy')));delta=np.sum(z['weights'][ix]*ph[ix],axis=1);ec=eb-delta;r['new_hit_'+stage]={'metric':metric(ec),'incremental_CI_percent':boot(eb,ec,src),'gain_point_percent':100*(1-metric(ec)['w68']/metric(eb)['w68'])}
r['inner_momentum_change']=metric(np.log(a['mevmom1'][ix,13]/a['mevmom1'][ix,7]))
(H/'after_baseline.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
