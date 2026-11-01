import json
from pathlib import Path
import numpy as np,uproot
from analyze import metric
H=Path(__file__).resolve().parent;a=np.load(str(H/'aligned.npz'));z=np.load(str(H/'decomposition.npz'))
with uproot.open(str(H/'intervals.root')) as f:k=f['intervals'].arrays(library='np')
assert np.array_equal(k['event_index'],a['ref_event_index'])
tr=(a['ref_role']==b'train')&z['validtruth'];te=(a['ref_role']==b'test')&z['validtruth'];err=np.log(a['ref_R']/a['ref_R_true']);result={'intervals':[]}
for l in range(1,7):
 x=k['kick'][:,l,1];b=k['bend'][:,l,1]
 result['intervals'].append({'from_layer':l+1,'to_layer':l+2,'kick_y_mrad':metric(x[te]*1000),'median_magnetic_bend_y_mrad':float(np.median(b[te])*1000),'transport_y_um':metric(k['offset'][te,l,1]*1e4),'dlogp':metric(k['dlogp'][te,l]),'correlation_with_rawR':float(np.corrcoef(x[te]*a['ref_R'][te],err[te])[0,1])})
features=np.concatenate([k['kick'][:,1:7,1]*a['ref_R'][:,None],k['offset'][:,1:7,1]*a['ref_R'][:,None]],axis=1)
mu=features[tr].mean(axis=0);sd=features[tr].std(axis=0);X=(features-mu)/sd
coef=np.linalg.lstsq(X[tr],err[tr]-np.mean(err[tr]),rcond=None)[0];pred=X.dot(coef)
result['truth_innovations_diagnostic']={'test_after_removal':metric((err-pred)[te]),'test_correlation':float(np.corrcoef(pred[te],err[te])[0,1]),'coefficients':coef.tolist(),'feature_labels':['kick_L%d_L%d'%(i,i+1) for i in range(2,8)]+['offset_L%d_L%d'%(i,i+1) for i in range(2,8)],'warning':'Truth state innovations are diagnostic; unavailable in real data.'}
(H/'interval_results.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
