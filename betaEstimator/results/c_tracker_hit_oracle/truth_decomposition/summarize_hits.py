import json
from pathlib import Path
import numpy as np,uproot
from analyze import metric,boot
H=Path(__file__).resolve().parent;out={};allpred={}
for stage in ['base','adc']:
 with uproot.open(str(H/('hit_'+stage+'.root'))) as f:
  t=f['ds_hit_'+stage+'/TestTree'].arrays(library='np')
  tr=f['ds_hit_'+stage+'/TrainTree'].arrays(library='np')
 y=t['true_error_um'];pred=np.asarray(t['BDTG']).reshape(-1);va=t['role']==1;te=t['role']==2;src=t['source'][te]
 scan=[metric((y-alpha*pred)[va])['w68'] for alpha in [0,.25,.5,.75,1.]];alpha=[0,.25,.5,.75,1.][int(np.argmin(scan))];res=y-alpha*pred
 out[stage]={'chosen_alpha':alpha,'validation_w68_scan':scan,'test_original_um':metric(y[te]),'test_corrected_um':metric(res[te]),'test_gain_CI_percent':boot(y[te],res[te],src),'training_original_um':metric(tr['true_error_um']),'training_corrected_um':metric(tr['true_error_um']-alpha*np.asarray(tr['BDTG']).reshape(-1)),'test_correlation':float(np.corrcoef(y[te],pred[te])[0,1])}
 allpred[stage]=(t['hit_index'].astype(int),alpha*pred,res,t['role'],t['source'])
bi,bp,br,roles,src=allpred['base'];ai,ap,ar,aroles,asrc=allpred['adc'];assert np.array_equal(bi,ai)
te=roles==2;out['incremental_adc_vs_base_CI_percent']=boot(br[te],ar[te],src[te]);out['incremental_adc_vs_base_percent']=100*(1-metric(ar[te])['w68']/metric(br[te])['w68'])
a=np.load(str(H/'aligned.npz'));z=np.load(str(H/'decomposition.npz'));R=a['ref_R'];rt=a['ref_R_true'];test=(a['ref_role']==b'test')&z['validtruth'];valid=(a['ref_role']==b'valid')&z['validtruth'];err=np.log(R/rt)
# Response propagation only: this is not a rerun of the original Kalman fit.
for stage in ['base','adc']:
 ix,pred,_,_,_=allpred[stage];P=np.zeros(len(R)*7);P[ix]=pred;P=P.reshape(-1,7)*1e-4
 delta=np.sum(z['weights']*P,axis=1);corrected=err-delta
 out[stage]['linear_response_test']=metric(corrected[test]);out[stage]['linear_response_gain_CI_percent']=boot(err[test],corrected[test],a['ref_source'][test]);out[stage]['linear_response_note']='Measured-input native-field GLS response applied to baseline; not an exact Kalman replay.'
 np.save(str(H/('hit_'+stage+'_prediction_cm.npy')),P)
(H/'hit_results.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
