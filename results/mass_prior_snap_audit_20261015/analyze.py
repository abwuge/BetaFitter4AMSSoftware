#!/usr/bin/env python3
"""MC holdout audit: measured beta vs mass-peak projection without/with isotope labels.
All computations stay on sdiat; only compact JSON and TSV are saved.
"""
import json
import numpy as np
from pathlib import Path
def auc(label, score):
    order = np.argsort(score,kind='stable')
    values=score[order]
    positives=(label[order]==1).astype(float)
    negatives=1-positives
    starts=np.r_[0, np.flatnonzero(np.diff(values)!=0)+1]
    pg=np.add.reduceat(positives,starts)
    ng=np.add.reduceat(negatives,starts)
    return float(np.dot(pg,np.cumsum(ng)-0.5*ng)/(sum(pg)*sum(ng)))

ROOT = Path('/nas01/home/hxwu/AMSSoft/betaEstimator/results')
OUT = ROOT / 'mass_prior_snap_audit_20261015'
P = np.load(ROOT/'c_joint_beta_rigidity_closure/expanded_validation/predictions.npz')
E = np.load(ROOT/'c_joint_beta_rigidity_closure/expanded_validation/energy_comparison/energy_predictions.npz')
M = np.load(ROOT/'terminal_ceiling_20261015/mixed.npz')

assert np.array_equal(P['run'], E['run']) and np.array_equal(P['event'], E['event'])
assert np.array_equal(P['run'], M['run']) and np.array_equal(P['event'], M['event'])
assert np.array_equal(P['label'], M['label'])
# P and M use different source-index codings; run/event/label are exact.
sel = P['truth_valid'] & np.isfinite(E['u_anchored_energy'])
lab = P['label'][sel]
R = P['R'][sel]
Rtrue = P['R_true'][sel]
ut = P['u_true'][sel]
beta_true = 1/ut
u_raw = P['u_raw'][sel]
u_energy = E['u_anchored_energy'][sel]
mass_true = M['true_mass'][sel]
sources = M['source'][sel]
m12, m13 = [float(np.median(mass_true[lab==k])) for k in (0,1)]
mass_by_label = np.where(lab==0,m12,m13)
assert np.allclose(mass_by_label,mass_true,rtol=1e-7)
assert np.all(np.isfinite(R)) and np.all(R>0)
assert np.all(np.isfinite(ut)) and np.all(ut>=1)

def beta_from_mass(R, mass):
    pr = 6*R
    return pr/np.hypot(pr,mass)

def w68(x, weights=None):
    if weights is None:
        return float(np.diff(np.quantile(x,[0.16,0.84]))[0]/2)
    order = np.argsort(x)
    y=x[order]; w=weights[order]
    return float(np.diff(np.interp([.16,.84],(np.cumsum(w)-.5*w)/sum(w),y))[0]/2)

def prevalence_weights(pi):
    counts = np.bincount(lab,minlength=2)
    return np.where(lab==0, (1-pi)/counts[0], pi/counts[1])

def metrics(b, prior=None, score_override=None):
    d = b/beta_true-1
    u = 1/b
    mass2=(6*R)**2*(u*u-1) if score_override is None else score_override
    res = {}
    for k,name in [(0,'C12'),(1,'C13')]:
        mask=lab==k
        res[name]={'n':int(mask.sum()), 'median_dbeta_relative_pct':float(100*np.median(d[mask])),
            'w68_dbeta_relative_pct':100*w68(d[mask]),
            'rmse_dbeta_relative_pct':float(100*np.sqrt(np.mean(d[mask]**2))),
            'w68_du':w68(u[mask]-ut[mask]),
            'w68_mass2_GeV2':w68(mass2[mask])}
    for p in [0.1,0.5]:
        w=prevalence_weights(p)
        res['mix_pi13_'+str(p)]={'w68_dbeta_relative_pct':100*w68(d,w),
            'rmse_dbeta_relative_pct':float(100*np.sqrt(np.average(d*d,weights=w))),
            'mean_dbeta_relative_pct':float(100*np.average(d,weights=w)),
            'w68_du':w68(u-ut,w)}
    res['mass2_auc']=float(auc(lab,mass2))
    if prior is not None:
        # Actual inferred mass from the peak projection; use equality to classify.
        is_c13=np.isclose(mass2,m13*m13,rtol=1e-9)
        is_c12=np.isclose(mass2,m12*m12,rtol=1e-9)
        if np.all(is_c13 | is_c12):
            tpr=float(np.mean(is_c13[lab==1]))
            fpr=float(np.mean(is_c13[lab==0]))
            res['hard_classifier']={'tpr_c13':tpr,'fpr_c12':fpr,
              'balanced_accuracy':(tpr+1-fpr)/2,
              'accuracy_at_pi13_'+str(prior):(1-prior)*(1-fpr)+prior*tpr,
              'predicted_c13_fraction_for_pi13_'+str(prior):(1-prior)*fpr+prior*tpr}
    return res

x=(6*R)**2*(u_raw*u_raw-1)
mid=(m12*m12+m13*m13)/2
gap=m13*m13-m12*m12
sig=35.0   # fixed illustrative response width GeV^2, not tuned to the holdout
out={'N':len(lab),'sources':int(len(np.unique(sources))),
  'counts':{'C12':int((lab==0).sum()),'C13':int((lab==1).sum())},
  'mass_GeV':{'C12':m12,'C13':m13},
  'settings':{'sigma_m2_GeV2':sig,'pi13_scenarios':[.1,.5],
              'note':'pi13=0.10 is illustrative, not inferred from these MC events; Gaussian likelihood is a deliberately simple constructed prior test.'}}
methods={
 'raw_tof':1/u_raw,
 'anchored_energy':1/u_energy,
 'oracle_truth_label_plus_measured_R':beta_from_mass(R,mass_true),
 'oracle_truth_label_plus_truth_R':beta_from_mass(Rtrue,mass_true),
 'all_C12_mass_plus_measured_R':beta_from_mass(R,m12),
 'all_C13_mass_plus_measured_R':beta_from_mass(R,m13),
 'prior_weighted_mass_pi13_0.10':beta_from_mass(R,.9*m12+.1*m13),
}
for pi in [.1,.5]:
    logit=np.log(pi/(1-pi)) + gap*(x-mid)/(sig*sig)
    q13=1/(1+np.exp(-np.clip(logit,-700,700)))
    hard=(q13>=.5)
    pred_mass_hard=np.where(hard,m13,m12)
    pred_mass_soft=(1-q13)*m12+q13*m13
    methods['hard_MAP_pi13_'+str(pi)]=beta_from_mass(R,pred_mass_hard)
    methods['posterior_mean_mass_pi13_'+str(pi)]=beta_from_mass(R,pred_mass_soft)
    out['settings']['threshold_m2_pi13_'+str(pi)]=float(mid + (sig*sig/gap)*np.log((1-pi)/pi))
    out['settings']['hard_c13_frac_on_natural_sample_pi13_'+str(pi)]=float(hard.mean())
out['methods']={}
scores={
 'oracle_truth_label_plus_measured_R':mass_true**2,
 'oracle_truth_label_plus_truth_R':mass_true**2,
 'all_C12_mass_plus_measured_R':np.full(len(lab),m12*m12),
 'all_C13_mass_plus_measured_R':np.full(len(lab),m13*m13),
 'prior_weighted_mass_pi13_0.10':np.full(len(lab),(.9*m12+.1*m13)**2),
}
for pi in [.1,.5]:
    logit=np.log(pi/(1-pi))+gap*(x-mid)/(sig*sig)
    post=1/(1+np.exp(-np.clip(logit,-700,700)))
    scores['hard_MAP_pi13_'+str(pi)]=np.where(post>=.5,m13*m13,m12*m12)
    scores['posterior_mean_mass_pi13_'+str(pi)]=((1-post)*m12+post*m13)**2
for name,b in methods.items():
    pi=float(name.split('_')[-1]) if name.startswith('hard_MAP_pi13_') else None
    out['methods'][name]=metrics(b,pi,scores.get(name))
out['mass_AUC_raw_equivalent']=float(auc(lab,x))
out['R_bin_w68_beta_relative_pct']={}
for lo,hi in [(2.5,3),(3,4),(4,5),(5,6)]:
    mask=(R>=lo)&(R<hi)
    segment={}
    for name in ['raw_tof','anchored_energy','oracle_truth_label_plus_measured_R',
                 'hard_MAP_pi13_0.5','posterior_mean_mass_pi13_0.5']:
        d=methods[name][mask]/beta_true[mask]-1
        segment[name]=100*w68(d)
    out['R_bin_w68_beta_relative_pct'][str(lo)+'-'+str(hi)]=segment

(OUT/'summary.json').write_text(json.dumps(out,indent=2,ensure_ascii=False))
with (OUT/'summary.tsv').open('w') as f:
    f.write('method\tC12_betaW68_pct\tC13_betaW68_pct\tmix10_betaW68_pct\tmix10_RMSE_pct\tC12_m2W68\tC13_m2W68\tAUC\tC12_median_bias_pct\tC13_median_bias_pct\n')
    for name,v in out['methods'].items():
        f.write('\t'.join(map(str,[name,v['C12']['w68_dbeta_relative_pct'],v['C13']['w68_dbeta_relative_pct'],
         v['mix_pi13_0.1']['w68_dbeta_relative_pct'],v['mix_pi13_0.1']['rmse_dbeta_relative_pct'],
         v['C12']['w68_mass2_GeV2'],v['C13']['w68_mass2_GeV2'],v['mass2_auc'],
         v['C12']['median_dbeta_relative_pct'],v['C13']['median_dbeta_relative_pct']]))+'\n')
print('N',len(lab),'m12/m13',m12,m13,'sigma_m2',sig)
for name,v in out['methods'].items():
 print('%-40s C12_W68 %7.4f C13_W68 %7.4f 10%%W68 %7.4f 10%%RMSE %7.4f  mass2_AUC %.6f  m2_W68_C12/C13 %.3f/%.3f'%
   (name,v['C12']['w68_dbeta_relative_pct'],v['C13']['w68_dbeta_relative_pct'],
    v['mix_pi13_0.1']['w68_dbeta_relative_pct'],v['mix_pi13_0.1']['rmse_dbeta_relative_pct'],
    v['mass2_auc'],v['C12']['w68_mass2_GeV2'],v['C13']['w68_mass2_GeV2']))
for p in [.1,.5]:
 n='hard_MAP_pi13_'+str(p)
 print(n,'threshold',out['settings']['threshold_m2_pi13_'+str(p)],'confusion',out['methods'][n].get('hard_classifier'))
print('BIN_METRICS',out['R_bin_w68_beta_relative_pct'])
