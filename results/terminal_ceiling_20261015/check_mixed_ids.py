import csv,numpy as np
from pathlib import Path
r=Path(__file__).resolve().parent.parent
m=list(csv.DictReader((r/'c_joint_beta_rigidity_closure/expanded_validation/manifest.tsv').open(),delimiter='\t'))
s=list(csv.DictReader((r/'c_mass_toy_oracle/split_manifest.tsv').open(),delimiter='\t'))
p=np.load(str(r/'c_joint_beta_rigidity_closure/expanded_validation/predictions.npz'))
print('MANIFEST',*[m[i] for i in [0,255,256,-1]],sep='\n')
print('SOURCES',[(j,np.unique(p['source'][p['label']==j])[:8].tolist()) for j in [0,1]])
print('SOURCE0',np.unique(p['label'][p['source']==0],return_counts=True))
print('SPLIT13',[x for x in s if x['isotope']=='C13'][:2])
