import csv,collections,json
from pathlib import Path
import numpy as np,uproot
H=Path(__file__).resolve().parent;R=H.parent
m=list(csv.DictReader((R/'c_mass_toy_oracle/split_manifest.tsv').open(),delimiter='\t'));n=list(csv.DictReader((R/'c_joint_beta_rigidity_closure/expanded_validation/manifest.tsv').open(),delimiter='\t'))
print('MANIFEST',len(n),dict(collections.Counter(x['isotope'] for x in n)),flush=True)
print('SPLIT',len(m),dict(collections.Counter(x['role'] for x in m)),flush=True)
p=n[-1]['path'];t=uproot.open(p)['amstreea'];print('C13',p,t.num_entries,[k for k in ['mmass','mtrmom','mtrpri','mtrpar','mtrcoo1','tof_Aadc','tof_rtime','tk_edep','tk_ql'] if k in t.keys()],flush=True)
p=H/'native_parts/048.root';a=uproot.open(str(p))['reference'].arrays(entry_stop=3,library='np')
for i in range(3):
 out={k:float(a[k][i]) for k in ['run','event','gen_mass','gen_charge','gen_particle','gen_track','gen_parent']}
 out['tof']={k:a[k][i].tolist() for k in ['tof_layer','tof_bar','tof_particle','tof_track','tof_parent','tof_time','tof_beta','tof_edep','tof_step','tof_position']}
 out['plane10']=[(int(j),float(a['momentum'][i][j]),int(a['track'][i][j]),int(a['particle'][i][j]),int(a['parent'][i][j]),a['position'][i][3*j:3*j+3].tolist()) for j in range(len(a['plane'][i])) if a['plane'][i][j]==10]
 print('NATIVE',json.dumps(out),flush=True)
