from pathlib import Path
import numpy as np,uproot,json
root=Path('/nas01/home/hxwu/AMSSoft/betaEstimator')
rows=[l.split() for l in (root.parent/'betaFitter/config/hybrid_material_fine.sources').read_text().splitlines() if l and not l.startswith('#')]
for z in [1,2,6,8]:
 p=next(r[2] for r in rows if int(r[1])==z)
 t=uproot.open(p)['amstreea'];ks=list(t.keys())
 print('FILE',z,p,'entries',t.num_entries,flush=True)
 print('TRUTH',[(k,t[k].typename) for k in ks if k.startswith(('mtof','mev','mtr','tof_betahmc'))],flush=True)
 if z==6:
  f=['run','event','mpar','mch','mevmom1','mevcoo1','mtrmom','mtrpri','mtrpar','mtrcoo','mtrcoo1','mtrz','tof_tl','tof_leng','tof_pos','mtofdep','mtoftime','mtofcoo','tof_betahmc','mmass']
  a=t.arrays([k for k in f if k in ks],entry_stop=2,library='np')
  print('EXAMPLE',json.dumps({k:v.tolist() for k,v in a.items()}),flush=True)
 p0=root/'results/comparison/comparison_Z%d.root'%z if False else root/('results/comparison/comparison_Z%d.root'%z)
 q=uproot.open(str(p0));print('COMPARISON',q.keys(),flush=True)
 for k in q.keys():
  tr=q[k]
  if hasattr(tr,'num_entries'):print(k,tr.num_entries,list(tr.keys()),flush=True)
