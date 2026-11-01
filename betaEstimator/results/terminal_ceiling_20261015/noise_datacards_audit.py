#!/usr/bin/env python3
"""Read only generator metadata and retain source-level timing-noise evidence."""
from pathlib import Path
import json,hashlib,uproot
from audit_common import HERE,write_json
plan=json.loads((HERE/'native_plan.json').read_text());records=[]
for r in [plan['records'][0],plan['records'][-1]]:
 f=uproot.open(r['path']);text=str(f['DataCards']);lines=text.splitlines();chosen=[]
 for j,line in enumerate(lines):
  if any(k in line.upper() for k in ['TFMC','TOF','G4FF','GEANT','BUILD','VERSION']):chosen.append({'line':j+1,'text':line})
 records.append({'source':r['source'],'raw_path':r['path'],'cards_bytes':len(text.encode()),'cards_sha256':hashlib.sha256(text.encode()).hexdigest(),'relevant_cards':chosen})
 print('CARDS',records[-1],flush=True)
root=Path('/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/CC')
evidence=[]
for name,lo,hi in [('Tofsim02.C',474,485),('Tofsim02.C',672,683),('tofdbc02.C',104,124),('job.C',976,987),('Tofsim02.C',155,166),('Tofsim02.C',502,509)]:
 p=root/name;ls=p.read_text(errors='replace').splitlines();evidence.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'start_line':lo,'end_line':hi,'text':'\n'.join(ls[lo-1:hi])})
write_json(HERE/'noise_datacards_audit.json',{'raw_metadata':records,'source_evidence':evidence,'conclusion':'Local source has additive LT discriminator jitter with amplitude daqpwd(10)*dcnoise, but local default dcnoise is zero. The active source-build/default equivalence to B1308/B1315 is not established by these cards. Photon transit-time and waveform noise remain. No independent positive Gaussian timing floor is certified.'})
print('NOISE_SOURCE_AUDIT_COMPLETE',flush=True)
