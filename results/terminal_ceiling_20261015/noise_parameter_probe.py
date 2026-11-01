#!/usr/bin/env python3
from pathlib import Path
import re,json,hashlib
root=Path('/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9')
names=['CC/tofdbc02.C','CC/Tofdbc.C','CC/job.C','CC/amsvmc.C','include/tofdbc02.h','include/link.h','include/Tofdbc.h']
pat=re.compile(r'dcnoise|daqpwd|_daqpwd|DAQPWD|TimeSigma\s*=')
records=[]
for name in names:
 p=root/name
 if not p.exists():continue
 lines=p.read_text(errors='replace').splitlines();hits=[i for i,l in enumerate(lines) if pat.search(l)]
 for i in hits[:20]:
  lo=max(0,i-3);hi=min(len(lines),i+8);record={'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'start_line':lo+1,'end_line':hi,'text':'\n'.join(lines[lo:hi])};records.append(record)
  print(name,lo+1,record['text'],flush=True)
p=Path(__file__).resolve().parent/'noise_parameter_probe.json';p.write_text(json.dumps(records,indent=2))
import uproot
raw='/lustre03/data/eos/ams/MC/AMS02/2022/C.B1308/c12.pl1.l1.612000.6_05/82334649.00000001.root'
f=uproot.open(raw);print('RAW_KEYS',f.classnames(),flush=True)
