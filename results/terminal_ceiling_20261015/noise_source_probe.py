#!/usr/bin/env python3
"""Bounded local-source probe of the timing simulation."""
from pathlib import Path
import re
root=Path('/nas01/home/hxwu/AMSSoft/AMSVX_TOFBetaFit_el9/CC')
paths=sorted(p for p in root.iterdir() if p.is_file() and p.suffix in ['.C','.cxx','.cc'] and 'tof' in p.name.lower())
print('FILES',[(p.name,p.stat().st_size) for p in paths],flush=True)
pat=re.compile(r'gaus|rannor|smear|timres|tim.*sigma|sigma.*tim|rnorm',re.I)
for p in paths:
    lines=p.read_text(errors='replace').splitlines();hits=[i for i,l in enumerate(lines) if pat.search(l)]
    if not hits:continue
    print('\nFILE',p.name,'HITS',len(hits),flush=True)
    for i in hits[:24]:
        print('\n'.join('%d:%s'%(j+1,lines[j]) for j in range(max(0,i-2),min(len(lines),i+3))),flush=True)
