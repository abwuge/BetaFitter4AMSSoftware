#!/usr/bin/env python3
"""Inspect existing primary-reference planes and reconstructed-ring field names."""
import json
from pathlib import Path
import uproot
import awkward as ak
from research_flight_conditional import ROOT

base=ROOT/'results/flight_calibration_research/continuation_20261013/additional_mc'
path=sorted((base/'reference/parts').glob('*.root'))[0]
with uproot.open(path) as stream:
    row=ak.to_list(stream['reference'].arrays(['plane','position','momentum','mass','particle','track','parent','gen_particle','gen_track','gen_parent'],entry_stop=1,library='ak'))[0]
print('REFERENCE_FILE',str(path),flush=True)
print('PLANES',[(row['plane'][i],row['position'][3*i+2],row['momentum'][i],row['track'][i],row['parent'][i]) for i in range(len(row['plane']))],flush=True)
print('GENERATOR',row['gen_track'],row['gen_parent'],flush=True)
source=json.loads((base/'protocol.json').read_text())['records'][0]['path']
with uproot.open(source) as stream:
    print('RICH_FIELDS',[key for key in stream['amstreea'].keys() if key.startswith('rich_')],flush=True)
