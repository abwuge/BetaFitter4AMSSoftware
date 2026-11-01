#!/usr/bin/env python3
"""Record a field-name correction without changing the original MC cohort."""
import json
from pathlib import Path
import numpy as np
from audit_mc_rich_charge import OUT,POOL
from research_flight_conditional import ROOT,digest,load,write
from research_conditional_mc import stamp


def main():
    if (OUT/'scoring_amendment.json').exists():raise FileExistsError('MC scoring amendment already exists')
    original=json.loads((OUT/'protocol.json').read_text());manifest=json.loads((OUT/'readout_manifest.json').read_text())
    a=load(POOL/'measured.npz');chi=a['chi']
    valid=np.isfinite(chi).all(1)&((chi>0)&(chi<=10)).all(1)
    if not valid.all():raise ValueError('The MC cache does not already satisfy the documented chi-square domain')
    if digest(OUT/'readout.npz')!=manifest['readout_sha256']:raise ValueError('MC ring readout changed')
    write(OUT/'scoring_amendment.json',dict(created_utc=stamp(),original_source_sha256=original['source_sha256'],
        corrected_source_sha256=digest(ROOT/'analysis/audit_mc_rich_charge.py'),protocol_sha256=digest(OUT/'protocol.json'),
        measured_sha256=digest(POOL/'measured.npz'),readout_sha256=manifest['readout_sha256'],
        preparation_source_sha256=digest(ROOT/'analysis/research_conditional_mc.py'),
        change='Use a[chi][:,2], the stored tk_chis1[:,2,2,1,2] component, instead of the absent flight-only key native_chi.',
        events=len(valid),already_pass_all_three_chi_components=int(valid.sum()),membership_changed=False,
        prior_attempt='All 128 source and native-reference reads completed and were frozen. Scoring then raised KeyError before computing its beta diagnostics.',
        original_protocol_preserved=True,source_sha256=digest(Path(__file__))))
    print('MC_CHI_KEY_AMENDMENT',len(valid),'events; no membership change',flush=True)


if __name__=='__main__':main()
