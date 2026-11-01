#!/usr/bin/env python3
"""Verify all period fits and retire only their held, unstarted duplicates."""
import argparse
import json
import subprocess
from pathlib import Path
from research_periodic_adc import ROOT,OUT
from research_flight_conditional import write,digest
from research_conditional_mc import stamp
BATCH=OUT/'batch'


def main(retire=False):
    target=OUT/'execution_verification.json'
    if target.exists():raise FileExistsError('Period execution was already verified')
    protocol=json.loads((BATCH/'protocol.json').read_text());rows=[]
    for index,(dataset,family,variant,start) in enumerate(protocol['tasks']):
        receipt_path=BATCH/('task_%02d.json'%index);receipt=json.loads(receipt_path.read_text())
        path=OUT/dataset/(family+'_'+variant)/('start_%d.json'%start);result=json.loads(path.read_text())
        if receipt['exit_code'] or not receipt['numerically_converged'] or not result['numerically_converged']:
            raise ValueError('Period task not successful: %d'%index)
        if digest(path)!=receipt['result_sha256'] or receipt['protocol_sha256']!=digest(BATCH/'protocol.json'):
            raise ValueError('Period result provenance differs: %d'%index)
        if result['protocol_sha256']!=digest(OUT/dataset/'protocol.json'):raise ValueError('Scientific protocol changed')
        rows.append(dict(index=index,task=protocol['tasks'][index],host=receipt['host'],result_sha256=digest(path),receipt_sha256=digest(receipt_path)))
    completion=json.loads((OUT/'fit_completion.json').read_text())
    for name,fingerprint in completion['freeze_sha256'].items():
        if digest(ROOT/name)!=fingerprint:raise ValueError('Period freeze changed')
    old=json.loads((OUT.parent/'prospective_flight/verification.json').read_text());inherited=json.loads((OUT.parent/'inherited_status.json').read_text())
    protected=dict(old['files_sha256']);protected.update(inherited['protected'])
    for name,fingerprint in protected.items():
        if digest(ROOT/name)!=fingerprint:raise ValueError('Inherited protected artifact changed '+name)
    cluster=json.loads((BATCH/'submission.json').read_text())['cluster'];packed=json.loads((BATCH/'packed/protocol.json').read_text())
    query=['condor_q',str(cluster),'-json','-attributes','ClusterId,ProcId,JobStatus,RemoteHost']
    queue=json.loads(subprocess.check_output(query).decode() or '[]');expected=set(packed['tasks'])
    if any(row['ProcId'] not in expected or row['JobStatus']!=5 or row.get('RemoteHost') for row in queue):
        raise ValueError('Original queue contains something other than held packed-task duplicates')
    record=dict(complete=True,created_utc=stamp(),tasks=rows,models=15,starts=30,
        actual_unique_fit_hosts=sorted(set(row['host'] for row in rows)),
        distinct_host_claim='The original thirty-host dispatch was capacity-limited. Actual receipt hosts, not job count, are authoritative.',
        inherited_artifacts_unchanged=len(protected),fit_completion_sha256=digest(OUT/'fit_completion.json'),
        pending_original_queue=queue,retirement_requested=retire,source_sha256=digest(Path(__file__)))
    write(OUT/'execution_before_retirement.json',record)
    if retire and queue:
        record['retirement_reply']=subprocess.check_output(['condor_rm','-constraint','ClusterId==%d && JobStatus==5'%cluster]).decode()
        record['queue_after_retirement']=json.loads(subprocess.check_output(query).decode() or '[]')
        if record['queue_after_retirement']:raise ValueError('Held duplicate retirement incomplete')
    write(target,record)
    print('PERIOD_EXECUTION_VERIFIED',len(rows),'starts',len(record['actual_unique_fit_hosts']),'actual hosts',len(protected),'protected artifacts; held duplicates retired',bool(retire),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--retire',action='store_true');args=parser.parse_args();main(args.retire)
