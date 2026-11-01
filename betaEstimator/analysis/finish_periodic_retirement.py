#!/usr/bin/env python3
"""Confirm asynchronous scheduler retirement after the original verification."""
import json
import subprocess
from pathlib import Path
from research_periodic_adc import ROOT,OUT
from research_flight_conditional import write,digest
from research_conditional_mc import stamp


def main():
    path=OUT/'execution_verification.json'
    if path.exists():raise FileExistsError('Execution verification already completed')
    before=OUT/'execution_before_retirement.json';record=json.loads(before.read_text())
    if not record['retirement_requested'] or len(record['tasks'])!=30:raise ValueError('Unexpected retirement record')
    for row in record['tasks']:
        dataset,family,variant,start=row['task'];result=OUT/dataset/(family+'_'+variant)/('start_%d.json'%start)
        if digest(result)!=row['result_sha256']:raise ValueError('Result changed after retirement')
        if digest(OUT/'batch'/('task_%02d.json'%row['index']))!=row['receipt_sha256']:raise ValueError('Receipt changed')
    clusters=[json.loads((OUT/'batch'/folder/'submission.json').read_text())['cluster'] for folder in ['', 'packed']]
    queue=json.loads(subprocess.check_output(['condor_q','-constraint','ClusterId==%d || ClusterId==%d'%tuple(clusters),'-json',
        '-attributes','ClusterId,ProcId,JobStatus']).decode() or '[]')
    if queue:raise ValueError('The original or packed period jobs have not all retired')
    record.update(created_utc=stamp(),queue_after_retirement=queue,original_and_packed_clusters=clusters,
        asynchronous_retirement_note='The initial immediate post-condor_rm query still contained removal-in-progress records and the checker raised an error. A subsequent empty query is now verified; no fit was rerun and no remaining queue record was assumed complete.',
        before_record_sha256=digest(before),completion_source_sha256=digest(Path(__file__)))
    write(path,record)
    print('PERIOD_EXECUTION_FINAL',30,'starts',len(record['actual_unique_fit_hosts']),'actual hosts; original and packed clusters empty',flush=True)


if __name__=='__main__':main()
