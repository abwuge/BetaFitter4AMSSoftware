#!/usr/bin/python3
import csv
import os
import subprocess
import sys
PROD=os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(PROD,'repair_assignments.tsv')) as f:
    assignments={r['successor_input']:r for r in csv.DictReader(f,delimiter='\t')}
args=sys.argv[1:]
repair=assignments.get(args[0])
wrapper=os.path.join(PROD,'run_job_safe.sh')
if repair and not os.path.exists(repair['failed_output']):
    print('REPAIR',repair['failed_input'],flush=True)
    subprocess.check_call([wrapper,repair['failed_input'],repair['failed_output'],'-1'])
subprocess.check_call([wrapper]+args)
