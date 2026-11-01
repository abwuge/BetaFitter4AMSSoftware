#!/usr/bin/env bash
set -euo pipefail
export ROOTSYS=/nas01/home/hxwu/.local
compiler=/cvmfs/sft.cern.ch/lcg/contrib/gcc/9.2.0/x86_64-centos7
export PATH="$ROOTSYS/bin:$compiler/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:$compiler/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /nas01/home/hxwu/AMSSoft/betaEstimator/results/c_tracker_hit_oracle
stage=$1
rbin=$2
mkdir -p tmva_tkadc
"$ROOTSYS/bin/root" -l -b -q "train_tracker_bdt_list.C(\"tkadc_shape_features.root\",\"tmva_tkadc\",\"${stage}\",\"features_${stage}.txt\",${rbin},300)"
/usr/bin/python3 - "$stage" "$rbin" <<'PY'
import sys
import numpy as np
import uproot
stage, rbin = sys.argv[1], int(sys.argv[2])
with uproot.open('tkadc_shape_features.root') as f:
    data = f['events'].arrays(['event_index', 'role', 'bin'], library='np')
with uproot.open('tmva_tkadc/%s_bin%d.root' % (stage, rbin)) as f:
    for token, mask in [('TrainTree', data['role'] < 2), ('TestTree', data['role'] == 2)]:
        key = next(k for k in f.keys(recursive=True) if token in k)
        a = f[key].arrays(['event_index', 'BDTG'], library='np')
        expected = data['event_index'][mask & (data['bin'] == rbin)]
        if len(a['event_index']) != len(expected) or not np.array_equal(np.sort(a['event_index']), np.sort(expected)):
            raise RuntimeError('Prediction event coverage mismatch')
        if not np.all(np.isfinite(a['BDTG'])):
            raise RuntimeError('Nonfinite prediction')
print('VALIDATED_MODEL', stage, rbin)
PY
