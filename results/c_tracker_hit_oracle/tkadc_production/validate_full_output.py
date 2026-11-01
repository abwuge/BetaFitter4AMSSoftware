import re
import sys
import numpy as np
import uproot

root_path, log_path, ids_path, limit, input_path = sys.argv[1:]
with open(log_path, errors="replace") as source:
    log = source.read()
expected = re.findall(r"^num2=(\d+)$", log, re.M)
processed = re.findall(r"^   Event\((\d+)\)", log, re.M)
if len(expected) != 1 or len(processed) != 1 or int(processed[0]) < int(expected[0]):
    raise RuntimeError("Input processing count mismatch")
if "Process Run=1 Stat=3" not in log:
    raise RuntimeError("Missing successful reconstruction completion")
with uproot.open(input_path) as source:
    input_entries = source["AMSRoot"].num_entries
    expected_entries = input_entries if int(limit) == -1 else min(input_entries, int(limit))
    if int(expected[0]) != expected_entries:
        raise RuntimeError("Input entry limit mismatch")
with uproot.open(root_path) as source:
    tree = source["amstreea"]
    if tree.num_entries <= 0 or len(tree.keys()) != 361:
        raise RuntimeError("Unexpected output schema or empty tree")
    a = tree.arrays(["run", "event", "tk_adc", "tk_qlqs", "tk_stripx"], library="np")
    for name in ["tk_adc", "tk_qlqs"]:
        if not np.count_nonzero(a[name]):
            raise RuntimeError("Empty new branch " + name)
    if int(limit) == -1:
        wanted = np.loadtxt(ids_path, dtype=np.int64, ndmin=2)
        present = set(zip(a["run"].tolist(), a["event"].tolist()))
        missing = [tuple(pair) for pair in wanted.tolist() if tuple(pair) not in present]
        if missing:
            raise RuntimeError("Missing %d target events" % len(missing))
    print("VALIDATED full output: selection_calls=%s saved=%d branches=%d numpy=%s" %
          (processed[0], tree.num_entries, len(tree.keys()), np.__version__))
