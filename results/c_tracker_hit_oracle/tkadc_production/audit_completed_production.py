#!/usr/bin/env python3
"""Audit every published file and its original analysis event identities."""
import csv
import hashlib
import json
import multiprocessing as mp
import os
import re
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
CLUSTER = 392996


def check_file(task):
    proc, row, exit_code = task
    result = {"proc": proc, "source": int(row["source"]), "role": row["role"],
              "output_root": row["output_root"], "ok": False}
    try:
        path = Path(row["output_root"])
        log = (HERE / "logs_full" / ("{role}_{source}_{run}.out".format(**row))).read_text(errors="replace")
        if exit_code != 0:
            raise RuntimeError("Missing successful scheduler termination")
        if "DONE " + str(path) not in log or "Process Run=1 Stat=3" not in log:
            raise RuntimeError("Missing producer completion")
        counts = re.findall(r"^num2=(\d+)$", log, re.M)
        calls = re.findall(r"^   Event\((\d+)\)", log, re.M)
        saved = re.findall(r"VALIDATED full output: selection_calls=\d+ saved=(\d+) branches=361", log)
        if len(counts) != 1 or len(calls) != 1 or len(saved) != 1:
            raise RuntimeError("Ambiguous producer counts")
        with uproot.open(row["input_root"]) as f:
            input_count = f["AMSRoot"].num_entries
        if int(counts[0]) != input_count or int(calls[0]) < input_count:
            raise RuntimeError("Incomplete input processing")
        with uproot.open(str(path)) as f:
            tree = f["amstreea"]
            schema = [(name, tree[name].typename) for name in tree.keys()]
            entries = tree.num_entries
            if entries != int(saved[0]) or entries <= 0 or len(schema) != 361:
                raise RuntimeError("Output count/schema mismatch")
            a = tree.arrays(["run", "event", "tk_adc", "tk_qlqs", "tk_stripx", "tk_rigidityk"], library="np")
        ids = list(zip(a["run"].tolist(), a["event"].tolist()))
        byid = {key: i for i, key in enumerate(ids)}
        if len(byid) != len(ids):
            raise RuntimeError("Duplicate run/event identities")
        wanted = np.loadtxt(str(HERE / "event_lists" / (row["run"] + ".txt")), dtype=np.int64, ndmin=2)
        if any(tuple(pair) not in byid for pair in wanted.tolist()):
            raise RuntimeError("Missing original target identities")
        with np.load(str(HERE.parent.parent / "c_mass_toy_oracle" / "predictions.npz")) as z:
            mask = (z["label"] == 0) & (z["source"] == int(row["source"]))
            targets = list(zip(z["run"][mask].tolist(), z["event"][mask].tolist()))
            if set(targets) != set(map(tuple, wanted.tolist())):
                raise RuntimeError("Event list differs from frozen analysis")
            ix = np.asarray([byid[key] for key in targets])
            relative = np.abs(a["tk_rigidityk"][ix, 1, 1] - z["R"][mask]) / np.maximum(np.abs(z["R"][mask]), 1e-12)
        if np.max(relative) != 0:
            raise RuntimeError("Rigidity replay changed")
        for name, shape in [("tk_adc", (9, 2, 11)), ("tk_qlqs", (9, 2, 5)), ("tk_stripx", (9,))]:
            if a[name].shape[1:] != shape:
                raise RuntimeError("Unexpected shape for " + name + ": " + str(a[name].shape))
            if not np.all(np.isfinite(a[name])) or not np.count_nonzero(a[name]):
                raise RuntimeError("Nonfinite or empty branch " + name)
            result[name + "_nonzero"] = int(np.count_nonzero(a[name]))
        result.update(ok=True, input_events=input_count, saved_events=entries,
                      target_events=len(targets), max_relative_R_replay=float(np.max(relative)),
                      schema_sha256=hashlib.sha256(json.dumps(schema).encode()).hexdigest(),
                      bytes=path.stat().st_size)
    except Exception as exc:
        result["error"] = str(exc)
    return result


def main():
    rows = list(csv.DictReader(open(str(HERE / "jobs.tsv")), delimiter="\t"))
    exits = {}
    for block in (HERE / "logs_full" / "condor.log").read_text(errors="replace").split("..."):
        match = re.search(r"005 \(%d\.(\d+)\.\d+\)" % CLUSTER, block)
        code = re.search(r"Normal termination \(return value (\d+)\)", block)
        if match:
            exits[int(match[1])] = int(code[1]) if code else "signal"
    tasks = [(i, row, exits.get(i)) for i, row in enumerate(rows)]
    results = []
    with mp.Pool(int(os.environ.get("TKADC_WORKERS", "6"))) as pool:
        for done, result in enumerate(pool.imap_unordered(check_file, tasks), 1):
            results.append(result)
            if done % 32 == 0:
                print("AUDIT %d/%d errors=%d" % (done, len(rows), sum(not r["ok"] for r in results)), flush=True)
    results.sort(key=lambda r: r["proc"])
    good = [r for r in results if r["ok"]]
    partial = list(Path(rows[0]["output_root"]).parent.glob("*.part.*"))
    hashes = {}
    for name in ["analysis_amsd69n_ROOT5SLC6P8HPBNV11TKADCSAFE", "analysis_amsd69n_tkadc_safe.C", "run_job_full.sh", "validate_full_output.py"]:
        hashes[name] = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
    report = {"cluster": CLUSTER, "planned_files": len(rows), "checked_files": len(good),
              "input_events": sum(r["input_events"] for r in good),
              "saved_events": sum(r["saved_events"] for r in good),
              "target_events": sum(r["target_events"] for r in good),
              "output_GiB": sum(r["bytes"] for r in good) / 2**30,
              "schema_count": len(set(r["schema_sha256"] for r in good)),
              "max_relative_R_replay": max((r["max_relative_R_replay"] for r in good), default=None),
              "partial_files": [str(p) for p in partial], "producer_sha256": hashes,
              "errors": [r for r in results if not r["ok"]], "files": results}
    report["ok"] = len(good) == len(rows) and report["schema_count"] == 1 and not partial
    (HERE / "production_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "files"}), flush=True)
    if not report["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
