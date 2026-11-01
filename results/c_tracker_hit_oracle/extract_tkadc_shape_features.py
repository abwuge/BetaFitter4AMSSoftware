#!/usr/bin/env python3
"""Build a C12-only safe strip/cluster-morphology rigidity-oracle dataset.

The new features intentionally remove absolute charge/dE/dx scale.  They keep
only normalized strip-shape information, local strip-address metadata, and
ratios among layer charge estimators.  The baseline columns are exactly the
previous strict_geometry feature list.
"""
from __future__ import print_function

import csv
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import uproot

HERE = Path(__file__).resolve().parent
MASS = HERE.parent / "c_mass_toy_oracle" / "predictions.npz"
BASE_ROOT = HERE / "tracker_features.root"
BASE_LIST = HERE / "features_strict_geometry.txt"
PROD = HERE / "tkadc_production"
JOBS = PROD / "jobs.tsv"
OUT = HERE / "tkadc_shape_features.root"
META = HERE / "tkadc_shape_meta.json"

NEW_BRANCHES = ["run", "event", "tk_adc", "tk_qlqs", "tk_stripx", "tk_rigidityk"]


def log(message):
    print(message, flush=True)


def w68(values):
    q = np.percentile(np.asarray(values, dtype=float), [16.0, 84.0])
    return float((q[1] - q[0]) / 2.0)


def add_feature(columns, names, values, prefix):
    values = np.asarray(values, dtype=np.float64)
    flat = values.reshape(len(values), -1)
    columns.append(flat)
    shape = values.shape[1:]
    if not shape:
        names.append(prefix)
    else:
        for index in np.ndindex(shape):
            names.append(prefix + "_" + "_".join(str(i) for i in index))


def derive_safe_features(adc, qlqs, stripx):
    """Use strip morphology, not total ionization amplitude."""
    adc = np.asarray(adc, dtype=np.float64)
    qlqs = np.asarray(qlqs, dtype=np.float64)
    stripx = np.asarray(stripx, dtype=np.int64)
    columns = []
    names = []
    groups = {}

    # Strip ADC shape. Negative/non-finite signals do not carry usable shape.
    signal = np.where(np.isfinite(adc) & (adc > 0.0), adc, 0.0)
    total = np.sum(signal, axis=-1)
    good = total > 0.0
    frac = np.zeros_like(signal)
    frac[good] = signal[good] / total[good, None]
    start = sum(x.shape[1] for x in columns)
    add_feature(columns, names, frac, "adcfrac")
    coord = np.arange(-5.0, 6.0)
    centroid = np.sum(frac * coord, axis=-1)
    variance = np.sum(frac * (coord[None, None, None, :] - centroid[..., None]) ** 2, axis=-1)
    asymmetry = np.sum(frac[..., 6:], axis=-1) - np.sum(frac[..., :5], axis=-1)
    seed = frac[..., 5]
    nearest = np.sum(frac[..., 4:7], axis=-1)
    tails = np.sum(frac[..., [0, 1, 9, 10]], axis=-1)
    active = np.sum(signal > 0.0, axis=-1).astype(np.float64)
    present = good.astype(np.float64)
    for values, prefix in [
        (centroid, "adc_centroid"),
        (variance, "adc_variance"),
        (asymmetry, "adc_asymmetry"),
        (seed, "adc_seed_fraction"),
        (nearest, "adc_nearest_fraction"),
        (tails, "adc_tail_fraction"),
        (active, "adc_active_strips"),
        (present, "adc_present"),
    ]:
        add_feature(columns, names, values, prefix)
    groups["adcshape"] = [start, sum(x.shape[1] for x in columns)]

    # qlqs: remove the event/layer charge magnitude by centering log Q across
    # the five charge estimators for each layer and side.
    qgood = np.isfinite(qlqs) & (qlqs > 0.0)
    logq = np.full(qlqs.shape, np.nan, dtype=np.float64)
    logq[qgood] = np.log(qlqs[qgood])
    count = np.sum(qgood, axis=-1)
    summ = np.nansum(logq, axis=-1)
    mean = np.divide(summ, count, out=np.zeros_like(summ), where=count > 0)
    centered = logq - mean[..., None]
    centered[~qgood] = np.nan
    start = sum(x.shape[1] for x in columns)
    add_feature(columns, names, centered, "qlqs_centered_log")
    add_feature(columns, names, qgood.astype(np.float64), "qlqs_present")
    groups["qlqs_ratio"] = [start, sum(x.shape[1] for x in columns)]

    # Decode tk_stripx = ik + sensor*10 + sensor_cad*1000 + noverflow*1e6.
    # noverflow may be -1; floor division plus remainder recovers the address.
    x = stripx.copy()
    noverflow = np.floor_divide(x, 1000000)
    rest = x - noverflow * 1000000
    sensor_cad = np.floor_divide(rest, 1000)
    rest2 = rest - sensor_cad * 1000
    sensor = np.floor_divide(rest2, 10)
    ik = rest2 - sensor * 10
    plausible = (sensor_cad >= 0) & (sensor_cad < 300) & (sensor >= 0) & (sensor < 20) & (ik >= 0) & (ik < 10)
    start = sum(x.shape[1] for x in columns)
    for values, prefix in [
        (np.where(plausible, ik, np.nan), "strip_type"),
        (np.where(plausible, sensor, np.nan), "strip_sensor"),
        (np.where(plausible, sensor_cad, np.nan), "strip_local_address"),
        (np.where(plausible, noverflow, np.nan), "strip_overflow_code"),
        (plausible.astype(np.float64), "strip_decode_valid"),
    ]:
        add_feature(columns, names, values, prefix)
    groups["strip_address"] = [start, sum(x.shape[1] for x in columns)]
    return np.column_stack(columns).astype(np.float32), names, groups


def one_file(task):
    local_rows, path, target_run, target_event, target_R = task
    started = time.time()
    with uproot.open(path) as source:
        tree = source["amstreea"]
        missing = [name for name in NEW_BRANCHES if name not in tree.keys()]
        if missing:
            raise RuntimeError("%s missing %s" % (path, missing))
        arrays = tree.arrays(NEW_BRANCHES, library="np")
    by_id = {key: i for i, key in enumerate(zip(arrays["run"].tolist(), arrays["event"].tolist()))}
    if len(by_id) != len(arrays["run"]):
        raise RuntimeError("duplicate run/event identities in " + path)
    targets = list(zip(target_run.tolist(), target_event.tolist()))
    positions = np.asarray([i for i, key in enumerate(targets) if key in by_id], dtype=np.int64)
    entries = np.asarray([by_id[targets[i]] for i in positions], dtype=np.int64)
    rows = local_rows[positions]
    adc = arrays["tk_adc"][entries]
    qlqs = arrays["tk_qlqs"][entries]
    safe, names, groups = derive_safe_features(adc, qlqs, arrays["tk_stripx"][entries])
    # Omitted absolute charge scales are retained only as diagnostics/spectators.
    adc_signal = np.where(np.isfinite(adc) & (adc > 0.0), adc, 0.0)
    adc_total = np.sum(adc_signal, axis=-1)
    adc_logs = np.where(adc_total > 0.0, np.log(np.maximum(adc_total, 1e-12)), np.nan)
    adc_log_scale = np.nanmedian(adc_logs.reshape(len(adc_logs), -1), axis=1)
    q_good = np.isfinite(qlqs) & (qlqs > 0.0)
    q_logs = np.where(q_good, np.log(np.maximum(qlqs, 1e-12)), np.nan)
    qlqs_log_scale = np.nanmedian(q_logs.reshape(len(q_logs), -1), axis=1)
    spectators = {
        "adc_log_scale": adc_log_scale.astype(np.float32),
        "qlqs_log_scale": qlqs_log_scale.astype(np.float32),
    }
    rnew = arrays["tk_rigidityk"][entries, 1, 1].astype(np.float64)
    rold = target_R[positions].astype(np.float64)
    relative = np.abs(rnew - rold) / np.maximum(np.abs(rold), 1e-12)
    audit = {
        "path": path,
        "target": int(len(local_rows)),
        "matched": int(len(rows)),
        "coverage": float(len(rows) / float(max(1, len(local_rows)))),
        "median_relative_R_replay": float(np.median(relative)) if len(relative) else float("nan"),
        "max_relative_R_replay": float(np.max(relative)) if len(relative) else float("nan"),
        "seconds": time.time() - started,
    }
    return rows, safe, names, groups, spectators, audit


def fit_preprocessor(raw, names, train):
    values = np.asarray(raw[train], dtype=np.float64)
    finite_fraction = np.mean(np.isfinite(values), axis=0)
    keep = finite_fraction >= 0.01
    indices = np.flatnonzero(keep)
    values = values[:, keep]
    kept_names = [name for name, selected in zip(names, keep) if selected]
    finite_fraction = finite_fraction[keep]
    low = np.nanpercentile(values, 0.5, axis=0)
    high = np.nanpercentile(values, 99.5, axis=0)
    clipped = np.minimum(np.maximum(values, low), high)
    center = np.nanmedian(clipped, axis=0)
    q25, q75 = np.nanpercentile(clipped, [25.0, 75.0], axis=0)
    robust = 0.741301109 * (q75 - q25)
    standard = np.nanstd(clipped, axis=0)
    scale = np.maximum(robust, 0.10 * standard)
    informative = np.isfinite(center) & np.isfinite(scale) & (scale > 1e-7)
    indices = indices[informative]
    kept_names = [name for name, selected in zip(kept_names, informative) if selected]
    low = low[informative]
    high = high[informative]
    center = center[informative]
    scale = scale[informative]
    finite_fraction = finite_fraction[informative]
    missing_indicator = (1.0 - finite_fraction) >= 0.01
    state = {
        "indices": indices.tolist(),
        "low": low.tolist(), "high": high.tolist(),
        "center": center.tolist(), "scale": scale.tolist(),
        "missing_indicator": missing_indicator.tolist(),
        "feature_names": kept_names,
    }
    return state


def apply_preprocessor(raw, state):
    indices = np.asarray(state["indices"], dtype=np.int32)
    low = np.asarray(state["low"])
    high = np.asarray(state["high"])
    center = np.asarray(state["center"])
    scale = np.asarray(state["scale"])
    miss_ind = np.asarray(state["missing_indicator"], dtype=bool)
    values = np.asarray(raw[:, indices], dtype=np.float64)
    missing = ~np.isfinite(values)
    values = np.minimum(np.maximum(values, low), high)
    values = np.where(missing, center, values)
    z = np.clip((values - center) / scale, -6.0, 6.0)
    columns = [z.astype(np.float32)]
    names = ["s%04d" % i for i in range(z.shape[1])]
    if np.any(miss_ind):
        columns.append(missing[:, miss_ind].astype(np.float32))
        names += ["s%04d" % i for i in range(z.shape[1], z.shape[1] + int(np.sum(miss_ind)))]
    return np.column_stack(columns), names


def main():
    started = time.time()
    with np.load(str(MASS), allow_pickle=False) as z:
        mask = z["label"] == 0
        global_index = np.flatnonzero(mask)
        data = {key: z[key][mask] for key in ["run","event","source","label","bin","R","R_true","u","u_true","role"]}
    roles = np.char.decode(data["role"], "ascii")
    train = roles == "train"
    valid = roles == "valid"
    test = roles == "test"
    base_names = [line.strip() for line in BASE_LIST.read_text().splitlines() if line.strip()]
    with uproot.open(str(BASE_ROOT)) as source:
        tree = source["events"]
        base_arrays = tree.arrays(base_names + ["event_index"], library="np")
    if not np.array_equal(base_arrays["event_index"][global_index], global_index.astype(base_arrays["event_index"].dtype)):
        raise RuntimeError("base feature event_index alignment failure")
    base = np.column_stack([base_arrays[name][global_index] for name in base_names]).astype(np.float32)

    jobs = list(csv.DictReader(open(str(JOBS)), delimiter="\t"))
    by_source = {int(row["source"]): row for row in jobs}
    raw_new = None
    new_names = None
    groups = None
    matched = np.zeros(len(global_index), dtype=bool)
    audits = []
    diagnostic_spectators = {
        "adc_log_scale": np.full(len(global_index), np.nan, dtype=np.float32),
        "qlqs_log_scale": np.full(len(global_index), np.nan, dtype=np.float32),
    }
    tasks = []
    for source_id in sorted(int(x) for x in np.unique(data["source"])):
        rows = np.flatnonzero(data["source"] == source_id)
        row = by_source[source_id]
        path = row["output_root"]
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            raise RuntimeError("missing production output " + path)
        tasks.append((rows, path, data["run"][rows], data["event"][rows], data["R"][rows]))
    with mp.Pool(processes=int(os.environ.get("TKADC_WORKERS", "6"))) as pool:
        for done, (rows, local, names, local_groups, spectators, audit) in enumerate(pool.imap_unordered(one_file, tasks, chunksize=1), 1):
            if raw_new is None:
                raw_new = np.full((len(global_index), local.shape[1]), np.nan, dtype=np.float32)
                new_names = names
                groups = local_groups
            if names != new_names:
                raise RuntimeError("new feature schema changed")
            raw_new[rows] = local
            for key in diagnostic_spectators:
                diagnostic_spectators[key][rows] = spectators[key]
            matched[rows] = True
            audits.append(audit)
            if done % 24 == 0 or done == len(tasks):
                log("extract %d/%d files matched=%d/%d elapsed=%.1fs" % (
                    done, len(tasks), int(np.sum(matched)), len(matched), time.time()-started))
    if raw_new is None:
        raise RuntimeError("no new features")
    matched_train = matched & train
    state = fit_preprocessor(raw_new, new_names, matched_train)
    new, s_names = apply_preprocessor(raw_new, state)

    selected = matched
    # Keep only events reproduced by the new ntuple; report the selection loss.
    global_sel = global_index[selected]
    role_code = np.where(roles[selected] == "train", 0, np.where(roles[selected] == "valid", 1, 2)).astype(np.int8)
    payload = {
        "event_index": global_sel.astype(np.int32),
        "role": role_code,
        "label": np.zeros(np.sum(selected), dtype=np.int8),
        "source": data["source"][selected].astype(np.int32),
        "bin": data["bin"][selected].astype(np.int8),
        "R": data["R"][selected].astype(np.float32),
        "R_true": data["R_true"][selected].astype(np.float32),
        "u": data["u"][selected].astype(np.float32),
        "u_true": data["u_true"][selected].astype(np.float32),
        "target": np.log(data["R"][selected] / data["R_true"][selected]).astype(np.float32),
        "adc_log_scale": diagnostic_spectators["adc_log_scale"][selected],
        "qlqs_log_scale": diagnostic_spectators["qlqs_log_scale"][selected],
    }
    base_sel = base[selected]
    new_sel = new[selected]
    for i, name in enumerate(base_names):
        payload[name] = base_sel[:, i]
    for i, name in enumerate(s_names):
        payload[name] = new_sel[:, i]
    with uproot.recreate(str(OUT), compression=uproot.ZLIB(4)) as output:
        output["events"] = payload

    # Feature-list stages. All new inputs are safe morphology features with
    # charge magnitude removed.
    adc_end_raw = groups["adcshape"][1]
    q_end_raw = groups["qlqs_ratio"][1]
    # Map raw-group boundaries onto preprocessed retained columns by original index.
    retained = np.asarray(state["indices"], dtype=np.int32)
    base_s = ["s%04d" % i for i in range(len(state["feature_names"]))]
    miss = np.asarray(state["missing_indicator"], dtype=bool)
    miss_s = ["s%04d" % i for i in range(len(base_s), len(s_names))]
    def select_new(max_raw_index):
        keep = [base_s[i] for i, original in enumerate(retained) if original < max_raw_index]
        miss_positions = np.flatnonzero(miss)
        keep += [miss_s[j] for j, pos in enumerate(miss_positions) if retained[pos] < max_raw_index]
        return keep
    adc_features = select_new(adc_end_raw)
    adc_set = set(adc_features)
    qlqs_features = [name for name in select_new(q_end_raw) if name not in adc_set]
    qlqs_set = set(qlqs_features)
    strip_features = [name for name in s_names if name not in adc_set and name not in qlqs_set]
    all_features = list(s_names)
    lists = {
        "features_tkadc_base.txt": base_names,
        "features_tkadc_adcshape.txt": base_names + adc_features,
        "features_tkadc_adcstrip.txt": base_names + adc_features + strip_features,
        "features_tkadc_morph.txt": base_names + all_features,
        "features_tkadc_shape_only.txt": all_features,
    }
    for filename, features in lists.items():
        (HERE / filename).write_text("\n".join(features) + "\n")

    replay = [x["max_relative_R_replay"] for x in audits if np.isfinite(x["max_relative_R_replay"])]
    summary = {
        "events_c12": int(len(global_index)),
        "matched": int(np.sum(matched)),
        "coverage": float(np.mean(matched)),
        "role_counts_before": {r:int(np.sum(roles==r)) for r in ["train","valid","test"]},
        "role_counts_matched": {r:int(np.sum(matched & (roles==r))) for r in ["train","valid","test"]},
        "raw_new_features": int(raw_new.shape[1]),
        "retained_new_features": int(new.shape[1]),
        "base_features": int(base.shape[1]),
        "raw_group_boundaries": groups,
        "max_relative_R_replay": float(max(replay)) if replay else float("nan"),
        "median_file_coverage": float(np.median([x["coverage"] for x in audits])),
        "min_file_coverage": float(np.min([x["coverage"] for x in audits])),
        "preprocessor": state,
        "audits": audits,
        "definition": "New features exclude absolute ADC/charge magnitude. tk_adc is normalized within each layer/side; tk_qlqs enters only as centered log ratios; tk_stripx contributes decoded local address/type/overflow metadata.",
        "elapsed_seconds": time.time() - started,
    }
    META.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    log("DONE matched %d/%d (%.4f), base=%d new=%d max_R_replay=%.3g elapsed %.1fs" % (
        np.sum(matched), len(matched), np.mean(matched), base.shape[1], new.shape[1],
        summary["max_relative_R_replay"], time.time()-started))

if __name__ == "__main__":
    main()
