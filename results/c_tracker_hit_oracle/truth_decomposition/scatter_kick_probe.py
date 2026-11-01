#!/usr/bin/env python3
import json
from pathlib import Path

import numpy as np
import uproot

H = Path(__file__).resolve().parent
A = np.load(str(H / "aligned.npz"))
Z = np.load(str(H / "decomposition.npz"))
with uproot.open(str(H / "native.root")) as f:
    N = f["probe"].arrays(["event_index", "idealrec", "jacrec"], library="np")
with uproot.open(str(H / "reverse.root")) as f:
    E = f["reverse"].arrays(["event_index", "yref"], library="np")
with uproot.open(str(H / "intervals.root")) as f:
    K = f["intervals"].arrays(["event_index", "kick", "offset"], library="np")

assert np.array_equal(N["event_index"], A["ref_event_index"])
assert np.array_equal(E["event_index"], A["ref_event_index"])
assert np.array_equal(K["event_index"], A["ref_event_index"])

role = A["ref_role"]
valid_truth = Z["validtruth"]
train = (role == b"train") & valid_truth
valid = (role == b"valid") & valid_truth
test = (role == b"test") & valid_truth

R = A["ref_R"].astype(float)
err = np.log(R / A["ref_R_true"].astype(float))
source = A["ref_source"].astype(np.int32)
y = A["tk_hitc"][:, 1:8, 1].astype(float)
z = A["tk_hitc"][:, 1:8, 2].astype(float)
ideal = N["idealrec"][:, 1:8, 1].astype(float)
jac = N["jacrec"][:, 1:8, 1, :][:, :, [1, 3, 4]].astype(float)

offset = K["offset"][:, 1:7, 1].astype(float)
kick = K["kick"][:, 1:7, 1].astype(float)
process = np.concatenate([offset * R[:, None], kick * R[:, None]], axis=1)


def metric(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    q = np.percentile(x, [2.5, 16, 50, 84, 97.5])
    return {
        "n": int(len(x)),
        "median": float(q[2]),
        "w68": float((q[3] - q[1]) / 2),
        "w95": float((q[4] - q[0]) / 2),
        "std": float(np.std(x)),
    }


def w68(x):
    return metric(x)["w68"]


def robust_cov(x, mask, p=1.0):
    lo, hi = np.percentile(x[mask], [p, 100 - p], axis=0)
    t = np.minimum(np.maximum(x[mask], lo), hi)
    c = np.cov(t, rowvar=False)
    scale = np.trace(c) / c.shape[0]
    return c + np.eye(c.shape[0]) * (1e-10 * scale + 1e-16)


def ridge_multi(x, target, mtrain, mvalid, lambdas):
    good = np.all(np.isfinite(x), axis=1) & np.all(np.isfinite(target), axis=1)
    mt = mtrain & good
    mv = mvalid & good
    mu_x = np.mean(x[mt], axis=0)
    sd_x = np.std(x[mt], axis=0)
    sd_x[sd_x < 1e-9] = 1.0
    xx = np.clip((x - mu_x) / sd_x, -10, 10)
    mu_y = np.mean(target[mt], axis=0)
    yy = target - mu_y

    best = None
    gram = xx[mt].T.dot(xx[mt]) / mt.sum()
    rhs = xx[mt].T.dot(yy[mt]) / mt.sum()
    for lam in lambdas:
        beta = np.linalg.solve(gram + lam * np.eye(gram.shape[0]), rhs)
        pred = xx.dot(beta) + mu_y
        score = np.mean([
            w68(target[mv, j] - pred[mv, j]) / w68(target[mv, j])
            for j in range(target.shape[1])
        ])
        if best is None or score < best["score"]:
            best = {
                "score": float(score),
                "lambda": float(lam),
                "beta": beta,
                "prediction": pred,
            }
    return best


def interval_metrics(pred, mask):
    rows = []
    for j in range(6):
        t = process[mask, 6 + j]
        p = pred[mask, 6 + j]
        rows.append({
            "from_layer": j + 2,
            "to_layer": j + 3,
            "corr": float(np.corrcoef(p, t)[0, 1]),
            "truth": metric(t),
            "residual": metric(t - p),
            "w68_reduction_pct": float(100 * (1 - w68(t - p) / w68(t))),
            "prediction_std_fraction": float(np.std(p) / np.std(t)),
        })
    return rows


def bootstrap_gain(before, after, src, mask, n=600, seed=20261002):
    b = np.asarray(before)[mask]
    a = np.asarray(after)[mask]
    s = np.asarray(src)[mask]
    ids = np.unique(s)
    groups = [np.flatnonzero(s == k) for k in ids]
    rng = np.random.RandomState(seed)
    vals = []
    for _ in range(n):
        ix = np.concatenate([groups[i] for i in rng.randint(len(groups), size=len(groups))])
        vals.append(100 * (1 - w68(a[ix]) / w68(b[ix])))
    return [float(x) for x in np.percentile(vals, [2.5, 50, 97.5])]


# GBL-like linear-Gaussian posterior.
measurement_cov = robust_cov(Z["measurement_hit_cm"].astype(float), train)
process_cov = robust_cov(process, train)
residual_mean = np.median((y - ideal)[train], axis=0)

gbl_process = np.zeros_like(process)
gbl_dlogr = np.zeros(len(R))
for start in range(0, len(R), 1500):
    ii = np.arange(start, min(start + 1500, len(R)))
    g = np.zeros((len(ii), 7, 12))
    for j in range(7):
        for k in range(6):
            if k < j:
                g[:, j, k] = 1.0 / R[ii]
                if j > k + 1:
                    g[:, j, 6 + k] = (z[ii, j] - z[ii, k + 1]) / R[ii]

    c = measurement_cov[None, :, :] + np.einsum(
        "nij,jk,nlk->nil", g, process_cov, g
    )
    inv_c = np.linalg.inv(c)
    jj = jac[ii]
    rr = (y - ideal)[ii] - residual_mean
    jw = np.einsum("npi,npq->niq", jj, inv_c)
    info = np.einsum("nip,npj->nij", jw, jj)
    rhs = np.einsum("niq,nq->ni", jw, rr)
    aa = np.linalg.solve(info + np.eye(3)[None, :, :] * 1e-12, rhs)
    left = rr - np.einsum("nij,nj->ni", jj, aa)
    weighted = np.einsum("nij,nj->ni", inv_c, left)
    gbl_process[ii] = np.einsum("ab,nib,ni->na", process_cov, g, weighted)
    gbl_dlogr[ii] = aa[:, 2]

# Direct triplet observable.
reverse_residual = y - E["yref"].astype(float)
segment_slope = np.diff(reverse_residual, axis=1) / np.diff(z, axis=1)
triplet_kink = np.diff(segment_slope, axis=1) * R[:, None]

mc_index = np.array([7, 8, 9, 10, 11, 12, 13])
truth_y = A["mevcoo1"][:, mc_index, 1].astype(float)
truth_dir = A["mevdir1"][:, mc_index].astype(float)
truth_y += (
    truth_dir[:, :, 1] / truth_dir[:, :, 2]
) * (z - A["mevcoo1"][:, mc_index, 2])
truth_residual = truth_y - E["yref"].astype(float)
truth_segment_slope = np.diff(truth_residual, axis=1) / np.diff(z, axis=1)
truth_triplet_kink = np.diff(truth_segment_slope, axis=1) * R[:, None]

visible_triplets = []
for j in range(5):
    t = truth_triplet_kink[test, j]
    p = triplet_kink[test, j]
    visible_triplets.append({
        "center_layer": j + 3,
        "corr": float(np.corrcoef(p, t)[0, 1]),
        "truth": metric(t),
        "measurement_residual": metric(t - p),
        "w68_reduction_pct": float(100 * (1 - w68(t - p) / w68(t))),
    })

triplet_fit = ridge_multi(
    triplet_kink,
    process,
    train,
    valid,
    [0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1],
)
triplet_process = triplet_fit["prediction"]

# A small supervised process model using the existing 95 measured-only reverse features.
feature_names = [
    q.strip()
    for q in (H / "features_reverse_noeloss.txt").read_text().splitlines()
    if q.strip()
]
with uproot.open(str(H / "reverse_features.root")) as f:
    B = f["events"].arrays(["event_index"] + feature_names, library="np")
lookup = {int(v): i for i, v in enumerate(A["ref_event_index"])}
feature_index = np.array([lookup[int(v)] for v in B["event_index"]], dtype=np.int64)
x_full = np.column_stack([B[n].astype(float) for n in feature_names])
target_full = process[feature_index]
train_full = train[feature_index]
valid_full = valid[feature_index]
test_full = test[feature_index]
full_fit = ridge_multi(
    x_full,
    target_full,
    train_full,
    valid_full,
    [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10],
)
full_process = np.full_like(process, np.nan)
full_process[feature_index] = full_fit["prediction"]

# Train one fixed physical map from true process innovations to the truth transport projection.
mu_process = np.mean(process[train], axis=0)
xp = process - mu_process
transport_target = Z["transport"].astype(float)
transport_beta = np.linalg.lstsq(
    xp[train],
    transport_target[train] - np.mean(transport_target[train]),
    rcond=None,
)[0]


def process_to_rigidity(pred, label):
    pred_transport = (pred - mu_process).dot(transport_beta)
    finite = np.isfinite(pred_transport)
    mv = valid & finite
    mt = test & finite
    best = None
    for alpha in np.arange(-1.5, 1.5001, 0.1):
        corrected = err - alpha * pred_transport
        score = w68(corrected[mv])
        if best is None or score < best[0]:
            best = (score, float(alpha))
    alpha = best[1]
    corrected = err - alpha * pred_transport
    return {
        "label": label,
        "alpha_selected_on_validation": alpha,
        "corr_pred_transport_raw_error_test": float(
            np.corrcoef(pred_transport[mt], err[mt])[0, 1]
        ),
        "corr_pred_transport_truth_transport_test": float(
            np.corrcoef(pred_transport[mt], transport_target[mt])[0, 1]
        ),
        "raw_test": metric(err[mt]),
        "corrected_test": metric(corrected[mt]),
        "gain_point_pct": float(100 * (1 - w68(corrected[mt]) / w68(err[mt]))),
        "gain_bootstrap_CI_pct": bootstrap_gain(err, corrected, source, mt),
        "prediction": pred_transport,
        "corrected": corrected,
        "finite_test_mask": mt,
    }


rig_gbl = process_to_rigidity(gbl_process, "gbl_like")
rig_triplet = process_to_rigidity(triplet_process, "triplet_linear")
rig_full = process_to_rigidity(full_process, "supervised_95")

# Direct global curvature correction from the GBL-like joint fit.
best_direct = None
for alpha in np.arange(-1.5, 1.5001, 0.1):
    corrected = err + alpha * gbl_dlogr
    score = w68(corrected[valid])
    if best_direct is None or score < best_direct[0]:
        best_direct = (score, float(alpha))
gbl_direct_alpha = best_direct[1]
gbl_direct_corrected = err + gbl_direct_alpha * gbl_dlogr

# Incremental closure against the fair existing-geometry retraining used by reverse_results.json.
def load_stage(stage):
    rows = []
    for b in range(4):
        with uproot.open(str(H / "tmva_reverse" / ("%s_bin%d.root" % (stage, b)))) as f:
            key = next(k for k in f.keys(recursive=True) if "TestTree" in k)
            rows.append(
                f[key].arrays(
                    ["event_index", "role", "source", "target", "BDTG"],
                    library="np",
                )
            )
    return {k: np.concatenate([r[k] for r in rows]) for k in rows[0]}


base = load_stage("base_fair")
order = np.argsort(base["event_index"])
base_ids = base["event_index"][order].astype(np.int64)
base_role = base["role"][order].astype(np.int8)
base_source = base["source"][order].astype(np.int32)
base_raw = base["target"][order].astype(float)
base_pred = np.asarray(base["BDTG"][order]).reshape(-1).astype(float)
base_alpha = json.loads((H / "reverse_results.json").read_text())[
    "alpha_selected_on_validation"
]["base_fair"]
base_error = base_raw - base_alpha * base_pred
a_index = np.array([lookup[int(v)] for v in base_ids], dtype=np.int64)


def add_to_base(rig):
    extra = rig["prediction"][a_index]
    mv = base_role == 1
    mt = base_role == 2
    best = None
    for alpha in np.arange(-1.5, 1.5001, 0.1):
        corrected = base_error - alpha * extra
        score = w68(corrected[mv])
        if best is None or score < best[0]:
            best = (score, float(alpha))
    alpha = best[1]
    corrected = base_error - alpha * extra
    gain = bootstrap_gain(base_error, corrected, base_source, mt)
    return {
        "alpha_selected_on_validation": alpha,
        "base_test": metric(base_error[mt]),
        "augmented_test": metric(corrected[mt]),
        "increment_point_pct": float(
            100 * (1 - w68(corrected[mt]) / w68(base_error[mt]))
        ),
        "increment_bootstrap_CI_pct": gain,
    }


result = {
    "definition": {
        "sample": "C12, reconstructed rigidity 2.5-6 GV, source-file-disjoint train/valid/test",
        "process_target": "six interval endpoint offsets times R plus six bending-plane angular kicks times R; MC truth is used only for model training/evaluation and train-derived priors",
        "gbl_like": "linear-Gaussian broken-line posterior around the reconstructed-state native-field reference; 3 global parameters plus 12 process variables",
        "triplet": "five measured local residual-slope kinks; a train-only ridge map predicts the 12 true process variables",
        "supervised_95": "ridge regression from the existing 95 measured-only reverse-trajectory features to the 12 true process variables",
        "rigidity_projection": "one train-only linear map from true process variables to the truth transport projection; a single scale is selected on validation before test evaluation",
    },
    "counts": {
        "train": int(train.sum()),
        "valid": int(valid.sum()),
        "test": int(test.sum()),
    },
    "visible_triplet_truth_position_kinks": visible_triplets,
    "gbl_like": {
        "interval_kicks": interval_metrics(gbl_process, test),
        "direct_dlogR": {
            "alpha_selected_on_validation": gbl_direct_alpha,
            "corr_dlogR_raw_error_test": float(
                np.corrcoef(gbl_dlogr[test], err[test])[0, 1]
            ),
            "corrected_test": metric(gbl_direct_corrected[test]),
            "gain_point_pct": float(
                100 * (1 - w68(gbl_direct_corrected[test]) / w68(err[test]))
            ),
            "gain_bootstrap_CI_pct": bootstrap_gain(
                err, gbl_direct_corrected, source, test
            ),
        },
        "process_projected_to_rigidity": {
            k: v for k, v in rig_gbl.items()
            if k not in ("prediction", "corrected", "finite_test_mask")
        },
    },
    "triplet_linear": {
        "lambda_selected_on_validation": triplet_fit["lambda"],
        "validation_mean_process_residual_ratio": triplet_fit["score"],
        "interval_kicks": interval_metrics(triplet_process, test),
        "process_projected_to_rigidity": {
            k: v for k, v in rig_triplet.items()
            if k not in ("prediction", "corrected", "finite_test_mask")
        },
    },
    "supervised_95": {
        "feature_count": len(feature_names),
        "lambda_selected_on_validation": full_fit["lambda"],
        "validation_mean_process_residual_ratio": full_fit["score"],
        "interval_kicks": interval_metrics(full_process, test),
        "process_projected_to_rigidity": {
            k: v for k, v in rig_full.items()
            if k not in ("prediction", "corrected", "finite_test_mask")
        },
    },
    "incremental_over_fair_existing_baseline": {
        "gbl_like": add_to_base(rig_gbl),
        "triplet_linear": add_to_base(rig_triplet),
        "supervised_95": add_to_base(rig_full),
        "note": "The fair baseline is the source-disjoint retraining in reverse_results.json (test W68 about 0.1001), used because it has both validation and test predictions for honest incremental scale selection.",
    },
}

(H / "scatter_kick_results.json").write_text(json.dumps(result, indent=2))
np.savez_compressed(
    str(H / "scatter_kick_predictions.npz"),
    event_index=A["ref_event_index"],
    gbl_process=gbl_process.astype(np.float32),
    triplet_process=triplet_process.astype(np.float32),
    supervised_process=full_process.astype(np.float32),
    true_process=process.astype(np.float32),
    gbl_dlogr=gbl_dlogr.astype(np.float32),
    triplet_visible=triplet_kink.astype(np.float32),
    truth_triplet_visible=truth_triplet_kink.astype(np.float32),
    pred_transport_gbl=rig_gbl["prediction"].astype(np.float32),
    pred_transport_triplet=rig_triplet["prediction"].astype(np.float32),
    pred_transport_supervised=rig_full["prediction"].astype(np.float32),
)
print(json.dumps(result, indent=2))
