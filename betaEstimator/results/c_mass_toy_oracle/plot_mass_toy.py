#!/usr/bin/env python3
"""Publication-style diagnostic plots for the supervised mass toy."""
from __future__ import print_function

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import mass_toy_common as mt

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
OUT.mkdir(exist_ok=True)

COLORS = {
    "blue": "#5790fc",
    "yellow": "#f89c20",
    "red": "#e42536",
    "grape": "#964a8b",
    "gray": "#9c9ca1",
    "violet": "#7a21dd",
    "black": "#222222",
}
METHOD_COLORS = {
    "raw_linear": COLORS["gray"],
    "raw_hybrid": COLORS["yellow"],
    "mass_proxy": COLORS["blue"],
    "mass_energy": COLORS["grape"],
    "all_physics": COLORS["red"],
    "all_detector": COLORS["violet"],
    "nuisance_control": COLORS["black"],
    "shuffle_control": "#c9c9c9",
}
DISPLAY = {
    "raw_linear": r"raw $m^2(\beta_{\rm lin},R)$",
    "raw_hybrid": r"raw $m^2(\beta_{\rm hybrid},R)$",
    "mass_proxy": "all mass proxies",
    "mass_energy": "+ energy loss",
    "all_physics": "+ RICH (physics-only)",
    "all_detector": "+ geometry/quality",
    "nuisance_control": "geometry/quality only",
    "shuffle_control": "shuffled labels",
}


def setup():
    plt.rcParams.update({
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "legend.fontsize": 9.5,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "figure.figsize": (7.4, 4.9),
        "figure.dpi": 150,
        "savefig.dpi": 220,
        "axes.linewidth": 0.9,
        "lines.linewidth": 2.0,
    })


def clean(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", length=4, width=0.9)
    ax.grid(axis="y", alpha=0.18, linewidth=0.7)


def save(fig, name):
    fig.tight_layout()
    fig.savefig(str(OUT / (name + ".png")), bbox_inches="tight")
    fig.savefig(str(OUT / (name + ".pdf")), bbox_inches="tight")
    plt.close(fig)


def read_tsv(path):
    with open(str(path), newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    for row in rows:
        for key, value in list(row.items()):
            if key in ("method", "feature", "group"):
                continue
            try:
                row[key] = float(value)
            except (TypeError, ValueError):
                pass
    return rows


def density_step(ax, values, bins, label, color, linestyle="-", linewidth=2.0):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    counts, edges = np.histogram(values, bins=bins)
    widths = np.diff(edges)
    total = float(np.sum(counts))
    density = counts / (total * widths) if total else counts.astype(float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    ax.step(centers, density, where="mid", label=label, color=color,
            linestyle=linestyle, linewidth=linewidth)


def roc_curve(score, label):
    score = np.asarray(score, dtype=np.float64)
    label = np.asarray(label, dtype=np.int8)
    order = np.argsort(-score, kind="mergesort")
    y = label[order]
    tp = np.cumsum(y == 1).astype(float)
    fp = np.cumsum(y == 0).astype(float)
    tp /= max(1.0, tp[-1])
    fp /= max(1.0, fp[-1])
    return np.r_[0.0, fp], np.r_[0.0, tp]


def main():
    setup()
    with np.load(str(HERE / "predictions.npz"), allow_pickle=False) as source:
        data = {key: source[key] for key in source.files}
    role = np.char.decode(data["role"], "ascii")
    test = role == "test"
    label = data["label"].astype(np.int8)
    event_metrics = read_tsv(HERE / "event_metrics.tsv")
    physical_metrics = read_tsv(HERE / "physical_mass_metrics.tsv")
    calibration = read_tsv(HERE / "calibration.tsv")
    with open(str(HERE / "bootstrap.json")) as source:
        bootstrap = json.load(source)

    overall = {
        row["method"]: row for row in event_metrics if int(row["bin"]) == -1
    }

    # 1. Raw physical mass peaks.
    fig, ax = plt.subplots()
    mass2 = data["mass2_raw_linear"][test]
    y = label[test]
    mass = np.full(len(mass2), np.nan)
    mass[mass2 > 0.0] = np.sqrt(mass2[mass2 > 0.0])
    bins = np.linspace(7.0, 16.5, 110)
    density_step(ax, mass[y == 0], bins, r"true $^{12}$C", COLORS["blue"])
    density_step(ax, mass[y == 1], bins, r"true $^{13}$C", COLORS["red"])
    ax.axvline(mt.MASSES[0], color=COLORS["blue"], linestyle=":", linewidth=1.5)
    ax.axvline(mt.MASSES[1], color=COLORS["red"], linestyle=":", linewidth=1.5)
    ax.set_xlabel(r"physical reconstructed mass $\sqrt{(ZR)^2(1/\beta^2-1)}$ [GeV/$c^2$]")
    ax.set_ylabel("normalized density")
    ax.set_title("Raw physical mass peaks on held-out source files")
    ax.text(0.98, 0.95, "AUC = %.3f" % overall["raw_linear"]["auc"],
            transform=ax.transAxes, ha="right", va="top")
    ax.legend(frameon=False)
    clean(ax)
    save(fig, "01_raw_physical_mass_peaks")

    # 2. Continuous posterior-mean mass, no hard snapping.
    fig, ax = plt.subplots()
    soft = data["mass_soft_all_physics"][test]
    bins = np.linspace(mt.MASSES[0] - 0.015, mt.MASSES[1] + 0.015, 125)
    density_step(ax, soft[y == 0], bins, r"true $^{12}$C", COLORS["blue"])
    density_step(ax, soft[y == 1], bins, r"true $^{13}$C", COLORS["red"])
    ax.axvline(mt.MASSES[0], color=COLORS["blue"], linestyle=":", linewidth=1.5)
    ax.axvline(mt.MASSES[1], color=COLORS["red"], linestyle=":", linewidth=1.5)
    ax.set_xlabel(r"posterior-mean mass $m_{12}+P(^{13}{\rm C}|x)(m_{13}-m_{12})$ [GeV/$c^2$]")
    ax.set_ylabel("normalized density")
    ax.set_title("Toy-model mass output: physics-only feature set")
    ax.text(0.02, 0.95, "AUC %.3f   balanced accuracy %.3f" % (
        overall["all_physics"]["auc"], overall["all_physics"]["balanced_accuracy"]),
        transform=ax.transAxes, ha="left", va="top")
    ax.legend(frameon=False)
    clean(ax)
    save(fig, "02_posterior_mean_mass_peaks")

    # 3. Posterior probability reveals ambiguity directly.
    fig, ax = plt.subplots()
    probability = data["p_c13_all_physics"][test]
    bins = np.linspace(0.0, 1.0, 81)
    density_step(ax, probability[y == 0], bins, r"true $^{12}$C", COLORS["blue"])
    density_step(ax, probability[y == 1], bins, r"true $^{13}$C", COLORS["red"])
    ax.axvline(0.5, color=COLORS["gray"], linestyle="--", linewidth=1.3)
    ax.set_xlabel(r"equal-prior posterior $P(^{13}{\rm C}|x)$")
    ax.set_ylabel("normalized density")
    ax.set_title("Event-level posterior, not just the aggregate peak")
    ax.legend(frameon=False)
    clean(ax)
    save(fig, "03_event_posterior")

    # 4. ROC hierarchy.
    fig, ax = plt.subplots()
    for method in ["raw_linear", "mass_proxy", "mass_energy", "all_physics", "all_detector", "nuisance_control"]:
        score = data["p_c13_" + method][test]
        fpr, tpr = roc_curve(score, y)
        ax.plot(fpr, tpr, color=METHOD_COLORS[method], label="%s (%.3f)" % (DISPLAY[method], overall[method]["auc"]))
    ax.plot([0, 1], [0, 1], color="#777777", linestyle=":", linewidth=1.2)
    ax.set_xlabel(r"$^{12}$C acceptance")
    ax.set_ylabel(r"$^{13}$C efficiency")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Held-out-source ROC hierarchy")
    ax.legend(frameon=False, loc="lower right")
    clean(ax)
    save(fig, "04_roc_hierarchy")

    # 5. AUC versus rigidity.
    fig, ax = plt.subplots()
    centers = 0.5 * (mt.R_EDGES[:-1] + mt.R_EDGES[1:])
    for method in ["raw_linear", "mass_proxy", "mass_energy", "all_physics", "all_detector", "nuisance_control"]:
        values = []
        for bin_id in range(4):
            row = next(item for item in event_metrics if item["method"] == method and int(item["bin"]) == bin_id)
            values.append(row["auc"])
        ax.plot(centers, values, marker="o", color=METHOD_COLORS[method], label=DISPLAY[method])
    ax.set_xlabel("reconstructed rigidity [GV]")
    ax.set_ylabel("AUC")
    ax.set_xticks(centers)
    ax.set_xticklabels(["2.5–3", "3–4", "4–5", "5–6"])
    ax.set_ylim(0.45, 1.01)
    ax.set_title("Microscopic isotope discrimination by rigidity bin")
    ax.legend(frameon=False, ncol=2)
    clean(ax)
    save(fig, "05_auc_vs_rigidity")

    # 6. Physical error decomposition: beta versus rigidity bottleneck.
    fig, ax = plt.subplots()
    physical_styles = [
        ("linear_raw", "raw β, raw R", COLORS["gray"], "o"),
        ("perfect_beta_raw_R", "truth β, raw R", COLORS["blue"], "s"),
        ("raw_beta_perfect_R", "raw β, truth R", COLORS["red"], "^"),
    ]
    for method, text, color, marker in physical_styles:
        values = []
        for bin_id in range(4):
            row = next(item for item in physical_metrics if item["method"] == method and int(item["bin"]) == bin_id)
            values.append(row["mass2_error_w68"])
        ax.plot(centers, values, marker=marker, color=color, label=text)
    ax.set_xlabel("reconstructed rigidity [GV]")
    ax.set_ylabel(r"$W_{68}(m^2_{\rm reco}-m^2_{\rm true})$ [GeV$^2$]")
    ax.set_xticks(centers)
    ax.set_xticklabels(["2.5–3", "3–4", "4–5", "5–6"])
    ax.set_title("What actually limits the physical mass peak?")
    ax.legend(frameon=False)
    clean(ax)
    save(fig, "06_mass2_error_decomposition")

    # 7. Equal-prior probability calibration.
    fig, ax = plt.subplots()
    for method in ["raw_linear", "all_physics", "all_detector", "nuisance_control"]:
        rows = [row for row in calibration if row["method"] == method]
        rows.sort(key=lambda row: row["mean_prediction"])
        ax.plot([row["mean_prediction"] for row in rows],
                [row["observed_c13_fraction"] for row in rows],
                marker="o", color=METHOD_COLORS[method], label=DISPLAY[method])
    ax.plot([0, 1], [0, 1], color="#777777", linestyle=":", linewidth=1.2)
    ax.set_xlabel("mean predicted probability")
    ax.set_ylabel(r"observed $^{13}$C fraction")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_title("Test-set calibration under equal isotope priors")
    ax.legend(frameon=False)
    clean(ax)
    save(fig, "07_probability_calibration")

    # 8. Overall AUC with source-file bootstrap intervals.
    fig, ax = plt.subplots()
    methods = ["raw_linear", "raw_hybrid", "mass_proxy", "mass_energy", "all_physics", "all_detector", "nuisance_control"]
    positions = np.arange(len(methods))
    central = np.asarray([overall[method]["auc"] for method in methods])
    low = np.asarray([bootstrap["methods"][method]["auc"]["low95"] for method in methods])
    high = np.asarray([bootstrap["methods"][method]["auc"]["high95"] for method in methods])
    colors = [METHOD_COLORS[method] for method in methods]
    ax.barh(positions, central - 0.5, left=0.5, color=colors, alpha=0.88)
    ax.errorbar(central, positions, xerr=np.vstack([central - low, high - central]),
                fmt="none", ecolor=COLORS["black"], capsize=3, linewidth=1.2)
    ax.set_yticks(positions)
    ax.set_yticklabels([DISPLAY[method] for method in methods])
    ax.invert_yaxis()
    ax.set_xlim(0.5, 1.0)
    ax.set_xlabel("AUC (95% source-file bootstrap interval)")
    ax.set_title("Information-source ladder on independent files")
    clean(ax)
    save(fig, "08_auc_bootstrap_summary")

    print("wrote figures to", OUT)


if __name__ == "__main__":
    main()
