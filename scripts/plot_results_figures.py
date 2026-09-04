"""Generates README-facing result figures from already-committed CSVs --
no new measurement, just visualizing numbers that were previously only
presented as markdown tables. Five figures:

1. `fixed_vs_log_tiecount_plot.png` -- headline result 1's flat-tie-count
   assumption vs. the more realistic log-scaled alternative (does cost
   really decrease with scale? depends which you assume).
2. `ladder_cx_plot.png` -- the escalating ladder's two default conditions
   (short_log, long_log), cost-aware bounded-witness mixer, CX vs. size.
3. `construction_progression_plot.png` -- the three-stage fix, CX vs.
   size, whole-graph -> flat decomposition -> cost-capped decomposition,
   against the 500-CX target line.
4. `nisq_feasibility_plot.png` -- fidelity vs. CX count (published
   two-qubit gate error rates), with the two real networks' three
   construction variants plotted directly on the curve.
5. `real_network_comparison_plot.png` -- grouped bar chart, CX cost per
   construction, both real networks, against the 500-CX target line.

Reads only `results/*.csv`; writes only `results/*_plot.png`. Regenerate
with `python scripts/plot_results_figures.py`.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CX_THRESHOLD = 500


def _rows(name: str) -> list[dict]:
    with open(RESULTS_DIR / name, newline="") as f:
        return list(csv.DictReader(f))


def _series(rows: list[dict], key_col: str, key_val: str, x_col: str, y_col: str):
    pts = sorted(
        ((float(r[x_col]), float(r[y_col])) for r in rows if r[key_col] == key_val),
        key=lambda p: p[0],
    )
    xs, ys = zip(*pts)
    return list(xs), list(ys)


def plot_fixed_vs_log_tiecount() -> None:
    rows = _rows("scaling_log_ties_summary.csv")
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for cond, color, label in (
        ("fixed", "#4472C4", "fixed $k_{ties}=3$ (headline result 1's own assumption)"),
        ("log", "#C0392B", "log-scaled $k_{ties}(n) = round(1.43\\ln n)$ (real-data-calibrated)"),
    ):
        xs, ys = _series(rows, "condition", cond, "n_nodes", "cx_mean")
        ax.plot(xs, ys, marker="o", color=color, label=label)
    ax.set_xlabel("n_nodes")
    ax.set_ylabel("transpiled CX count (mean over seeds)")
    ax.set_title("Same exact construction, two tie-count growth assumptions:\ncost trend inverts depending which one is true")
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = RESULTS_DIR / "fixed_vs_log_tiecount_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def plot_ladder_cx() -> None:
    rows = _rows("cost_aware_scaling_ladder_summary.csv")
    fig, ax = plt.subplots(figsize=(7, 4.8))
    for cond, color, label in (
        ("short_log", "#4472C4", "short-range ties, log growth"),
        ("long_log", "#C0392B", "long-range ties, log growth"),
    ):
        cond_rows = [r for r in rows if r["condition"] == cond]
        cond_rows.sort(key=lambda r: float(r["n_nodes"]))
        xs = [float(r["n_nodes"]) for r in cond_rows]
        ys = [float(r["cx_mean"]) for r in cond_rows]
        ymin = [float(r["cx_min"]) for r in cond_rows]
        ymax = [float(r["cx_max"]) for r in cond_rows]
        ax.plot(xs, ys, marker="o", color=color, label=label)
        ax.fill_between(xs, ymin, ymax, color=color, alpha=0.15)
    ax.set_yscale("log")
    ax.set_xlabel("n_nodes")
    ax.set_ylabel("transpiled CX count (mean, min-max band)")
    ax.set_title("Cost-aware bounded-witness mixer, escalating realism ladder\n(whole-graph, no decomposition yet)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    out = RESULTS_DIR / "ladder_cx_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def plot_construction_progression() -> None:
    whole = _rows("cost_aware_scaling_ladder_summary.csv")
    flat = _rows("decomposed_cost_aware_ladder_summary.csv")
    capped = _rows("cost_capped_decomposition_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, cond, title in ((axes[0], "short_log", "short-range ties, log growth"), (axes[1], "long_log", "long-range ties, log growth")):
        xs_w, ys_w = _series(whole, "condition", cond, "n_nodes", "cx_mean")
        xs_f, ys_f = _series(flat, "condition", cond, "n_nodes", "cx_mean")
        xs_c, ys_c = _series(capped, "condition", cond, "n_nodes", "cx_mean")
        ax.plot(xs_w, ys_w, marker="o", color="#C0392B", label="whole-graph (technique 3+4)")
        ax.plot(xs_f, ys_f, marker="s", color="#D9822B", label="+ flat zone decomposition (technique 2)")
        ax.plot(xs_c, ys_c, marker="^", color="#2E8B57", label="+ cost-capped decomposition (technique 6)")
        ax.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
        ax.set_yscale("log")
        ax.set_xlabel("n_nodes")
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3, which="both")
    axes[0].set_ylabel("transpiled CX count (mean over seeds)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("The three-stage fix: each technique attacks a different part of the cost", fontsize=12)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    out = RESULTS_DIR / "construction_progression_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_nisq_feasibility() -> None:
    cx = np.logspace(1, 4.3, 300)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for p, color, label in ((0.001, "#4472C4", "best-case trapped-ion (p=0.001/CX)"), (0.005, "#C0392B", "typical superconducting (p=0.005/CX)")):
        ax.plot(cx, (1 - p) ** cx, color=color, label=label, zorder=3)

    real = _rows("real_networks_hierarchical_results.csv")

    def cx_of(network: str, method: str) -> float | None:
        for r in real:
            if r["network"] == network and r["method"] == method:
                return float(r["cx_count"])
        return None

    # Vertical reference lines, not curve-placed markers: a marker sitting
    # ON the trapped-ion curve would visually imply that specific fidelity
    # value even for the superconducting case, which isn't true -- a
    # vertical line at the CX count lets both curves' actual values be read
    # off directly, for either error rate.
    lines = [
        ("CIGRE MV exact", cx_of("CIGRE_MV", "exact_whole_graph"), "#888888", "dotted"),
        ("CIGRE MV decomposed", cx_of("CIGRE_MV", "decomposed"), "#D9822B", "dashdot"),
        ("CIGRE MV cost-capped", cx_of("CIGRE_MV", "cost_capped"), "#2E8B57", "solid"),
        ("IEEE33 exact (disconnected)", 96, "#B0392B", "dotted"),
        ("IEEE33 decomposed", cx_of("IEEE33", "decomposed"), "#F0B27A", "dashdot"),
        ("IEEE33 cost-capped", cx_of("IEEE33", "cost_capped"), "#66C2A5", "solid"),
    ]
    for label, xval, color, ls in lines:
        if xval is None:
            continue
        ax.axvline(xval, color=color, ls=ls, lw=1.6, alpha=0.85, label=f"{label} ({int(xval)} CX)", zorder=2)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("transpiled CX count")
    ax.set_ylabel(r"circuit survival, $fidelity \approx (1-p)^{N_{CX}}$")
    ax.set_title("Where each construction lands relative to NISQ feasibility\n(read a network's fidelity off either curve at its vertical line)")
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(fontsize=8, loc="lower left", ncol=1)
    fig.tight_layout()
    out = RESULTS_DIR / "nisq_feasibility_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def plot_real_network_comparison() -> None:
    real = _rows("real_networks_hierarchical_results.csv")

    def cx_of(network: str, method: str) -> float | None:
        for r in real:
            if r["network"] == network and r["method"] == method:
                return float(r["cx_count"])
        return None

    networks = ["CIGRE_MV", "IEEE33"]
    methods = [("exact_whole_graph", "exact whole-graph", "#4472C4"), ("decomposed", "flat decomposed", "#D9822B"), ("cost_capped", "cost-capped decomposed", "#2E8B57")]
    exact_vals = {"CIGRE_MV": cx_of("CIGRE_MV", "exact_whole_graph"), "IEEE33": 96}

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(networks))
    width = 0.25
    for i, (method, label, color) in enumerate(methods):
        vals = []
        for net in networks:
            v = exact_vals[net] if method == "exact_whole_graph" else cx_of(net, method)
            vals.append(v if v is not None else 0)
        bars = ax.bar(x + (i - 1) * width, vals, width, label=label, color=color)
        for b, v, net in zip(bars, vals, networks):
            note = ""
            if method == "exact_whole_graph" and net == "IEEE33":
                note = "\n(disconnected)"
            ax.text(b.get_x() + b.get_width() / 2, v, f"{int(v)}{note}", ha="center", va="bottom", fontsize=8)

    ax.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["CIGRE MV (15 bus, 3 ties)", "IEEE33 (33 bus, 5 ties)"])
    ax.set_ylabel("transpiled CX count")
    ax.set_title("Real networks: exact vs. decomposed vs. cost-capped construction")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    out = RESULTS_DIR / "real_network_comparison_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    plot_fixed_vs_log_tiecount()
    plot_ladder_cx()
    plot_construction_progression()
    plot_nisq_feasibility()
    plot_real_network_comparison()
