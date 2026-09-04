"""Generates README-facing result figures from already-committed CSVs --
no new measurement, just visualizing numbers that were previously only
presented as markdown tables. Four figures, two per side of the
synthetic/real split the README's Results section is built around:

1. `construction_progression_plot.png` -- synthetic ladder, CX vs. size,
   whole-graph (technique 2) -> zone decomposition -> its cost-capped
   refinement (technique 3, both stages), against the 500-CX target line.
2. `synthetic_nisq_feasibility_plot.png` -- the same ladder's whole-graph
   vs. cost-capped numbers at its hardest tested size (n_nodes=150),
   against published NISQ fidelity curves.
3. `real_network_comparison_plot.png` -- grouped bar chart, CX cost per
   construction, both real networks, against the 500-CX target line.
4. `real_nisq_feasibility_plot.png` -- the two real networks' three
   construction variants, against the same fidelity curves.

(Two earlier figures were removed as redundant once the README's own
sections were tightened: `ladder_cx_plot.png` -- its one line was already
reproduced exactly as construction_progression_plot.png's whole-graph
curve -- and `fixed_vs_log_tiecount_plot.png`, which existed only to
support a since-removed section. `scaling_log_ties_summary.csv` and the
script that produces it are kept regardless -- they remain the
evidentiary basis for `docs/scaling-ladder-and-decomposition.md`'s
section 1.)

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


def plot_construction_progression() -> None:
    whole = _rows("cost_aware_scaling_ladder_summary.csv")
    flat = _rows("decomposed_cost_aware_ladder_summary.csv")
    capped = _rows("cost_capped_decomposition_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, cond, title in ((axes[0], "short_log", "short-range ties, log growth"), (axes[1], "long_log", "long-range ties, log growth")):
        xs_w, ys_w = _series(whole, "condition", cond, "n_nodes", "cx_mean")
        xs_f, ys_f = _series(flat, "condition", cond, "n_nodes", "cx_mean")
        xs_c, ys_c = _series(capped, "condition", cond, "n_nodes", "cx_mean")
        ax.plot(xs_w, ys_w, marker="o", color="#C0392B", label="whole-graph, no decomposition (technique 2)")
        ax.plot(xs_f, ys_f, marker="s", color="#D9822B", label="+ zone decomposition (technique 3)")
        ax.plot(xs_c, ys_c, marker="^", color="#2E8B57", label="+ cost-capped refinement (technique 3)")
        ax.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
        ax.set_yscale("log")
        ax.set_xlabel("n_nodes")
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3, which="both")
    axes[0].set_ylabel("transpiled CX count (mean over seeds)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("Synthetic ladder: each stage attacks a different part of the cost", fontsize=12)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    out = RESULTS_DIR / "construction_progression_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def _fidelity_axes(ax, title: str) -> None:
    cx = np.logspace(1, 4.3, 300)
    for p, color, label in ((0.001, "#4472C4", "best-case trapped-ion (p=0.001/CX)"), (0.005, "#C0392B", "typical superconducting (p=0.005/CX)")):
        ax.plot(cx, (1 - p) ** cx, color=color, label=label, zorder=3)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("transpiled CX count")
    ax.set_ylabel(r"circuit survival, $fidelity \approx (1-p)^{N_{CX}}$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3, which="both")


def plot_synthetic_nisq_feasibility() -> None:
    """Vertical reference lines, not curve-placed markers -- same reasoning
    as the real-network version: a marker sitting ON one curve would
    misleadingly imply that curve's fidelity even where the other curve's
    value is what matters. Uses the ladder's hardest tested size
    (n_nodes=150) as the representative point for each condition."""
    whole = _rows("cost_aware_scaling_ladder_summary.csv")
    capped = _rows("cost_capped_decomposition_summary.csv")

    def cx_at_150(rows: list[dict], condition: str) -> float:
        return next(float(r["cx_mean"]) for r in rows if r["condition"] == condition and r["n_nodes"] == "150")

    fig, ax = plt.subplots(figsize=(9, 5.5))
    _fidelity_axes(ax, "Where the synthetic ladder lands relative to NISQ feasibility\n(hardest tested size, n_nodes=150)")

    lines = [
        ("short-range, whole-graph", cx_at_150(whole, "short_log"), "#D9822B", "dashdot"),
        ("short-range, cost-capped", cx_at_150(capped, "short_log"), "#2E8B57", "solid"),
        ("long-range, whole-graph", cx_at_150(whole, "long_log"), "#B0392B", "dashdot"),
        ("long-range, cost-capped", cx_at_150(capped, "long_log"), "#66C2A5", "solid"),
    ]
    for label, xval, color, ls in lines:
        ax.axvline(xval, color=color, ls=ls, lw=1.6, alpha=0.85, label=f"{label} ({int(round(xval))} CX)", zorder=2)

    ax.legend(fontsize=8, loc="lower left", ncol=1)
    fig.tight_layout()
    out = RESULTS_DIR / "synthetic_nisq_feasibility_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def plot_real_nisq_feasibility() -> None:
    real = _rows("real_networks_hierarchical_results.csv")

    def cx_of(network: str, method: str) -> float | None:
        for r in real:
            if r["network"] == network and r["method"] == method:
                return float(r["cx_count"])
        return None

    fig, ax = plt.subplots(figsize=(9, 5.5))
    _fidelity_axes(ax, "Where each real-network construction lands relative to NISQ feasibility\n(read a network's fidelity off either curve at its vertical line)")

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

    ax.legend(fontsize=8, loc="lower left", ncol=1)
    fig.tight_layout()
    out = RESULTS_DIR / "real_nisq_feasibility_plot.png"
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
    methods = [("exact_whole_graph", "exact whole-graph", "#4472C4"), ("decomposed", "zone decomposition", "#D9822B"), ("cost_capped", "cost-capped refinement", "#2E8B57")]
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
    ax.set_title("Real networks: exact vs. zone decomposition vs. cost-capped refinement")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    out = RESULTS_DIR / "real_network_comparison_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    plot_construction_progression()
    plot_synthetic_nisq_feasibility()
    plot_real_nisq_feasibility()
    plot_real_network_comparison()
