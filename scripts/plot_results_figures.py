"""Generates the result figures from already-committed CSVs -- no new
measurement, just visualizing numbers that were previously only
presented as markdown tables. Five figures; the README embeds three of
them (1, 4, 5) and the docs reference the other two:

Every figure shows the SAME TWO TIERS and nothing else: whole-graph with
no decomposition (the fault-tolerant tier) and cost-capped decomposition
(the NISQ tier). Zone decomposition on its own (technique 3a) and the
exact construction (technique 1) are deliberately absent everywhere --
3a is an intermediate step strictly dominated by 3b on both cost and
safety at every size tested, and technique 1's failure mode is dropping
candidates outright rather than getting expensive, which makes it look
artificially cheap on a cost axis and artificially perfect on a safety
axis. Plotting either invites a comparison between things that aren't
alternatives.

1. `construction_progression_plot.png` -- synthetic ladder, CX vs. size,
   both tiers, against the 500-CX target line. Long-range-tie condition
   only (the one that models real feeders; the short-range control stays
   in the docs and in figure 2).
2. `synthetic_mass_progression_plot.png` -- same ladder, same two tiers,
   mean feasible mass (safety) instead of CX.
3. `synthetic_nisq_feasibility_plot.png` -- both tiers at the ladder's
   hardest tested size (n_nodes=150), against published NISQ fidelity
   curves.
4. `real_network_comparison_plot.png` -- grouped bar chart, both tiers,
   both real networks, against the 500-CX target line.
5. `real_nisq_feasibility_plot.png` -- both tiers on both real networks,
   against the same fidelity curves.

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
    """Two tiers, two lines -- deliberately NOT three. Zone decomposition
    on its own (technique 3a) used to be plotted as a middle line, but it
    is an intermediate step toward the cost-capped construction, not a
    tier anyone would actually deploy: it is strictly worse than 3b on
    both axes (cost AND safety) at every size tested. Showing it invited
    readers to compare three things when only two are on offer."""
    # technique 2 = FIXED cost_alpha=0.01, this construction's validated
    # default -- NOT cost_aware_scaling_ladder_summary.csv, which
    # (per that script's own docstring) measures the ADAPTIVE variant.
    # Long-range ties only. The ladder's short-range condition is a
    # synthetic control (docs/scaling-ladder-and-decomposition.md keeps
    # it), not a claim about real feeders: real tie switches are
    # long-range by design (33-45% of network diameter), so the
    # README-facing figure shows only the condition that models them.
    # Per-seed rows, not the summary CSVs: the summaries carry mean/std or
    # mean/max only, and an evaluator reading a means-only line has no way
    # to see that decomposed variants reach 99% coefficient of variation
    # at some sizes (docs/scaling-ladder-and-decomposition.md §11). Min-max
    # bars over the 3 seeds make the spread visible on the headline figure.
    whole = _rows("fixed_alpha_ladder_results.csv")
    capped = _rows("cost_capped_decomposition_results.csv")

    def seed_stats(rows, cx_col):
        by_n: dict[int, list[float]] = {}
        for r in rows:
            if r["condition"] == "long_log":
                by_n.setdefault(int(r["n_nodes"]), []).append(float(r[cx_col]))
        xs = sorted(by_n)
        means = [sum(by_n[x]) / len(by_n[x]) for x in xs]
        lo = [m - min(by_n[x]) for m, x in zip(means, xs)]
        hi = [max(by_n[x]) - m for m, x in zip(means, xs)]
        return xs, means, [lo, hi]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    xs_w, ys_w, err_w = seed_stats(whole, "cx_count")
    xs_c, ys_c, err_c = seed_stats(capped, "total_cx")
    ax.errorbar(xs_w, ys_w, yerr=err_w, marker="o", capsize=3, color="#C0392B",
                label="whole-graph, no decomposition — fault-tolerant tier")
    ax.errorbar(xs_c, ys_c, yerr=err_c, marker="^", capsize=3, color="#2E8B57",
                label="cost-capped decomposition — NISQ tier")
    ax.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
    ax.set_yscale("log")
    ax.set_xlabel("network size (n_nodes)")
    ax.set_ylabel("transpiled CX count")
    ax.set_title("Synthetic feeders with real-topology tie statistics: the two tiers\n(mean over 3 seeds, bars = min–max)", fontsize=11)
    ax.grid(True, alpha=0.3, which="both")
    # Legend below the axes: every inside corner is occupied by one of
    # the two curves (red rises top-left, plateaus top-right; green
    # rises bottom-right), so any in-axes placement hides data.
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    out = RESULTS_DIR / "construction_progression_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_synthetic_mass_progression() -> None:
    """Companion to plot_construction_progression -- same two tiers,
    same ladder, but feasible mass (safety) instead of CX
    (cost). Technique 1 (exact) is deliberately absent from both: its
    failure mode is dropped candidates / disconnection, not leakage --
    `measure_exact` reports mean_feasible_mass=1.0 by construction
    whenever it runs at all, so a leakage-axis plot would show it as a
    flat, misleadingly perfect line that hides its actual failure mode
    (already told honestly via `dropped_candidates`/`fully_connected`
    elsewhere). Exact construction also doesn't scale to n_nodes=150
    long-range at all -- brute-force enumeration is the reason this
    repo's exact-construction sweeps stay at small sizes everywhere else."""
    # technique 2 = FIXED cost_alpha=0.01, this construction's validated
    # default -- NOT cost_aware_scaling_ladder_summary.csv, which
    # (per that script's own docstring) measures the ADAPTIVE variant.
    whole = _rows("fixed_alpha_ladder_summary.csv")
    capped = _rows("cost_capped_decomposition_summary.csv")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, cond, title in ((axes[0], "short_log", "short-range ties, log growth"), (axes[1], "long_log", "long-range ties, log growth")):
        xs_w, ys_w = _series(whole, "condition", cond, "n_nodes", "mean_feasible_mass_mean")
        xs_c, ys_c = _series(capped, "condition", cond, "n_nodes", "mean_feasible_mass_mean")
        ax.plot(xs_w, ys_w, marker="o", color="#C0392B", label="whole-graph, no decomposition — fault-tolerant tier")
        ax.plot(xs_c, ys_c, marker="^", color="#2E8B57", label="cost-capped decomposition — NISQ tier")
        ax.axhline(1.0, color="black", ls="dashed", lw=1, label="perfect (no leakage)")
        ax.set_xlabel("n_nodes")
        ax.set_title(title, fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0.85, 1.02)
    axes[0].set_ylabel("mean feasible mass (1.0 = no leakage)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.08))
    fig.suptitle("The same two tiers, measured for safety instead of cost", fontsize=12)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    out = RESULTS_DIR / "synthetic_mass_progression_plot.png"
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
    # technique 2 = FIXED cost_alpha=0.01, this construction's validated
    # default -- NOT cost_aware_scaling_ladder_summary.csv, which
    # (per that script's own docstring) measures the ADAPTIVE variant.
    whole = _rows("fixed_alpha_ladder_summary.csv")
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
    """Uses real_networks_hierarchical_summary.csv (mean over
    SEEDS_PER_NETWORK seeds for decomposed/cost_capped -- see that
    script's own docstring for why averaging matters here: CIGRE MV's
    flat-decomposed result varies 3,668-9,178 CX depending on seed, even
    though its zones use fully deterministic exact tree enumeration
    (the randomness is in the witness search's restart order, not the
    tree set)."""
    real = _rows("real_networks_hierarchical_summary.csv")

    def cx_of(network: str, method: str) -> float | None:
        for r in real:
            if r["network"] == network and r["method"] == method:
                return float(r["cx_mean"])
        return None

    fig, ax = plt.subplots(figsize=(9, 5.5))
    _fidelity_axes(ax, "Where each real network's two tiers land relative to NISQ feasibility\n(read a network's fidelity off either curve at its vertical line; mean over 5 seeds)")

    # Two tiers only, matching the synthetic figures: whole-graph
    # (fault-tolerant) and cost-capped decomposition (NISQ). Zone
    # decomposition alone and the exact construction are both omitted for
    # the same reason they are omitted from the synthetic plots -- neither
    # is a tier anyone would deploy (3a is strictly dominated by 3b; the
    # exact construction's IEEE33 result is cheap only because it is
    # mostly incomplete, 573/597 candidates dropped and disconnected).
    lines = [
        ("CIGRE MV whole-graph (fault-tolerant)", cx_of("CIGRE_MV", "truncated_whole_graph"), "#C0392B", "dashdot"),
        ("CIGRE MV cost-capped (NISQ)", cx_of("CIGRE_MV", "cost_capped"), "#2E8B57", "solid"),
        ("IEEE33 whole-graph (fault-tolerant)", cx_of("IEEE33", "truncated_whole_graph"), "#E8746A", "dashdot"),
        ("IEEE33 cost-capped (NISQ)", cx_of("IEEE33", "cost_capped"), "#66C2A5", "solid"),
    ]
    for label, xval, color, ls in lines:
        if xval is None:
            continue
        ax.axvline(xval, color=color, ls=ls, lw=1.6, alpha=0.85, label=f"{label} ({int(round(xval))} CX)", zorder=2)

    ax.legend(fontsize=8, loc="lower left", ncol=1)
    fig.tight_layout()
    out = RESULTS_DIR / "real_nisq_feasibility_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def plot_real_network_comparison() -> None:
    """Two tiers per network, error bars = min-max across
    SEEDS_PER_NETWORK seeds. Same two-line restriction as the synthetic
    figures (see plot_construction_progression's docstring): zone
    decomposition alone and the exact construction are omitted because
    neither is a deployable tier."""
    real = _rows("real_networks_hierarchical_summary.csv")

    def stats_of(network: str, method: str) -> tuple[float, float, float] | None:
        for r in real:
            if r["network"] == network and r["method"] == method:
                return float(r["cx_mean"]), float(r["cx_min"]), float(r["cx_max"])
        return None

    networks = ["CIGRE_MV", "IEEE33"]
    methods = [
        ("truncated_whole_graph", "whole-graph, no decomposition — fault-tolerant tier", "#C0392B"),
        ("cost_capped", "cost-capped decomposition — NISQ tier", "#2E8B57"),
    ]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(networks))
    width = 0.32
    for i, (method, label, color) in enumerate(methods):
        vals, err_lo, err_hi = [], [], []
        for net in networks:
            v, mn, mx = stats_of(net, method)
            vals.append(v)
            err_lo.append(v - mn)
            err_hi.append(mx - v)
        bars = ax.bar(x + (i - 0.5) * width, vals, width, label=label, color=color,
                       yerr=[err_lo, err_hi], capsize=4, ecolor="black")
        for b, v, hi in zip(bars, vals, err_hi):
            # label above the error bar's upper cap, not on top of it
            ax.text(b.get_x() + b.get_width() / 2, (v + hi) * 1.08, f"{int(round(v))}",
                    ha="center", va="bottom", fontsize=8)

    ax.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(["CIGRE MV (15 bus, 3 ties)", "IEEE33 (33 bus, 5 ties)"])
    ax.set_ylabel("transpiled CX count")
    ax.set_title("Real networks: the two deployment tiers")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis="y", which="both")
    fig.tight_layout()
    out = RESULTS_DIR / "real_network_comparison_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    plot_construction_progression()
    plot_synthetic_mass_progression()
    plot_synthetic_nisq_feasibility()
    plot_real_nisq_feasibility()
    plot_real_network_comparison()
