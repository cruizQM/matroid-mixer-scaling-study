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

Every synthetic figure shows the long-range-tie condition only: real
ties are long-range by design, and the short-range condition is a
synthetic control reported in the docs' tables, not drawn.

1. `construction_progression_plot.png` -- synthetic ladder, CX vs. size,
   both tiers, against the 500-CX target line. In the README.
2. `synthetic_mass_progression_plot.png` -- same ladder, same two tiers,
   mean feasible mass (safety) instead of CX. Docs only.
3. `real_networks_plot.png` -- both real networks, two panels: CX per
   tier with min-max bars, and the same numbers on published NISQ
   fidelity curves. In the README.

For the hard-instance case study (question 2), all in the README:

4. `hard_instance_proof_time_plot.png` -- CP-SAT proof time vs. gadget
   size, timeouts drawn as a wall.
5. `mps_convergence_plot.png` -- MPS accuracy vs. exact where
   computable, and simulation time vs. size at fixed bond dimension.
6. `large_b_grid_plot.png` -- the two dials (m, B) as a grid, one panel
   per classical route, DP cells carrying their state counts.

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
    # Long-range ties only, matching plot_construction_progression: the
    # short-range condition is a synthetic control, not a claim about
    # real feeders, and is reported in the docs' tables rather than drawn.
    whole = _rows("fixed_alpha_ladder_summary.csv")
    capped = _rows("cost_capped_decomposition_summary.csv")

    fig, ax = plt.subplots(figsize=(8, 4.8))
    xs_w, ys_w = _series(whole, "condition", "long_log", "n_nodes", "mean_feasible_mass_mean")
    xs_c, ys_c = _series(capped, "condition", "long_log", "n_nodes", "mean_feasible_mass_mean")
    ax.plot(xs_w, ys_w, marker="o", color="#C0392B", label="whole-graph, no decomposition — fault-tolerant tier")
    ax.plot(xs_c, ys_c, marker="^", color="#2E8B57", label="cost-capped decomposition — NISQ tier")
    ax.axhline(1.0, color="black", ls="dashed", lw=1, label="perfect (no leakage)")
    ax.set_xlabel("network size (n_nodes)")
    ax.set_ylabel("mean feasible mass (1.0 = no leakage)")
    ax.set_title("The same two tiers, measured for safety instead of cost", fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.85, 1.02)
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
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


def plot_real_networks() -> None:
    """The two real networks on one figure, two panels: left, transpiled
    CX per tier with min-max bars over SEEDS_PER_NETWORK seeds; right,
    where those same numbers land on published NISQ fidelity curves,
    drawn as vertical reference lines (a marker sitting ON one curve
    would misleadingly imply that curve's fidelity where the other's
    value is what matters). Replaces two separate figures that said the
    same four numbers twice at full width.

    Two tiers only, matching the synthetic figures: whole-graph
    (fault-tolerant) and cost-capped decomposition (NISQ). Zone
    decomposition alone and the exact construction are omitted for the
    same reason as everywhere else -- neither is a tier anyone would
    deploy (3a is dominated by 3b; the exact construction's IEEE33 result
    is cheap only because it is mostly incomplete, 573/597 candidates
    dropped and disconnected). Means over 5 seeds -- CIGRE MV's
    whole-graph tier varies 274-432 CX with seed."""
    real = _rows("real_networks_hierarchical_summary.csv")

    def stats_of(network: str, method: str) -> tuple[float, float, float]:
        r = next(r for r in real if r["network"] == network and r["method"] == method)
        return float(r["cx_mean"]), float(r["cx_min"]), float(r["cx_max"])

    networks = ["CIGRE_MV", "IEEE33"]
    tiers = [
        ("truncated_whole_graph", "whole-graph, no decomposition — fault-tolerant tier", "#C0392B", "#E8746A", "dashdot"),
        ("cost_capped", "cost-capped decomposition — NISQ tier", "#2E8B57", "#66C2A5", "solid"),
    ]

    fig, (ax_bar, ax_fid) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1, 1.25]})

    # --- left: bars
    x = np.arange(len(networks))
    width = 0.32
    for i, (method, label, color, _, _) in enumerate(tiers):
        vals, err_lo, err_hi = [], [], []
        for net in networks:
            v, mn, mx = stats_of(net, method)
            vals.append(v)
            err_lo.append(v - mn)
            err_hi.append(mx - v)
        bars = ax_bar.bar(x + (i - 0.5) * width, vals, width, label=label, color=color,
                          yerr=[err_lo, err_hi], capsize=4, ecolor="black")
        for b, v, hi in zip(bars, vals, err_hi):
            ax_bar.text(b.get_x() + b.get_width() / 2, (v + hi) * 1.08, f"{int(round(v))}",
                        ha="center", va="bottom", fontsize=8)
    ax_bar.axhline(CX_THRESHOLD, color="black", ls="dashed", lw=1, label=f"{CX_THRESHOLD} CX target")
    ax_bar.set_yscale("log")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(["CIGRE MV\n(15 bus, 3 ties)", "IEEE33\n(33 bus, 5 ties)"])
    ax_bar.set_ylabel("transpiled CX count")
    ax_bar.set_title("Circuit cost per tier (mean over 5 seeds, bars = min–max)", fontsize=10.5)
    ax_bar.grid(True, alpha=0.3, axis="y", which="both")

    # --- right: the same numbers on the fidelity curves
    _fidelity_axes(ax_fid, "Where those costs land on today's hardware")
    for net, net_label in (("CIGRE_MV", "CIGRE MV"), ("IEEE33", "IEEE33")):
        for method, _, color, light, ls in tiers:
            v, _, _ = stats_of(net, method)
            c = color if net == "CIGRE_MV" else light
            tier_word = "whole-graph" if method == "truncated_whole_graph" else "cost-capped"
            ax_fid.axvline(v, color=c, ls=ls, lw=1.6, alpha=0.85, label=f"{net_label} {tier_word} ({int(round(v))} CX)", zorder=2)
    ax_fid.legend(fontsize=7.5, loc="lower left")

    handles, labels = ax_bar.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Real networks: the two deployment tiers", fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95))
    out = RESULTS_DIR / "real_networks_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_hard_instance_proof_time() -> None:
    """Question 2's headline: CP-SAT's time to PROVE optimality on the
    hard-instance gadget vs. instance size, from
    results/hard_instance_hardness_sweep.csv. Timeouts (status != OPTIMAL)
    are drawn hollow at the cap so the curve's end is visibly a wall, not
    a data point. NO-instances (no perfect partition exists) are the
    exhaustive direction and the one that matters."""
    rows = _rows("hard_instance_hardness_sweep.csv")
    cap = 1800.0
    fig, ax = plt.subplots(figsize=(8, 4.6))
    for suffix, color, marker, label in (("_yes", "#4472C4", "o", "perfect partition exists (YES)"),
                                         ("_no", "#C0392B", "s", "no perfect partition (NO) — solver must exhaust")):
        pts = sorted((int(r["n_nodes"]), float(r["wall_clock_s"]), r["status"] == "OPTIMAL")
                     for r in rows if r["label"].endswith(suffix))
        xs = [p[0] for p in pts]
        ys = [p[1] for p in pts]
        ax.plot(xs, ys, color=color, marker=marker, label=label)
        for x, y, proved in pts:
            if not proved:
                ax.plot([x], [y], marker=marker, markersize=12, markerfacecolor="white", markeredgecolor=color,
                        markeredgewidth=2, linestyle="none", zorder=4)
                ax.annotate("unproved\nat cap", (x, y), textcoords="offset points", xytext=(-54, -26),
                            fontsize=8, color=color)
    ax.axhline(cap, color="black", ls="dashed", lw=1, label="30-minute cap")
    ax.set_yscale("log")
    ax.set_xlabel("gadget size (nodes)")
    ax.set_ylabel("CP-SAT time to prove optimality (s)")
    ax.set_title("A general-purpose solver gives up early on the hard instance", fontsize=11)
    ax.grid(True, alpha=0.3, which="both")
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    out = RESULTS_DIR / "hard_instance_proof_time_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_mps_convergence() -> None:
    """Can a classical tensor-network (MPS) simulator stand in for the
    hard-instance circuit? Two panels, one per half of the finding, from
    results/mps_scaling_check.csv (p=2, fixed representative angles).

    Left -- accuracy, only where it can be judged: the sizes with an
    exact statevector reference (m=2 and m=3). Relative error vs. bond
    dimension, in percent, so both lines share a meaningful axis. m=2
    converges; m=3 plateaus ~8% off and does not move.

    Right -- simulation time: wall-clock at a fixed bond dimension (64)
    vs. size, against a dashed baseline of what the time would be if it
    grew only with the circuit's own gate count (computed from the
    circuit's structure, scaled to match at the smallest size). The gap
    between the measured line and the baseline is the price of
    representing entanglement, not circuit size. The m=7 run that
    exceeded its budget is drawn hollow at the time it had consumed when
    discarded -- a lower bound. Deliberately never called "cost": in
    this repo cost means circuit cost (CX count), and this axis is
    seconds.

    An earlier single-panel version plotted expected cost / optimum for
    every m on one axis; that number differs per m even for the exact
    answer, so the lines were not comparable and, for m>=4, there was
    no reference to converge to -- unreadable."""
    from math import comb

    rows = _rows("mps_scaling_check.csv")
    fig, (ax_err, ax_cost) = plt.subplots(1, 2, figsize=(12, 4.6))

    # --- left: relative error where exact is known
    for m, color in ((2, "#4472C4"), (3, "#C0392B")):
        exact = next(float(r["expected_cost"]) for r in rows
                     if int(r["m"]) == m and r["method"] == "statevector" and r["expected_cost"])
        pts = sorted((int(r["bond_dim"]), 100.0 * abs(float(r["expected_cost"]) - exact) / exact)
                     for r in rows if int(r["m"]) == m and r["method"] == "mps" and r["expected_cost"])
        xs, ys = zip(*pts)
        ax_err.plot(xs, ys, marker="o", color=color, label=f"m={m} ({3 * m * m} qubits)")
        ax_err.annotate(f"{ys[-1]:.1f}%", (xs[-1], ys[-1]), textcoords="offset points", xytext=(6, 0),
                        fontsize=8.5, color=color, va="center")
    ax_err.set_xscale("log", base=2)
    ax_err.set_xlabel("MPS bond dimension")
    ax_err.set_ylabel("error vs. exact statevector (%)")
    ax_err.set_title("Accuracy, where exact is computable", fontsize=11)
    ax_err.grid(True, alpha=0.3, which="both")
    ax_err.legend(fontsize=8.5)

    # --- right: wall-clock at bond dim 64 vs. size, against a circuit-size baseline
    def two_qubit_gates(m: int, p: int = 2) -> int:
        k = 3 * m
        return p * (2 * comb(k, 2) * m + 4 * k * comb(m, 2))  # rzz -> 2 CX, exchange_swap -> 4 CX

    bd64 = sorted((int(r["m"]), float(r["wall_clock_s"]), bool(r["error"]))
                  for r in rows if r["method"] == "mps" and int(r["bond_dim"]) == 64)
    ms = [b[0] for b in bd64]
    ts = [b[1] for b in bd64]
    ok = [(m, t) for m, t, err in bd64 if not err]
    bad = [(m, t) for m, t, err in bd64 if err]
    ax_cost.plot([m for m, _ in ok], [t for _, t in ok], marker="o", color="#C0392B", label="measured, bond dimension 64")
    for m, t in bad:
        ax_cost.plot([m], [t], marker="o", markersize=11, markerfacecolor="white", markeredgecolor="#C0392B",
                     markeredgewidth=2, linestyle="none")
        ax_cost.annotate("over budget\n(lower bound)", (m, t), textcoords="offset points", xytext=(-62, -22),
                         fontsize=8, color="#C0392B")
    m0, t0 = ok[1] if len(ok) > 1 else ok[0]  # anchor the baseline at m=3, the first non-trivial size
    baseline = [t0 * two_qubit_gates(m) / two_qubit_gates(m0) for m in ms]
    ax_cost.plot(ms, baseline, ls="dashed", color="grey", label="if time grew only with circuit size")
    ax_cost.set_yscale("log")
    ax_cost.set_xlabel("m  (qubits = 3m²)")
    ax_cost.set_ylabel("MPS simulation time (s)")
    ax_cost.set_title("Simulation time at fixed bond dimension (64)", fontsize=11)
    ax_cost.grid(True, alpha=0.3, which="both")
    ax_cost.legend(fontsize=8.5, loc="upper left")

    fig.suptitle("MPS on the hard-instance circuit (p=2): accurate only at the smallest size, and increasingly slow",
                 fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = RESULTS_DIR / "mps_convergence_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


def plot_large_b_grid() -> None:
    """The two dials of the hard instance on one figure: m (rows, sets
    the size -- qubits = 3m^2) and B (columns, sets the magnitude of the
    weights and nothing else). One panel per classical route, from
    results/large_b_hardness_sweep.csv; a cell is green if that route
    succeeded (DP finished / CP-SAT proved) and red if it failed (timed
    out / unproved at the cap), labeled with the time. DP cells also
    carry the number of partial-sum states explored, from
    results/large_b_dp_states.csv -- the mechanism: B doesn't grow the
    instance, it grows the DP's memo table, because random large
    weights produce sums that stop colliding."""
    sweep = _rows("large_b_hardness_sweep.csv")
    states_rows = {(int(r["m"]), int(r["B"])): r for r in _rows("large_b_dp_states.csv")} \
        if (RESULTS_DIR / "large_b_dp_states.csv").exists() else {}
    ms = sorted({int(r["m"]) for r in sweep})
    bs = sorted({int(r["B"]) for r in sweep})
    cell = {(int(r["m"]), int(r["B"])): r for r in sweep}

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    ok_color, fail_color = "#CDEBD6", "#F4C7C3"
    for ax, title, which in ((axes[0], "Structure-aware route: the exact DP", "dp"),
                             (axes[1], "General-purpose route: CP-SAT on the natural MILP", "cpsat")):
        for i, m in enumerate(ms):
            for j, B in enumerate(bs):
                r = cell.get((m, B))
                if r is None:
                    ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="white", edgecolor="#999999"))
                    ax.text(j + 0.5, i + 0.5, "not run", ha="center", va="center", fontsize=8, color="#999999")
                    continue
                if which == "dp":
                    failed = r["dp_timed_out"] == "True"
                    t = float(r["dp_time_s"])
                    line1 = f">{t / 60:.0f} min, timed out" if failed else f"{t:.2f}s"
                    s = states_rows.get((m, B))
                    line2 = ""
                    if s is not None:
                        n = int(s["dp_states"])
                        line2 = f"\n{n:,}{'+' if s['dp_timed_out'] == 'True' else ''} states"
                    text = line1 + line2
                else:
                    failed = r["cpsat_proved_optimal"] != "True"
                    t = float(r["cpsat_time_s"])
                    text = f">{t / 60:.0f} min, unproved" if failed else f"{t:.1f}s, proved"
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=fail_color if failed else ok_color,
                                           edgecolor="#999999"))
                ax.text(j + 0.5, i + 0.5, text, ha="center", va="center", fontsize=8.5)
        ax.set_xlim(0, len(bs))
        ax.set_ylim(0, len(ms))
        ax.set_xticks([j + 0.5 for j in range(len(bs))])
        ax.set_xticklabels([f"B = {B:,}" for B in bs])
        ax.set_yticks([i + 0.5 for i in range(len(ms))])
        ax.set_yticklabels([f"m = {m}\n({3 * m * m} qubits)" for m in ms])
        ax.set_xlabel("B — magnitude of the weights (graph and circuit unchanged)")
        if ax is axes[0]:
            ax.set_ylabel("m — size of the instance")
        ax.set_title(title, fontsize=10.5)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, facecolor=ok_color, edgecolor="#999999", label="route succeeded"),
               plt.Rectangle((0, 0), 1, 1, facecolor=fail_color, edgecolor="#999999", label="route failed at its cap")]
    fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("Two dials: m makes the instance bigger, B makes it harder without making it bigger", fontsize=11.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.94))
    out = RESULTS_DIR / "large_b_grid_plot.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out}")


if __name__ == "__main__":
    plot_construction_progression()
    plot_synthetic_mass_progression()
    plot_real_networks()
    plot_hard_instance_proof_time()
    plot_mps_convergence()
    plot_large_b_grid()
