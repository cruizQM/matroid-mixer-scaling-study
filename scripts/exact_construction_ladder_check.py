"""Checks a claim `run_cost_aware_scaling_ladder.py`'s own figures make by
omission, not measurement: that technique 1 (the exact construction)
can't be added to the escalating ladder's CX/mass progression plots
without misleading on both axes at once. Builds it directly on the SAME
graphs (`CONDITIONS`, same generator + k_ties functions + seeds) at the
two sizes where brute-force spanning-tree enumeration stays tractable
(`n_nodes` 10 and 30 -- 60 already exceeds a 30s budget on this hardware
for both conditions, confirmed by direct timing before writing this
script, not assumed).

## What this finds, and why it matters

On `long_log` at `n_nodes=30`, the exact construction drops 206-278 of
561 candidate exchanges (37-50%) and is fully disconnected -- yet reports
only 72-76 CX, far BELOW `run_cost_aware_scaling_ladder.py`'s own
whole-graph number at the same size (17,225 CX) -- because a dropped
candidate costs nothing, not because the construction is cheap. Worse,
`mean_feasible_mass` reports exactly 1.0 on every single row, including
that catastrophically incomplete one: leakage is only measured among
terms that exist, so a dropped candidate isn't leaky, it's just absent.
Technique 1 would look artificially CHEAP on a cost plot and artificially
PERFECT on a safety plot, in exactly the two places this document's own
progression figures compare techniques 2/3a/3b -- confirming, with real
numbers rather than just the general argument already made in the
Techniques section, why it's correctly excluded from both.

Writes results/exact_construction_ladder_check_results.csv and two
figures: results/exact_construction_ladder_check_cx_plot.png,
results/exact_construction_ladder_check_mass_plot.png.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import numpy as np
from qiskit import transpile

from graphs import enumerate_spanning_trees
from leakage_trace import final_feasible_mass
from measure import TRANSPILE_BASIS
from mixer import build_matroid_mixer, mixer_circuit
from random_trees import random_spanning_tree
from run_cost_aware_scaling_ladder import CONDITIONS

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "exact_construction_ladder_check_results.csv"

BETA = 0.37
N_NODES_TRACTABLE = [10, 30]  # 60+ exceeds a 30s brute-force-enumeration budget, both conditions
SEEDS_PER_SIZE = 3
N_STARTING_TREES = 50


def measure_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)

    trees = enumerate_spanning_trees(graph)
    construction = build_matroid_mixer(graph, trees)
    qc = mixer_circuit(construction, beta=BETA)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    cx = tqc.count_ops().get("cx", 0)

    rng = np.random.default_rng(seed * 1000 + 2)
    masses = [
        final_feasible_mass(graph, construction, random_spanning_tree(graph, rng), BETA, sparse=True)
        for _ in range(N_STARTING_TREES)
    ]
    n_candidates = graph.n_edges * (graph.n_edges - 1) // 2

    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits": graph.n_edges, "n_trees": len(trees), "n_candidates": n_candidates,
        "n_terms": len(construction.terms), "dropped_candidates": construction.dropped_candidates,
        "fully_connected": construction.fully_connected,
        "cx_count": cx, "depth": tqc.depth(),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
    }


def make_plots(rows: list[dict]) -> None:
    def series(cond: str, n_nodes_list: list[int], field: str) -> tuple[list[int], list[float]]:
        out = []
        for n in n_nodes_list:
            group = [r[field] for r in rows if r["condition"] == cond and r["n_nodes"] == n]
            out.append(statistics.mean(group))
        return n_nodes_list, out

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, cond, title in ((axes[0], "short_log", "short-range"), (axes[1], "long_log", "long-range")):
        xs, ys = series(cond, N_NODES_TRACTABLE, "cx_count")
        ax.plot(xs, ys, marker="x", color="#7D3C98", markersize=10, label="exact (technique 1)")
        ax.set_yscale("log")
        ax.set_xlabel("n_nodes")
        ax.set_title(f"{title}, log growth", fontsize=11)
        ax.grid(True, alpha=0.3, which="both")
        ax.set_xlim(0, 40)
    axes[0].set_ylabel("transpiled CX count (mean over seeds)")
    fig.suptitle("Exact construction (technique 1) on the ladder's own graphs: CX drops when candidates get dropped", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / "exact_construction_ladder_check_cx_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for ax, cond, title in ((axes[0], "short_log", "short-range"), (axes[1], "long_log", "long-range")):
        xs, ys = series(cond, N_NODES_TRACTABLE, "mean_feasible_mass")
        ax.plot(xs, ys, marker="x", color="#7D3C98", markersize=10, label="exact (technique 1)")
        ax.axhline(1.0, color="black", ls="dashed", lw=1)
        ax.set_ylim(0.85, 1.02)
        ax.set_xlabel("n_nodes")
        ax.set_title(f"{title}, log growth", fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, 40)
    axes[0].set_ylabel("mean feasible mass (1.0 = no leakage)")
    fig.suptitle("Same data, safety axis: looks PERFECT even where 37-50% of candidates were dropped", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    out = RESULTS_DIR / "exact_construction_ladder_check_mass_plot.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in ["short_log", "long_log"]:
        for n_nodes in N_NODES_TRACTABLE:
            for seed in range(SEEDS_PER_SIZE):
                row = measure_one(condition, n_nodes, seed)
                rows.append(row)
                print(f"  {row}", flush=True)

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    make_plots(rows)


if __name__ == "__main__":
    main()
