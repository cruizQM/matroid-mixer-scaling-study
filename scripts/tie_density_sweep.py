"""Tests whether tie DENSITY alone -- holding network size and the
nearest-tie (short-range) selection rule fixed, just adding more ties --
breaks the whole-graph witness-conditioned mixer, independent of the
long-range-tie failure mode `docs/circuit-validity.md` already documents.

## Why this is a different question from what's already documented

`docs/circuit-validity.md`'s real-feeder section traces the whole-graph
witness blowup to tie RANGE: real tie switches connect distant parts of a
feeder, giving long fundamental cycles, and `generate_feeder_graph_long_range_ties`
reproduces that specific failure mode at controllable size (see
`run_decomposition_scaling_study.py`). That leaves open a different
question: does `generate_feeder_graph` (the repo's own nearest-tie model,
whose SHORT fundamental cycles are exactly why the module docstring in
`mixer.py` says the witness bound "never exceeded 2 qubits" on it) stay
safe as tie density grows, with tie RANGE never changing? If witness size
tracks density on its own, that's a second, independent way the
bounded-witness assumption breaks -- one zone decomposition doesn't
obviously fix, since it targets long-range ties specifically, not density.

## What this measures

Fixes `n_nodes` (small enough that `enumerate_spanning_trees` -- brute
force over C(n_edges, n_nodes-1) subsets -- stays tractable across the
whole sweep) and increases `k_ties` from sparse to dense, using
`generate_feeder_graph` (nearest-tie, unchanged) at each point. Reports,
per instance: the max witness size `build_matroid_mixer`'s brute-force-only
search (`prefer_structural=False`, matching this repo's own default) finds
across all selected terms, and `dropped_candidates` -- exchange pairs that
connect new spanning trees but for which no witness was found within
`MAX_WITNESS_SEARCH_SIZE`, meaning the exact construction cannot represent
that move at all without raising the cap (which itself gets combinatorially
expensive -- see `mixer.find_witness_set`'s own docstring).

Writes results/tie_density_sweep_results.csv (per k_ties, per seed) and
results/tie_density_sweep_summary.csv (per k_ties, mean/median/max), plus
results/tie_density_sweep_plot.png.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from graphs import enumerate_spanning_trees, generate_feeder_graph
from mixer import build_matroid_mixer

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "tie_density_sweep_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "tie_density_sweep_summary.csv"
PLOT_PATH = RESULTS_DIR / "tie_density_sweep_plot.png"

N_NODES = 8
K_TIES_RANGE = [2, 4, 6, 8, 10]
SEEDS_PER_K = 4


def run_one(k_ties: int, seed: int) -> dict:
    t0 = time.perf_counter()
    graph = generate_feeder_graph(N_NODES, k_ties, seed)
    trees = enumerate_spanning_trees(graph)
    construction = build_matroid_mixer(graph, trees)
    elapsed = time.perf_counter() - t0

    witness_sizes = [t.control_count for t in construction.terms]
    return {
        "k_ties": k_ties,
        "seed": seed,
        "n_nodes": N_NODES,
        "n_qubits": graph.n_edges,
        "n_spanning_trees": len(trees),
        "n_candidates": graph.n_edges * (graph.n_edges - 1) // 2,
        "n_terms": len(construction.terms),
        "max_witness_size": max(witness_sizes, default=0),
        "mean_witness_size": round(statistics.mean(witness_sizes), 2) if witness_sizes else 0.0,
        "dropped_candidates": construction.dropped_candidates,
        "fully_connected": construction.fully_connected,
        "elapsed_s": round(elapsed, 2),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for k_ties in K_TIES_RANGE:
        for seed in range(SEEDS_PER_K):
            row = run_one(k_ties, seed)
            rows.append(row)
            print(f"  {row}", flush=True)

    fieldnames = list(rows[0].keys())
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {CSV_PATH}")

    def stats(group, key):
        vals = [r[key] for r in group]
        return {
            "mean": round(statistics.mean(vals), 2),
            "median": round(statistics.median(vals), 2),
            "max": max(vals),
        }

    summary_rows = []
    for k_ties in K_TIES_RANGE:
        group = [r for r in rows if r["k_ties"] == k_ties]
        row = {"k_ties": k_ties, "n_qubits": group[0]["n_qubits"], "n_seeds": len(group)}
        for key in ("max_witness_size", "n_spanning_trees", "dropped_candidates", "n_terms"):
            s = stats(group, key)
            row[f"{key}_mean"] = s["mean"]
            row[f"{key}_median"] = s["median"]
            row[f"{key}_max"] = s["max"]
        row["any_disconnected"] = any(not r["fully_connected"] for r in group)
        summary_rows.append(row)

    summary_fields = list(summary_rows[0].keys())
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")

    try:
        import matplotlib.pyplot as plt

        ks = [r["k_ties"] for r in summary_rows]
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        ax1.plot(ks, [r["max_witness_size_mean"] for r in summary_rows], "o-", color="crimson", label="max witness size (mean)")
        ax1.plot(ks, [r["max_witness_size_max"] for r in summary_rows], "o--", color="crimson", alpha=0.4, label="max witness size (max)")
        ax1.set_xlabel(f"k_ties (n_nodes={N_NODES} fixed, nearest-tie/short-range model)")
        ax1.set_ylabel("qubits")
        ax1.set_title("Tie density vs. required witness size")
        ax1.legend(fontsize=8)
        ax1.grid(alpha=0.3)

        ax2.plot(ks, [r["dropped_candidates_mean"] for r in summary_rows], "s-", color="steelblue", label="dropped candidates (mean)")
        ax2.set_xlabel(f"k_ties (n_nodes={N_NODES} fixed)")
        ax2.set_ylabel("candidates with no witness found within cap")
        ax2.set_title("Tie density vs. exact-construction coverage")
        ax2.legend(fontsize=8)
        ax2.grid(alpha=0.3)

        fig.tight_layout()
        fig.savefig(PLOT_PATH, dpi=150)
        print(f"wrote {PLOT_PATH}")
    except ImportError:
        print("matplotlib not available, skipped plot")

    print("\nSummary:")
    for row in summary_rows:
        print(f"  k_ties={row['k_ties']} n_qubits={row['n_qubits']} "
              f"max_witness(mean/max)={row['max_witness_size_mean']}/{row['max_witness_size_max']} "
              f"dropped(mean/max)={row['dropped_candidates_mean']}/{row['dropped_candidates_max']}")


if __name__ == "__main__":
    main()
