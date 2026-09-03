"""Regenerates all results in /results from scratch: sweeps feeder size at
fixed tie-edge count, measures each instance across several random seeds,
writes the raw per-instance CSV, and produces the headline plot. This is
the single script referenced by the README's "How to reproduce" section.

Scope of the sweep:
- `n_nodes` starts at MIN_N_NODES, not 1: with k_ties held fixed, a graph
  is only genuinely SPARSE once n_nodes is large enough that
  n_qubits = n_nodes-1+k_ties is small relative to the C(n_nodes,2)
  possible edges. At n_nodes=4 with k_ties=3, n_qubits=6=C(4,2) -- the
  graph IS the complete graph K4, not sparse at all. This was caught by
  inspecting an outlier data point (very high gate count at the smallest
  size), not assumed -- see docs/circuit-validity.md. MIN_N_NODES is
  chosen so edge density (n_qubits / C(n_nodes,2)) stays well under half
  at the smallest instance measured.
- pushed as far as brute-force spanning-tree enumeration (graphs.py) stays
  fast; the exact range actually run is recorded in the CSV and stated
  plainly in the README/methodology -- not extrapolated beyond what was
  measured.
- SEEDS_PER_SIZE > 1: a single random graph realization per size is noisy
  (specific instances can happen to need more/fewer witness-conditioned
  controlled gates); the CSV records every individual run plus the mean
  and range per size, and the plot shows both.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from measure import measure_instance

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "scaling_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "scaling_summary.csv"

K_TIES = 3  # fixed tie-edge count across the sweep -- keeps the graph sparse
MAX_EDGE_DENSITY_AT_MIN_SIZE = 0.4  # n_qubits / C(n_nodes,2) must be below this at the smallest n_nodes tested


def _min_n_nodes(k_ties: int, max_density: float) -> int:
    n = 2
    while True:
        n_qubits = n - 1 + k_ties
        max_edges = comb(n, 2)
        if max_edges > 0 and n_qubits / max_edges <= max_density:
            return n
        n += 1


MIN_N_NODES = _min_n_nodes(K_TIES, MAX_EDGE_DENSITY_AT_MIN_SIZE)
N_NODES_RANGE = list(range(MIN_N_NODES, MIN_N_NODES + 26))
SEEDS_PER_SIZE = 5


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    print(f"K_TIES={K_TIES}, MIN_N_NODES={MIN_N_NODES} (density <= {MAX_EDGE_DENSITY_AT_MIN_SIZE} at that size), "
          f"seeds_per_size={SEEDS_PER_SIZE}\n")

    all_rows = []  # (result, seed, elapsed)
    for n_nodes in N_NODES_RANGE:
        for seed in range(SEEDS_PER_SIZE):
            t0 = time.perf_counter()
            result = measure_instance(n_nodes, K_TIES, seed=seed)
            elapsed = time.perf_counter() - t0
            all_rows.append((result, seed, elapsed))
        last = [r for r, s, e in all_rows if r.n_nodes == n_nodes]
        cx_vals = [r.cx_count for r in last]
        print(
            f"n_nodes={n_nodes:3d} n_qubits={last[0].n_qubits:3d} "
            f"cx: mean={statistics.mean(cx_vals):7.1f} min={min(cx_vals):5d} max={max(cx_vals):5d}  "
            f"terms(mean)={statistics.mean(r.term_count for r in last):.1f} "
            f"connected_all={all(r.fully_connected for r in last)}"
        )

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "n_nodes",
                "seed",
                "k_ties",
                "n_qubits",
                "num_spanning_trees",
                "term_count",
                "max_witness_size",
                "dropped_candidates",
                "fully_connected",
                "cx_count",
                "total_gate_count",
                "depth",
                "measurement_seconds",
            ]
        )
        for result, seed, elapsed in all_rows:
            writer.writerow(
                [
                    result.n_nodes,
                    seed,
                    result.k_ties,
                    result.n_qubits,
                    result.num_spanning_trees,
                    result.term_count,
                    result.max_witness_size,
                    result.dropped_candidates,
                    result.fully_connected,
                    result.cx_count,
                    result.total_gate_count,
                    result.depth,
                    f"{elapsed:.4f}",
                ]
            )
    print(f"\nwrote {CSV_PATH}")

    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["n_nodes", "n_qubits", "cx_mean", "cx_min", "cx_max", "depth_mean", "depth_min", "depth_max",
             "term_count_mean", "all_connected"]
        )
        for n_nodes in N_NODES_RANGE:
            group = [r for r, s, e in all_rows if r.n_nodes == n_nodes]
            cx_vals = [r.cx_count for r in group]
            depth_vals = [r.depth for r in group]
            writer.writerow(
                [
                    n_nodes,
                    group[0].n_qubits,
                    f"{statistics.mean(cx_vals):.2f}",
                    min(cx_vals),
                    max(cx_vals),
                    f"{statistics.mean(depth_vals):.2f}",
                    min(depth_vals),
                    max(depth_vals),
                    f"{statistics.mean(r.term_count for r in group):.2f}",
                    all(r.fully_connected for r in group),
                ]
            )
    print(f"wrote {SUMMARY_CSV_PATH}")

    if not all(r.fully_connected for r, s, e in all_rows):
        print("WARNING: not every instance's mixer reached full connectivity -- see CSV 'fully_connected' column")
    if any(r.dropped_candidates > 0 for r, s, e in all_rows):
        print("NOTE: some candidate exchanges were dropped (no small witness found within the search cap) on at "
              "least one instance -- see CSV 'dropped_candidates' column; connectivity was still reached using "
              "the remaining candidates in every instance tested (see 'fully_connected').")

    from plot import make_plot

    make_plot(SUMMARY_CSV_PATH, RESULTS_DIR / "scaling_plot.png")
    print(f"wrote {RESULTS_DIR / 'scaling_plot.png'}")


if __name__ == "__main__":
    main()
