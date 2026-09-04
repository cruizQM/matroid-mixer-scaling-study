"""Tests an assumption headline result 1 (`run_scaling_study.py`) makes
and states plainly but doesn't stress-test: `k_ties` held FIXED across
the whole size sweep. That choice is what produces headline result 1's
actual punchline (cost trends flat-to-*decreasing* with size) --
`graphs.py`'s own docstring notes larger networks get progressively
SPARSER when tie count doesn't grow with them. Whether that's a
realistic model of how tie-switch count actually scales with feeder size
is a separate question from whether the graph stays sparse -- this repo
has exactly one real anchor point (IEEE 33-bus, 5 tie switches at
n_nodes=33), not a real growth curve, so any functional form is somewhat
arbitrary; a slowly-growing count (real utilities plausibly add some tie
switches as a feeder gets more complex, even if not proportionally) is
at least as defensible as perfectly flat, and lets `|E|` grow a bit
faster than headline result 1's model does.

This compares, at IDENTICAL sizes and seeds, `run_scaling_study.py`'s
own fixed `k_ties=3` against `k_ties_log(n_nodes) = round(1.43 *
ln(n_nodes))` -- calibrated so `k_ties_log(33) = 5`, matching the one
real anchor point exactly, not chosen to produce any particular result.
Both conditions reuse `measure_instance` (identical construction,
transpilation, and measurement code as headline result 1) unchanged --
only which `k_ties` gets passed in differs.

Extends the size range only as far as whole-graph brute-force
enumeration (`graphs.enumerate_spanning_trees`, headline result 1's own
tractability ceiling -- "pushed as far as brute-force spanning-tree
enumeration stays fast") remains fast under the LOG condition, which
needs a slightly smaller n_nodes ceiling than headline result 1's own
fixed-k=3 sweep did for the same reason: C(n_qubits, k_ties) grows in
k_ties too, and k_ties_log keeps climbing (3 to 5 across n_nodes=8-40)
where headline result 1's k_ties never does.

Writes results/scaling_log_ties_results.csv (per size, per seed, per
condition) and results/scaling_log_ties_summary.csv (per size, per
condition, mean/min/max).
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from measure import measure_instance

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "scaling_log_ties_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "scaling_log_ties_summary.csv"

FIXED_K_TIES = 3  # matches run_scaling_study.py's own K_TIES exactly
LOG_CALIBRATION_C = 1.4300  # calibrated so k_ties_log(33) == 5 (the real IEEE 33-bus anchor)
N_NODES_RANGE = list(range(8, 41))  # extends past headline result 1's own n_nodes=33 ceiling
SEEDS_PER_SIZE = 3


def k_ties_log(n_nodes: int) -> int:
    return max(1, round(LOG_CALIBRATION_C * math.log(n_nodes)))


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    print(f"FIXED_K_TIES={FIXED_K_TIES}, LOG condition k_ties_log(n)=round({LOG_CALIBRATION_C}*ln(n)), "
          f"n_nodes={N_NODES_RANGE[0]}-{N_NODES_RANGE[-1]}, seeds_per_size={SEEDS_PER_SIZE}\n")
    print("k_ties_log by n_nodes:", {n: k_ties_log(n) for n in N_NODES_RANGE}, "\n")

    all_rows = []  # (condition, result, seed, elapsed)
    for n_nodes in N_NODES_RANGE:
        k_log = k_ties_log(n_nodes)
        for condition, k_ties in (("fixed", FIXED_K_TIES), ("log", k_log)):
            for seed in range(SEEDS_PER_SIZE):
                t0 = time.perf_counter()
                result = measure_instance(n_nodes, k_ties, seed=seed)
                elapsed = time.perf_counter() - t0
                all_rows.append((condition, result, seed, elapsed))
            group = [r for c, r, s, e in all_rows if c == condition and r.n_nodes == n_nodes]
            cx_vals = [r.cx_count for r in group]
            print(
                f"[{condition:5s}] n_nodes={n_nodes:3d} k_ties={k_ties} n_qubits={group[0].n_qubits:3d} "
                f"cx: mean={statistics.mean(cx_vals):8.1f} min={min(cx_vals):6d} max={max(cx_vals):6d}  "
                f"connected_all={all(r.fully_connected for r in group)}",
                flush=True,
            )

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "condition", "n_nodes", "seed", "k_ties", "n_qubits", "num_spanning_trees", "term_count",
            "max_witness_size", "dropped_candidates", "fully_connected", "cx_count", "total_gate_count",
            "depth", "measurement_seconds",
        ])
        for condition, result, seed, elapsed in all_rows:
            writer.writerow([
                condition, result.n_nodes, seed, result.k_ties, result.n_qubits, result.num_spanning_trees,
                result.term_count, result.max_witness_size, result.dropped_candidates, result.fully_connected,
                result.cx_count, result.total_gate_count, result.depth, f"{elapsed:.4f}",
            ])
    print(f"\nwrote {CSV_PATH}")

    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "condition", "n_nodes", "k_ties", "n_qubits", "cx_mean", "cx_min", "cx_max",
            "depth_mean", "depth_min", "depth_max", "all_connected",
        ])
        for n_nodes in N_NODES_RANGE:
            for condition in ("fixed", "log"):
                group = [r for c, r, s, e in all_rows if c == condition and r.n_nodes == n_nodes]
                cx_vals = [r.cx_count for r in group]
                depth_vals = [r.depth for r in group]
                writer.writerow([
                    condition, n_nodes, group[0].k_ties, group[0].n_qubits,
                    f"{statistics.mean(cx_vals):.2f}", min(cx_vals), max(cx_vals),
                    f"{statistics.mean(depth_vals):.2f}", min(depth_vals), max(depth_vals),
                    all(r.fully_connected for r in group),
                ])
    print(f"wrote {SUMMARY_CSV_PATH}")


if __name__ == "__main__":
    main()
