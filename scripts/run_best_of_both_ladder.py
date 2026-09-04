"""Combines the two constructions that turned out to trade wins and
losses with each other (run_cost_aware_scaling_ladder.py's whole-graph
adaptive-alpha mixer vs. run_decomposed_cost_aware_ladder.py's decomposed
version): build BOTH at every (condition, n_nodes, seed), transpile
both, and keep whichever costs fewer CX gates. Neither construction
dominated the other -- decomposition won by up to 7.2x on hard,
large-scale (long-range, real-scale) instances but LOST by up to 2.2x on
easy or small instances, where its own assembly-graph overhead isn't
earned back. This is a purely mechanical "try both, keep the cheaper
one" selector -- not a smarter a-priori rule for predicting which will
win (that's a natural follow-up once this establishes there IS a
consistent pattern to predict), so it costs roughly 2x the compute of
either method alone, which is fine for CONSTRUCTION time (seconds) even
if it wouldn't be for something expensive to redo repeatedly.

Reports, per (condition, n_nodes): which method won each seed, and the
resulting best-of-both CX/safety numbers, so the actual achievable
numbers (not just "decomposition sometimes wins") are in hand.

Writes results/best_of_both_ladder_results.csv (per condition, size,
seed) and _summary.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_cost_aware_scaling_ladder import run_one as run_whole_graph
from run_decomposed_cost_aware_ladder import run_one as run_decomposed
from run_cost_aware_scaling_ladder import CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "best_of_both_ladder_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "best_of_both_ladder_summary.csv"


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    whole = run_whole_graph(condition, n_nodes, seed)
    decomp = run_decomposed(condition, n_nodes, seed)

    whole_cx = whole["cx_count"]
    decomp_cx = decomp["total_cx"]
    method = "whole_graph" if whole_cx <= decomp_cx else "decomposed"
    chosen = whole if method == "whole_graph" else decomp

    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed,
        "method": method,
        "whole_graph_cx": whole_cx, "decomposed_cx": decomp_cx,
        "chosen_cx": chosen["cx_count"] if method == "whole_graph" else chosen["total_cx"],
        "chosen_unsafe_rate": chosen["unsafe_rate"],
        "chosen_mean_feasible_mass": chosen["mean_feasible_mass"],
        "chosen_worst_feasible_mass": chosen["worst_feasible_mass"],
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            for seed in range(SEEDS_PER_SIZE):
                row = run_one(condition, n_nodes, seed)
                rows.append(row)
                print(f"  {row}", flush=True)

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    summary_rows = []
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            group = [r for r in rows if r["condition"] == condition and r["n_nodes"] == n_nodes]
            cx_vals = [r["chosen_cx"] for r in group]
            n_decomposed_wins = sum(1 for r in group if r["method"] == "decomposed")
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes,
                "n_seeds_decomposed_won": n_decomposed_wins, "n_seeds_total": len(group),
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                "cx_min": min(cx_vals), "cx_max": max(cx_vals),
                "unsafe_rate_mean": round(statistics.mean(r["chosen_unsafe_rate"] for r in group), 4),
                "mean_feasible_mass_mean": round(statistics.mean(r["chosen_mean_feasible_mass"] for r in group), 4),
                "worst_feasible_mass_min": round(min(r["chosen_worst_feasible_mass"] for r in group), 4),
            })
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")

    print("\nSummary:")
    for row in summary_rows:
        print(f"  {row}")


if __name__ == "__main__":
    main()
