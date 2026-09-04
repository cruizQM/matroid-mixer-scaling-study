"""Pushes the adaptive-alpha construction (run_cost_aware_scaling_ladder.py)
much harder on cost, explicitly trading safety for a circuit size in the
50-500 CX range that's closer to current NISQ hardware capability (see
the fidelity-estimate discussion this branch's writeup includes: even
5,000 CX already lands around 0.7% survival probability at best-case
trapped-ion error rates, let alone the 10,000-30,000+ CX the earlier
harder conditions cost).

Two levers, not one: `gain_price` alone saturates fast once every
approximate term already collapses to its cheapest option -- on the
hardest condition (long_linear, n_nodes=150), CX stayed FLAT at ~10,100
across gain_price 0.1 to 10.0, because ~75% of terms there are found via
`exact_search_max_size`'s brute-force path, which `gain_price` cannot
touch at all. Lowering `exact_search_max_size` to 0 pushes those
candidates through the cost-aware path instead, where `gain_price` can
actually act on them.

**A real failure mode found and guarded against here, not just
mentioned**: pushed far enough (exact_search_max_size=0 combined with a
large gain_price), most approximate terms independently conclude "never
fire" is their best majority-vote answer -- on one instance, 146 of 151
terms went inert, giving a circuit that looks perfect (cx=20,
unsafe_rate=0.0) for the wrong reason: it does almost nothing, not
because it's a good mixer. `fully_connected` does NOT catch this (it's
computed from which candidate pairs the exchange-graph SAMPLE says could
connect, not from which terms the compiled circuit actually fires) --
so this script tracks `active_terms` (nonzero minimized_gate_count)
explicitly and reports it alongside CX/safety at every point, instead of
trusting `fully_connected` alone.

Same conditions/sizes/seeds as run_cost_aware_scaling_ladder.py.
`exact_search_max_size=0`, `gain_price=0.05` (calibrated on the
long_linear/n_nodes=150 instance specifically -- the hardest condition --
to land in the 50-500 CX range there without collapsing to inert; other
conditions are NOT re-calibrated, so seeing how this ONE choice performs
across the whole ladder, including where it might over- or under-shoot,
is the point).

Writes results/cost_aware_scaling_ladder_aggressive_results.csv and
_summary.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer
from run_cost_aware_scaling_ladder import CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE, WALK_STEPS, N_STARTING_TREES, UNSAFE_TOL, BETA

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_aggressive_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_aggressive_summary.csv"

MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 0
GAIN_PRICE = 0.05


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)

    walked = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    truncated = build_truncated_witness_mixer(
        graph, walked, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
        adaptive=True, gain_price=GAIN_PRICE,
    )
    construction = truncated.construction
    active_terms = sum(1 for t in construction.terms if t.minimized_gate_count > 0)

    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx_count, depth = 0, 0

    rng = np.random.default_rng(seed * 1000 + 2)
    unsafe = 0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1

    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits": graph.n_edges, "n_terms": len(construction.terms), "active_terms": active_terms,
        "fully_connected": construction.fully_connected,
        "cx_count": cx_count, "depth": depth,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "worst_feasible_mass": round(float(min(masses)), 4),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    print(f"exact_search_max_size={EXACT_SEARCH_MAX_SIZE}, gain_price={GAIN_PRICE}\n")
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
            cx_vals = [r["cx_count"] for r in group]
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes, "k_ties": group[0]["k_ties"],
                "n_qubits": group[0]["n_qubits"], "n_terms": group[0]["n_terms"],
                "active_terms_mean": round(statistics.mean(r["active_terms"] for r in group), 1),
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                "cx_min": min(cx_vals), "cx_max": max(cx_vals),
                "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
                "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
                "worst_feasible_mass_min": round(min(r["worst_feasible_mass"] for r in group), 4),
                "all_connected": all(r["fully_connected"] for r in group),
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
