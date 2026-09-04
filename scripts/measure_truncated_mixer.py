"""Answers a question `docs/bounded-witness-mixer.md` raised but didn't
measure: does the bounded-witness (truncated) mixer's actual CIRCUIT COST
(transpiled CX count, depth) -- the metric headline results 1 and 2 are
built around -- scale any better than the exact construction's? Staying
CONNECTED where the exact construction fails (verified in
`run_bounded_witness_safety_survey.py`, true by construction) is a
different claim from "the resulting circuit is cheap," and the two were
never checked together until this script.

Two separate comparisons, because "scales better" was ambiguous between
them:

## A. Does the truncated mixer's OWN circuit cost stay bounded with size?

Same methodology as `measure.py`/`run_scaling_study.py` (transpile at
`optimization_level=1`, `TRANSPILE_BASIS`), applied for the first time to
`truncated_mixer.py`'s output, on the same long-range-tie family and size
range `run_bounded_witness_safety_survey.py` already used (so the two
scripts' results sit next to each other directly).

## B. Head-to-head against the EXACT construction, on the DENSITY failure
family specifically

`tie_density_sweep.py` found tie density alone (not range) collapses the
exact whole-graph construction's connectivity -- and, unlike the
range-driven failure, zone decomposition was never shown to fix a
density-driven failure (explicitly flagged as an open question in
`docs/bounded-witness-mixer.md`'s "Honest scope" section). This is
therefore the case where the truncated mixer's differentiation from
decomposition actually matters, and where a head-to-head is fair: same
small `n_nodes=8` instances (full `enumerate_spanning_trees`, not a
sample, so both constructions see the identical tree set), same seeds,
reporting the exact construction's cost ONLY when it's actually
`fully_connected` (a disconnected exact mixer produces an incomplete,
not-comparable circuit) alongside the truncated construction's cost and
connectivity (always connected, by construction) at every point,
including where the exact construction has already failed.

Writes results/truncated_mixer_scaling_results.csv (part A) and
results/truncated_vs_exact_density_results.csv (part B).
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qiskit import transpile

from graphs import enumerate_spanning_trees, generate_feeder_graph, generate_feeder_graph_long_range_ties
from measure import TRANSPILE_BASIS
from mixer import build_matroid_mixer, mixer_circuit
from random_trees import random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
BETA = 0.37

# Part A: same family/sizes/build parameters as run_bounded_witness_safety_survey.py
N_NODES_RANGE_A = [10, 30, 60, 100, 150]
SEEDS_PER_SIZE_A = 3
K_TIES_A = 5
WALK_STEPS = 400
MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 2

# Part B: same family/sizes/seeds as tie_density_sweep.py
N_NODES_B = 8
K_TIES_RANGE_B = [2, 4, 6, 8, 10]
SEEDS_PER_K_B = 4


def measure_circuit(construction) -> dict:
    qc = mixer_circuit(construction, beta=BETA)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    op_counts = tqc.count_ops()
    return {"cx_count": op_counts.get("cx", 0), "total_gate_count": sum(op_counts.values()), "depth": tqc.depth()}


def run_part_a(n_nodes: int, seed: int) -> dict:
    graph = generate_feeder_graph_long_range_ties(n_nodes, K_TIES_A, seed)
    walked = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    truncated = build_truncated_witness_mixer(
        graph, walked, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
    )
    t0 = time.perf_counter()
    cost = measure_circuit(truncated.construction)
    elapsed = time.perf_counter() - t0
    return {
        "n_nodes": n_nodes, "seed": seed, "n_qubits": graph.n_edges,
        "n_terms": len(truncated.construction.terms),
        "n_exact_terms": truncated.n_exact_terms, "n_approximate_terms": truncated.n_approximate_terms,
        "fully_connected": truncated.construction.fully_connected,
        **cost, "transpile_elapsed_s": round(elapsed, 2),
    }


def run_part_b(k_ties: int, seed: int) -> dict:
    graph = generate_feeder_graph(N_NODES_B, k_ties, seed)
    trees = enumerate_spanning_trees(graph)

    exact = build_matroid_mixer(graph, trees)
    exact_cost = measure_circuit(exact) if exact.terms else {"cx_count": 0, "total_gate_count": 0, "depth": 0}

    truncated = build_truncated_witness_mixer(
        graph, trees, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
    )
    truncated_cost = measure_circuit(truncated.construction)

    return {
        "k_ties": k_ties, "seed": seed, "n_qubits": graph.n_edges, "n_spanning_trees": len(trees),
        "exact_n_terms": len(exact.terms), "exact_dropped": exact.dropped_candidates,
        "exact_fully_connected": exact.fully_connected,
        "exact_cx_count": exact_cost["cx_count"], "exact_depth": exact_cost["depth"],
        "truncated_n_terms": len(truncated.construction.terms),
        "truncated_n_approximate_terms": truncated.n_approximate_terms,
        "truncated_fully_connected": truncated.construction.fully_connected,
        "truncated_cx_count": truncated_cost["cx_count"], "truncated_depth": truncated_cost["depth"],
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    print("=== Part A: truncated mixer's own circuit-cost trend vs. network size ===")
    rows_a = []
    for n_nodes in N_NODES_RANGE_A:
        for seed in range(SEEDS_PER_SIZE_A):
            row = run_part_a(n_nodes, seed)
            rows_a.append(row)
            print(f"  {row}")
    with open(RESULTS_DIR / "truncated_mixer_scaling_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_a[0].keys()))
        writer.writeheader()
        writer.writerows(rows_a)
    print(f"wrote {RESULTS_DIR / 'truncated_mixer_scaling_results.csv'}")

    print("\n=== Part B: truncated vs. exact, head-to-head on the density-failure family ===")
    rows_b = []
    for k_ties in K_TIES_RANGE_B:
        for seed in range(SEEDS_PER_K_B):
            row = run_part_b(k_ties, seed)
            rows_b.append(row)
            print(f"  {row}")
    with open(RESULTS_DIR / "truncated_vs_exact_density_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_b[0].keys()))
        writer.writeheader()
        writer.writerows(rows_b)
    print(f"wrote {RESULTS_DIR / 'truncated_vs_exact_density_results.csv'}")


if __name__ == "__main__":
    main()
