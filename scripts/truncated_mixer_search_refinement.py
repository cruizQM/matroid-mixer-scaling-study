"""Two follow-ups to `truncated_witness_cap_sweep.py`'s finding (lowering
`max_witness_size` trades circuit cost for little safety loss):

## Part A: round out the low end of the cap curve

The original sweep tested `max_witness_size` in {2..6}. This adds {0, 1}
on the same density-failure instances (`n_nodes=8`, same `k_ties`/seeds),
to find where the curve actually knees, and because `max_witness_size=0`
exercises the exact code path
(`witness=(), valid_patterns=() or ((),)`) the
`exchange_term_circuit` bug (fixed earlier on this branch) was about --
worth confirming that case still behaves correctly under real use, not
just the synthetic check in `verify_leakage_trace.py`.

## Part B: does a COST-AWARE search objective beat blind capping?

`_search_truncated_witness` now takes `cost_alpha` (see `truncated_mixer.py`):
score candidates by `leak + cost_alpha * 2**size` instead of `leak` alone.
`cost_alpha=0.0` (blind capping's implicit setting) has no reason to
prefer a narrower witness even when it's nearly as good -- Part A/B of
`truncated_witness_cap_sweep.py` found this empirically (leakage doesn't
improve much as size grows). This sweeps `cost_alpha` at a FIXED
`max_witness_size=6` (giving the search room to choose ANY size up to 6,
rather than hard-capping it) on the same instances, to see whether letting
the search choose its own size per-term (informed by cost) finds better
(cheaper for the same safety, or safer for the same cost) points than
just picking one fixed cap for every term.

Writes results/truncated_mixer_cap_0_1_results.csv and
results/truncated_mixer_cost_aware_results.csv.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from graphs import enumerate_spanning_trees, generate_feeder_graph
from leakage_trace import final_feasible_mass
from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from random_trees import random_spanning_tree
from truncated_mixer import build_truncated_witness_mixer

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
BETA = 0.37

N_NODES = 8
K_TIES_RANGE = [6, 8, 10]
SEEDS = [0, 1, 2]
EXACT_SEARCH_MAX_SIZE = 2
N_STARTING_TREES = 100
UNSAFE_TOL = 1e-6


def measure_and_survey(graph, construction) -> dict:
    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx_count, depth = 0, 0

    rng = np.random.default_rng(42)
    unsafe = 0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1
    return {
        "cx_count": cx_count, "depth": depth,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
    }


def part_a() -> None:
    print("=== Part A: cap=0,1 on the density-failure family ===")
    rows = []
    for k_ties in K_TIES_RANGE:
        for seed in SEEDS:
            graph = generate_feeder_graph(N_NODES, k_ties, seed)
            trees = enumerate_spanning_trees(graph)
            for max_witness_size in [0, 1]:
                truncated = build_truncated_witness_mixer(
                    graph, trees, max_witness_size=max_witness_size,
                    exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
                )
                stats = measure_and_survey(graph, truncated.construction)
                row = {
                    "k_ties": k_ties, "seed": seed, "max_witness_size": max_witness_size,
                    "n_qubits": graph.n_edges, "n_terms": len(truncated.construction.terms),
                    "n_approximate_terms": truncated.n_approximate_terms,
                    "mean_term_leakage_rate": round(truncated.mean_leakage_rate, 4),
                    "fully_connected": truncated.construction.fully_connected,
                    **stats,
                }
                rows.append(row)
                print(f"  {row}")

    with open(RESULTS_DIR / "truncated_mixer_cap_0_1_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {RESULTS_DIR / 'truncated_mixer_cap_0_1_results.csv'}")


def part_b() -> None:
    print("\n=== Part B: cost-aware search (fixed max_witness_size=6, varying cost_alpha) ===")
    cost_alphas = [0.0, 0.0005, 0.001, 0.005, 0.01, 0.05]
    rows = []
    for k_ties in K_TIES_RANGE:
        for seed in SEEDS:
            graph = generate_feeder_graph(N_NODES, k_ties, seed)
            trees = enumerate_spanning_trees(graph)
            for cost_alpha in cost_alphas:
                truncated = build_truncated_witness_mixer(
                    graph, trees, max_witness_size=6,
                    exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed, cost_alpha=cost_alpha,
                )
                stats = measure_and_survey(graph, truncated.construction)
                sizes = [t.control_count for t in truncated.construction.terms]
                row = {
                    "k_ties": k_ties, "seed": seed, "cost_alpha": cost_alpha,
                    "n_qubits": graph.n_edges, "n_terms": len(truncated.construction.terms),
                    "n_approximate_terms": truncated.n_approximate_terms,
                    "mean_witness_size": round(sum(sizes) / len(sizes), 2) if sizes else 0.0,
                    "max_witness_size_used": max(sizes, default=0),
                    "mean_term_leakage_rate": round(truncated.mean_leakage_rate, 4),
                    "fully_connected": truncated.construction.fully_connected,
                    **stats,
                }
                rows.append(row)
                print(f"  {row}")

    with open(RESULTS_DIR / "truncated_mixer_cost_aware_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {RESULTS_DIR / 'truncated_mixer_cost_aware_results.csv'}")


if __name__ == "__main__":
    RESULTS_DIR.mkdir(exist_ok=True)
    part_a()
    part_b()
