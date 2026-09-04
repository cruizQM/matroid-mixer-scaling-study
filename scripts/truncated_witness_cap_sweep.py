"""Option 1 from the circuit-cost investigation (measure_truncated_mixer.py):
does capping `max_witness_size` LOWER trade circuit cost for leakage in a
way that's actually worth it, on the same density-failure instances where
the truncated mixer was found to cost 17K-32K CX at just 15-17 qubits?

No new algorithm -- `build_truncated_witness_mixer`'s search already only
walks up to `max_witness_size`; this just varies that one parameter and
measures both sides of the trade directly: transpiled circuit cost (as
`measure_truncated_mixer.py` did) AND actual safety (via
`leakage_trace.final_feasible_mass` on a sample of real starting trees,
not just the term-level majority-vote leakage rate `truncated_mixer.py`
already reports -- the same distinction `docs/bounded-witness-mixer.md`'s
"Does the declared leakage rate actually cost feasible mass?" section
makes, applied here across the witness-cap axis instead of just network
size).

Fixed on the k_ties in {6, 8, 10} density-failure instances
(`n_nodes=8`, from `tie_density_sweep.py`) where the exact construction
is already known to fail or become incomplete -- this is exactly the
regime the cost finding matters in.

Writes results/truncated_witness_cap_sweep_results.csv.
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
MAX_WITNESS_SIZE_RANGE = [2, 3, 4, 5, 6]
EXACT_SEARCH_MAX_SIZE = 2
N_STARTING_TREES = 100
UNSAFE_TOL = 1e-6


def run_one(k_ties: int, seed: int, max_witness_size: int) -> dict:
    graph = generate_feeder_graph(N_NODES, k_ties, seed)
    trees = enumerate_spanning_trees(graph)

    truncated = build_truncated_witness_mixer(
        graph, trees, max_witness_size=max_witness_size,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
    )
    construction = truncated.construction

    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx_count, depth = 0, 0

    rng = np.random.default_rng(seed * 1000 + max_witness_size)
    unsafe = 0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1

    return {
        "k_ties": k_ties, "seed": seed, "max_witness_size": max_witness_size,
        "n_qubits": graph.n_edges, "n_spanning_trees": len(trees),
        "n_terms": len(construction.terms),
        "n_approximate_terms": truncated.n_approximate_terms,
        "mean_term_leakage_rate": round(truncated.mean_leakage_rate, 4),
        "fully_connected": construction.fully_connected,
        "cx_count": cx_count, "depth": depth,
        "n_starting_trees_surveyed": N_STARTING_TREES,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "worst_feasible_mass": round(float(min(masses)), 4),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for k_ties in K_TIES_RANGE:
        for seed in SEEDS:
            for max_witness_size in MAX_WITNESS_SIZE_RANGE:
                row = run_one(k_ties, seed, max_witness_size)
                rows.append(row)
                print(f"  {row}")

    with open(RESULTS_DIR / "truncated_witness_cap_sweep_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {RESULTS_DIR / 'truncated_witness_cap_sweep_results.csv'}")

    print("\nMean over seeds, per (k_ties, max_witness_size):")
    for k_ties in K_TIES_RANGE:
        for max_witness_size in MAX_WITNESS_SIZE_RANGE:
            group = [r for r in rows if r["k_ties"] == k_ties and r["max_witness_size"] == max_witness_size]
            mean_cx = sum(r["cx_count"] for r in group) / len(group)
            mean_unsafe = sum(r["unsafe_rate"] for r in group) / len(group)
            mean_mass = sum(r["mean_feasible_mass"] for r in group) / len(group)
            all_connected = all(r["fully_connected"] for r in group)
            print(f"  k_ties={k_ties} cap={max_witness_size}: mean_cx={mean_cx:.0f} "
                  f"mean_unsafe_rate={mean_unsafe:.3f} mean_feasible_mass={mean_mass:.4f} "
                  f"all_connected={all_connected}")


if __name__ == "__main__":
    main()
