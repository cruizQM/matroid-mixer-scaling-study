"""Standalone correctness check: exact circuit-unitary verification that
the mixer never leaks probability outside the feasible subspace (Hadfield's
validity condition), on a handful of small representative instances --
independent of, and a prerequisite to trusting, the scaling numbers in
/results. Run this to reproduce the correctness claim itself, not just the
performance numbers.

Runs BOTH `verify_no_leakage` (checks the full composed circuit's
2^n_qubits x 2^n_qubits unitary directly -- the strongest possible check,
but grows expensive fast: ~2s already at n_qubits=10 in testing) and
`verify_all_terms_no_leakage` (checks each term exactly but only on its
own few acting qubits, exploiting that a term is provably identity
elsewhere -- stays under 20ms regardless of instance size in testing) on
every instance here, and asserts they agree. This is what establishes that
the cheap per-term check is a trustworthy substitute for real-feeder-scale
validation (scripts/run_real_feeder_validation.py), where the full check
is not just slow but physically infeasible (2^37 for the IEEE 33-bus
case). An earlier version of `verify_all_terms_no_leakage` had its own bug
(assumed a full deterministic swap rather than the partial rotation a
generic-angle exp(-i*beta*H) actually performs) -- caught by exactly this
cross-check disagreeing, before the per-term check was trusted anywhere.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from graphs import enumerate_spanning_trees, generate_feeder_graph
from mixer import build_matroid_mixer, verify_all_terms_no_leakage, verify_no_leakage

# Small enough for exact 2^n_qubits statevector simulation (the full check;
# the per-term check has no such limit, but is only trustworthy where cross-
# validated against the full check, which is why these stay small).
INSTANCES = [(6, 2, seed) for seed in range(3)] + [(8, 2, seed) for seed in range(3)] + [(8, 3, seed) for seed in range(3)]


def main() -> None:
    all_ok = True
    for n_nodes, k_ties, seed in INSTANCES:
        graph = generate_feeder_graph(n_nodes, k_ties, seed)
        trees = enumerate_spanning_trees(graph)
        construction = build_matroid_mixer(graph, trees)
        leak_free_full = verify_no_leakage(construction)
        leak_free_per_term = verify_all_terms_no_leakage(construction)
        max_witness = max((t.control_count for t in construction.terms), default=0)
        agree = leak_free_full == leak_free_per_term
        status = "PASS" if (leak_free_full and leak_free_per_term and construction.fully_connected and agree) else "FAIL"
        print(
            f"[{status}] n_nodes={n_nodes} k_ties={k_ties} seed={seed} n_qubits={graph.n_edges} "
            f"|B|={len(trees)} terms={len(construction.terms)} max_witness={max_witness} "
            f"connected={construction.fully_connected} "
            f"no_leakage(full)={leak_free_full} no_leakage(per_term)={leak_free_per_term} agree={agree}"
        )
        all_ok = all_ok and leak_free_full and leak_free_per_term and construction.fully_connected and agree

    print()
    if all_ok:
        print("ALL INSTANCES: mixer is exactly leakage-free and fully connected; "
              "the two independent verification methods agree everywhere.")
    else:
        print("FAILURE: at least one instance failed verification, or the two methods disagree "
              "-- do not trust /results or real_feeder_results.csv until fixed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
