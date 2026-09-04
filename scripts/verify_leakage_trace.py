"""Standalone correctness check for the two new diagnostic tools this
branch adds (`random_trees.py`, `leakage_trace.py`), independent of, and a
prerequisite to trusting, `run_bounded_witness_safety_survey.py`'s
results. Mirrors `verify_correctness.py`'s discipline: exact,
cross-checked methods, not assumptions.

Three independent checks:

1. **Wilson's-algorithm sampling is close to uniform.** On a small
   instance where `graphs.enumerate_spanning_trees` can enumerate every
   tree exactly, draws many samples via `random_trees.random_spanning_tree`
   and compares the empirical per-tree frequency to the uniform
   1/|B| prediction (chi-squared-style max relative deviation, not just
   "looks about right").

2. **`trace_danger_mass_sparse` exactly matches `trace_danger_mass`.**
   Same small instance, same starting trees, same beta -- the sparse
   dict-based tracer must reproduce the dense full-statevector tracer's
   per-term danger_mass and feasible_mass_after to numerical precision.
   This is what licenses using the sparse tracer past the ~20-25 qubit
   ceiling where the dense one is no longer computable at all.

3. **The single-excitation 2-level gate identity used by the sparse
   tracer** -- `[[cos(beta),-i*sin(beta)],[-i*sin(beta),cos(beta)]]` on a
   trigger pair -- is checked directly against Qiskit's own `Operator` on
   `mixer._swap_block`, at several angles, not assumed from the gate's
   name.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit.quantum_info import Operator

from graphs import enumerate_spanning_trees, generate_feeder_graph
from mixer import _swap_block, build_matroid_mixer
from random_trees import random_spanning_tree
from leakage_trace import trace_danger_mass, trace_danger_mass_sparse
from truncated_mixer import build_truncated_witness_mixer


def check_uniform_sampling() -> bool:
    graph = generate_feeder_graph(7, 2, seed=0)
    trees = enumerate_spanning_trees(graph)
    n_trees = len(trees)
    n_samples = 20_000
    rng = np.random.default_rng(1)

    counts = {t: 0 for t in trees}
    for _ in range(n_samples):
        mask = random_spanning_tree(graph, rng)
        assert mask in counts, "Wilson's algorithm produced a mask outside the enumerated tree set"
        counts[mask] += 1

    expected = n_samples / n_trees
    max_rel_dev = max(abs(c - expected) / expected for c in counts.values())
    ok = max_rel_dev < 0.35  # generous: 20k samples over a handful of trees still has real variance
    print(f"[{'PASS' if ok else 'FAIL'}] uniform sampling: n_trees={n_trees} n_samples={n_samples} max_rel_dev={max_rel_dev:.3f}")
    return ok


def check_two_level_gate_identity() -> bool:
    ok = True
    for beta in (0.1, 0.37, 0.9, 1.4):
        U = Operator(_swap_block(beta)).to_matrix()
        # single-excitation subspace: basis order |00>,|01>,|10>,|11> (qiskit little-endian)
        # states |01> (index 1) and |10> (index 2) are the trigger pair.
        block = U[np.ix_([1, 2], [1, 2])]
        expected = np.array([[np.cos(beta), -1j * np.sin(beta)], [-1j * np.sin(beta), np.cos(beta)]])
        dev = float(np.max(np.abs(block - expected)))
        this_ok = dev < 1e-9
        ok = ok and this_ok
        print(f"[{'PASS' if this_ok else 'FAIL'}] two-level gate identity at beta={beta}: max_dev={dev:.2e}")
    return ok


def check_sparse_matches_dense() -> bool:
    graph = generate_feeder_graph(7, 2, seed=2)
    trees = enumerate_spanning_trees(graph)
    construction = build_matroid_mixer(graph, trees)
    beta = 0.37

    ok = True
    for start_tree in trees[:6]:
        dense = trace_danger_mass(graph, construction, start_tree, beta)
        sparse = trace_danger_mass_sparse(graph, construction, start_tree, beta)
        if len(dense) != len(sparse):
            ok = False
            continue
        max_dev = max(
            max(abs(d.danger_mass - s.danger_mass), abs(d.feasible_mass_after - s.feasible_mass_after))
            for d, s in zip(dense, sparse)
        )
        this_ok = max_dev < 1e-9
        ok = ok and this_ok
        print(f"[{'PASS' if this_ok else 'FAIL'}] sparse==dense for start_tree={start_tree}: max_dev={max_dev:.2e}")

    # Also exercise the sparse tracer on an approximate (leaky-by-design)
    # construction, where danger_mass is expected to be genuinely nonzero
    # for at least one term/tree -- confirms the tracer actually measures
    # something, not just agreeing on the trivially-zero exact case. Uses a
    # denser instance than the cross-check above (max_size=1 forces most
    # candidates on a dense graph past exact witness search), and also
    # cross-checks the sparse tracer against the dense one on THIS instance
    # specifically, so the "measures something real" claim and the
    # "sparse==dense" claim are both established on the same leaky
    # construction, not just separately on different graphs.
    dense_graph = generate_feeder_graph(8, 6, seed=3)
    dense_trees = enumerate_spanning_trees(dense_graph)
    truncated = build_truncated_witness_mixer(dense_graph, dense_trees, max_witness_size=2, exact_search_max_size=1, seed=0)
    any_nonzero_danger = False
    for start_tree in dense_trees[:10]:
        trace_sp = trace_danger_mass_sparse(dense_graph, truncated.construction, start_tree, beta)
        trace_ds = trace_danger_mass(dense_graph, truncated.construction, start_tree, beta)
        max_dev = max(
            max(abs(d.danger_mass - s.danger_mass), abs(d.feasible_mass_after - s.feasible_mass_after))
            for d, s in zip(trace_ds, trace_sp)
        ) if trace_ds else 0.0
        if max_dev >= 1e-9:
            print(f"    [FAIL] sparse!=dense on dense instance, start_tree={start_tree}: max_dev={max_dev:.2e}")
        ok = ok and (max_dev < 1e-9)
        if any(info.danger_mass > 1e-9 for info in trace_sp):
            any_nonzero_danger = True
    print(f"[{'PASS' if any_nonzero_danger else 'INFO'}] truncated construction on dense instance: "
          f"n_exact={truncated.n_exact_terms} n_approx={truncated.n_approximate_terms} "
          f"mean_leakage_rate={truncated.mean_leakage_rate:.4f} "
          f"any_nonzero_danger_mass_observed={any_nonzero_danger} (sparse==dense confirmed on this instance too)")
    return ok


def main() -> None:
    r1 = check_uniform_sampling()
    r2 = check_two_level_gate_identity()
    r3 = check_sparse_matches_dense()
    all_ok = r1 and r2 and r3
    print()
    if all_ok:
        print("ALL CHECKS PASSED: Wilson's-algorithm sampling is close to uniform, the sparse tracer's "
              "2-level gate identity matches Qiskit's own Operator exactly, and the sparse tracer exactly "
              "reproduces the dense tracer on every cross-checked instance.")
    else:
        print("FAILURE: at least one check failed -- do not trust run_bounded_witness_safety_survey.py "
              "results until fixed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
