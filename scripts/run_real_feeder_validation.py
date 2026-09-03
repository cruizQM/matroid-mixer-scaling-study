"""Validates the matroid basis-exchange mixer construction against real
published distribution-feeder topology (IEEE 33-bus, Baran & Wu 1989 --
see real_feeders.py for how it's loaded and why CATS is not used here),
rather than only the synthetic graphs in run_scaling_study.py. Runs the
same correctness check (exact circuit-unitary, zero leakage) and the same
measurement (real transpiled CX count/depth) as the rest of this repo.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qiskit import transpile

from measure import TRANSPILE_BASIS
from mixer import build_matroid_mixer, mixer_circuit, verify_all_terms_no_leakage
from real_feeders import load_ieee33

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "real_feeder_results.csv"


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    graph = load_ieee33()
    print(f"IEEE 33-bus (Baran & Wu 1989), via pandapower case33bw:")
    print(f"  n_nodes={graph.n_nodes}  k_ties={graph.k_ties}  n_qubits={graph.n_edges}")

    from graphs import enumerate_spanning_trees

    t0 = time.perf_counter()
    trees = enumerate_spanning_trees(graph)
    t_enum = time.perf_counter() - t0
    print(f"  spanning trees |B|={len(trees)}  (enumeration: {t_enum:.2f}s)")

    t0 = time.perf_counter()
    construction = build_matroid_mixer(graph, trees, verbose=True)
    t_construct = time.perf_counter() - t0
    max_witness = max((t.control_count for t in construction.terms), default=0)
    print(
        f"  terms={len(construction.terms)}  max_witness={max_witness}  "
        f"dropped={construction.dropped_candidates}  connected={construction.fully_connected}  "
        f"(construction: {t_construct:.2f}s)"
    )

    qc = mixer_circuit(construction, beta=0.37)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    op_counts = tqc.count_ops()
    cx_count = op_counts.get("cx", 0)
    total_gates = sum(op_counts.values())
    depth = tqc.depth()
    print(f"  transpiled: cx={cx_count}  total_gates={total_gates}  depth={depth}")

    print("\n  Verifying exact correctness per term (the full 2^37 circuit unitary is not just slow but")
    print("  physically infeasible to hold in memory -- see mixer.py's verify_term_no_leakage docstring")
    print("  and docs/circuit-validity.md for why per-term verification is an exact, not approximate,")
    print("  substitute, cross-validated against the full check on smaller instances)...")
    leak_free = verify_all_terms_no_leakage(construction)
    print(f"  no_leakage(per_term)={leak_free}")

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "source",
                "n_nodes",
                "k_ties",
                "n_qubits",
                "num_spanning_trees",
                "term_count",
                "max_witness_size",
                "dropped_candidates",
                "fully_connected",
                "no_leakage_verified",
                "cx_count",
                "total_gate_count",
                "depth",
            ]
        )
        writer.writerow(
            [
                "IEEE33_baran_wu_1989_pandapower_case33bw",
                graph.n_nodes,
                graph.k_ties,
                graph.n_edges,
                len(trees),
                len(construction.terms),
                max_witness,
                construction.dropped_candidates,
                construction.fully_connected,
                leak_free,
                cx_count,
                total_gates,
                depth,
            ]
        )
    print(f"\nwrote {CSV_PATH}")

    if not (construction.fully_connected and leak_free):
        print("WARNING: the real-feeder mixer failed a correctness check -- do not cite this result until fixed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
