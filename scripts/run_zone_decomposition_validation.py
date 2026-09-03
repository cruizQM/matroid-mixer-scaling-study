"""Full validation + measurement of the zone-decomposition approach
(`scripts/zone_decomposition.py`) on the real IEEE 33-bus feeder: for each
zone plus the assembly problem, builds the actual mixer circuit,
transpiles it, verifies exact leakage-freedom (`verify_all_terms_no_leakage`
-- not just the classical witness-size/connectivity checks
`zone_decomposition.py`'s `main()` prints), and writes every subproblem's
numbers to `results/zone_decomposition_results.csv`.

This exists because witness size and classical spanning-tree connectivity
are necessary but not sufficient: they don't by themselves prove the
compiled circuit doesn't leak. Every other measurement in this repo is
backed by an exact circuit-unitary check before being trusted (see
`docs/circuit-validity.md`); this is that check for the decomposition
result.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qiskit import transpile

from graphs import enumerate_spanning_trees
from measure import TRANSPILE_BASIS
from mixer import build_matroid_mixer, mixer_circuit, verify_all_terms_no_leakage
from real_feeders import load_ieee33
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "zone_decomposition_results.csv"

TARGET_ZONES = (3, 4, 5, 6)


def measure_subproblem(label: str, graph) -> dict:
    trees = enumerate_spanning_trees(graph)
    row = {
        "subproblem": label,
        "n_nodes": graph.n_nodes,
        "n_qubits": graph.n_edges,
        "num_spanning_trees": len(trees),
    }
    if not trees:
        row.update(
            term_count=0, max_witness_size=0, dropped_candidates=0, fully_connected=False,
            no_leakage_verified="N/A (no spanning tree)", cx_count=0, total_gate_count=0, depth=0,
        )
        return row

    construction = build_matroid_mixer(graph, trees)
    max_witness = max((t.control_count for t in construction.terms), default=0)
    leak_free = verify_all_terms_no_leakage(construction) if construction.terms else True

    if construction.terms:
        qc = mixer_circuit(construction, beta=0.37)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count = op_counts.get("cx", 0)
        total_gates = sum(op_counts.values())
        depth = tqc.depth()
    else:
        cx_count = total_gates = depth = 0

    row.update(
        term_count=len(construction.terms),
        max_witness_size=max_witness,
        dropped_candidates=construction.dropped_candidates,
        fully_connected=construction.fully_connected,
        no_leakage_verified=leak_free,
        cx_count=cx_count,
        total_gate_count=total_gates,
        depth=depth,
    )
    return row


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    graph = load_ieee33()
    print(f"IEEE 33-bus: n_nodes={graph.n_nodes}  n_qubits={graph.n_edges}")

    rows = []
    for target_zones in TARGET_ZONES:
        print(f"\n=== target_zones={target_zones} ===")
        zone_of = partition_zones(graph, target_zones)
        n_zones = len(set(zone_of.values()))

        for zid in sorted(set(zone_of.values())):
            zg = build_zone_subgraph(graph, zone_of, zid)
            if zg.n_nodes < 2:
                continue
            row = measure_subproblem(f"target_zones={target_zones}/zone_{zid}", zg)
            row["target_zones"] = target_zones
            rows.append(row)
            print(f"  {row}")

        assembly_graph, boundary_edges = build_assembly_graph(graph, zone_of, n_zones)
        row = measure_subproblem(f"target_zones={target_zones}/assembly", assembly_graph)
        row["target_zones"] = target_zones
        rows.append(row)
        print(f"  {row}")

    fieldnames = [
        "target_zones", "subproblem", "n_nodes", "n_qubits", "num_spanning_trees", "term_count",
        "max_witness_size", "dropped_candidates", "fully_connected", "no_leakage_verified",
        "cx_count", "total_gate_count", "depth",
    ]
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    all_leak_free = all(r["no_leakage_verified"] in (True, "N/A (no spanning tree)") for r in rows)
    all_connected = all(r["fully_connected"] or (r["term_count"] == 0 and r["num_spanning_trees"] <= 1) for r in rows)
    max_witness_overall = max(r["max_witness_size"] for r in rows)
    print(
        f"\nSummary: all subproblems leak-free={all_leak_free}, "
        f"max witness size across every zone/assembly/target_zones config={max_witness_overall}"
    )
    if not all_leak_free:
        print("WARNING: at least one subproblem failed leakage verification -- do not cite these results until fixed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
