"""Measurement of one instance: real transpiled gate count and depth, not
an analytic estimate -- consistent with this study's goal of reporting
what a real compiled circuit costs, not a theoretical formula.
"""

from __future__ import annotations

from dataclasses import dataclass

from qiskit import transpile

from graphs import enumerate_spanning_trees, generate_feeder_graph
from mixer import build_matroid_mixer, mixer_circuit

# A basis representative of gate-model hardware (single-qubit rotations +
# CX), used only to get a concrete, reproducible gate count/depth number --
# not a claim about any specific device.
TRANSPILE_BASIS = ["cx", "rz", "sx", "x"]


@dataclass
class InstanceResult:
    n_nodes: int
    k_ties: int
    n_qubits: int
    num_spanning_trees: int
    term_count: int
    max_witness_size: int
    dropped_candidates: int
    fully_connected: bool
    cx_count: int
    total_gate_count: int
    depth: int


def measure_instance(n_nodes: int, k_ties: int, seed: int, beta_value: float = 0.37) -> InstanceResult:
    graph = generate_feeder_graph(n_nodes, k_ties, seed)
    trees = enumerate_spanning_trees(graph)
    construction = build_matroid_mixer(graph, trees)

    qc = mixer_circuit(construction, beta=beta_value)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    op_counts = tqc.count_ops()

    return InstanceResult(
        n_nodes=n_nodes,
        k_ties=k_ties,
        n_qubits=construction.n_qubits,
        num_spanning_trees=len(trees),
        term_count=len(construction.terms),
        max_witness_size=max((t.control_count for t in construction.terms), default=0),
        dropped_candidates=construction.dropped_candidates,
        fully_connected=construction.fully_connected,
        cx_count=op_counts.get("cx", 0),
        total_gate_count=sum(op_counts.values()),
        depth=tqc.depth(),
    )
