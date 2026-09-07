"""Real transpiled gate count and depth for the partition-mixer QAOA
circuit -- same convention as `measure.py` uses for the general
basis-exchange mixer: a real compiled circuit's cost, not an analytic
estimate, transpiled to the same basis gate set for direct comparison.
"""

from __future__ import annotations

from dataclasses import dataclass

from qiskit import transpile

from partition_gadget import PartitionGadget
from partition_mixer import qaoa_circuit

TRANSPILE_BASIS = ["cx", "rz", "sx", "x"]


@dataclass
class InstanceResult:
    m: int
    k: int
    n_qubits: int
    p: int
    cx_count: int
    total_gate_count: int
    depth: int


def measure_instance(gadget: PartitionGadget, p: int) -> InstanceResult:
    qc = qaoa_circuit(gadget, p=p, measure=False)
    # Bind arbitrary concrete angles -- transpiled gate count/depth for a
    # parametrized circuit's structure doesn't depend on the angle values,
    # only real single/two-qubit gates need to be counted.
    bound = qc.assign_parameters({param: 0.37 for param in qc.parameters})
    tqc = transpile(bound, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    op_counts = tqc.count_ops()
    return InstanceResult(
        m=gadget.m,
        k=gadget.k,
        n_qubits=gadget.n_qubits,
        p=p,
        cx_count=op_counts.get("cx", 0),
        total_gate_count=sum(op_counts.values()),
        depth=tqc.depth(),
    )


def main() -> None:
    # Same "NO" instance as the classical hardness sweep (run_hardness_sweep.py)
    # -- CP-SAT needs ~6s to PROVE this one optimal, on a graph with only
    # 65 nodes / 100 edges (48 of them the qubit-bearing bipartite layer).
    gadget = PartitionGadget(items=(4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 7, 7), m=4, B=15)
    print(f"instance: m={gadget.m} k={gadget.k} n_qubits={gadget.n_qubits} "
          f"(classical CP-SAT proof time: ~6s)")
    for p in (1, 2, 3, 5):
        r = measure_instance(gadget, p)
        print(f"  p={p}: cx={r.cx_count:5d}  total_gates={r.total_gate_count:5d}  depth={r.depth:5d}")


if __name__ == "__main__":
    main()
