"""Correctness checks for `partition_mixer.py`, in the same exact-unitary
style `mixer.py`'s `verify_no_leakage`/`verify_term_no_leakage` use for
the general basis-exchange mixer -- never trust a sampled/empirical check
when an exact one is affordable.

Three checks, cheapest and most general first:

1. `verify_ising_coefficients`: pure Python, no circuit at all -- checks
   `PartitionGadget.linear_coeff`/`quad_coeff` reproduce the TRUE
   successors-squared cost (`bucket_cost`) exactly, up to a constant,
   across many random assignments.

2. `verify_weight_preservation`: the mathematical crux of why this mixer
   needs no witness search, proven exactly rather than assumed. RXX+RYY
   ("exchange_swap") preserves total Hamming weight of ANY 2 qubits it
   acts on; composing it across every pair in an m-qubit register
   (exactly what `mixer_layer` does, once per item, independently)
   therefore preserves that whole register's total weight, for ANY
   starting state -- not just the one-hot ones this gadget actually
   uses. Checked via the register's full 2^m x 2^m unitary, which is
   cheap (m is the bucket count, typically single digits) regardless of
   how many items or how large the overall gadget is: this one check,
   run once, covers every item's register at every gadget size.

3. `verify_full_gadget_no_leakage`: end-to-end exact check (full
   2^n_qubits x 2^n_qubits unitary) that starting from the gadget's
   actual one-hot initial state and applying the full mixer layer keeps
   all amplitude within the feasible (every item exactly one-hot)
   subspace. Only tractable for small `n_qubits` -- the same scale
   limitation `verify_no_leakage` documents for the general mixer.
"""

from __future__ import annotations

import itertools
import random

from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator

from partition_gadget import PartitionGadget
from partition_mixer import _swap_block, mixer_layer


def verify_ising_coefficients(gadget: PartitionGadget, max_assignments: int = 20000, seed: int = 0) -> bool:
    assignments = list(itertools.product(range(gadget.m), repeat=gadget.k))
    if len(assignments) > max_assignments:
        rng = random.Random(seed)
        assignments = rng.sample(assignments, max_assignments)
    consts = set()
    for a in assignments:
        bc = gadget.bucket_cost(a)
        hv = gadget.diag_hamiltonian_value(a)
        consts.add(round(bc - hv, 6))
    return len(consts) == 1


def verify_weight_preservation(m: int, theta_value: float = 0.37, tol: float = 1e-9) -> bool:
    if m < 2:
        return True  # nothing to mix
    qc = QuantumCircuit(m)
    block = _swap_block(theta_value).to_gate()
    for j in range(m):
        for j2 in range(j + 1, m):
            qc.append(block, [j, j2])
    U = Operator(qc).data
    for state in range(2**m):
        weight = bin(state).count("1")
        col = U[:, state]
        leaked = sum(abs(col[j]) ** 2 for j in range(len(col)) if bin(j).count("1") != weight)
        if leaked > tol:
            return False
    return True


def verify_full_gadget_no_leakage(gadget: PartitionGadget, theta_value: float = 0.37, tol: float = 1e-9) -> bool:
    n = gadget.n_qubits
    feasible_states: set[int] = set()
    for assignment in itertools.product(range(gadget.m), repeat=gadget.k):
        state = 0
        for i, j in enumerate(assignment):
            state |= 1 << gadget.qubit(i, j)
        feasible_states.add(state)

    qc = mixer_layer(gadget, theta_value)
    U = Operator(qc).data
    for state in feasible_states:
        col = U[:, state]
        leaked = sum(abs(col[j]) ** 2 for j in range(len(col)) if j not in feasible_states)
        if leaked > tol:
            return False
    return True


def main() -> None:
    gadgets = {
        "m2_tiny": PartitionGadget(items=(3,) * 6, m=2, B=9),
        "m3_symmetric": PartitionGadget(items=(4, 4, 5) * 3, m=3, B=13),
        "m4_no": PartitionGadget(items=(4, 4, 4, 4, 4, 4, 4, 6, 6, 6, 7, 7), m=4, B=15),
        "m4_mixed_yes": PartitionGadget(items=(4, 5, 6, 4, 5, 6, 4, 4, 7, 5, 5, 5), m=4, B=15),
    }
    for label, gadget in gadgets.items():
        ok = verify_ising_coefficients(gadget)
        print(f"[{label}] verify_ising_coefficients: {'OK' if ok else '*** FAILED ***'}")

    print()
    for m in range(2, 8):
        ok = verify_weight_preservation(m)
        print(f"[m={m}] verify_weight_preservation (2^{m}x2^{m} unitary): {'OK' if ok else '*** FAILED ***'}")

    print()
    for label in ("m2_tiny",):  # only size small enough for a full 2^n_qubits unitary
        gadget = gadgets[label]
        ok = verify_full_gadget_no_leakage(gadget)
        print(f"[{label}] verify_full_gadget_no_leakage (n_qubits={gadget.n_qubits}): {'OK' if ok else '*** FAILED ***'}")


if __name__ == "__main__":
    main()
