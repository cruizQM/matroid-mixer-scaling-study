"""QAOA circuit for `PartitionGadget`: a one-hot partition-matroid mixer
plus the pairwise cost oracle implementing its Ising Hamiltonian.

## Why this needs its own mixer, separate from mixer.py's basis-exchange one

`mixer.py`'s construction handles the GENERAL graphic matroid (any
spanning tree of any graph) by enumerating the feasible set and searching
for a witness -- other qubits whose value determines whether a given
exchange is valid -- because in general that validity genuinely depends
on the rest of the current tree (see mixer.py's own module docstring).

`partition_gadget.py`'s gadget is different: the only real freedom is a
PARTITION matroid (each item independently picks exactly one bucket out
of m). An exchange argument (removing bucket edge (u_j, v_i) from any
feasible tree disconnects exactly {v_i} union its leaves; the only edges
leaving that set are v_i's other m-1 bucket edges, since leaves have no
other neighbors and every bucket node is always root-connected) shows
swapping which bucket an item is assigned to is valid UNCONDITIONALLY,
for every feasible tree, regardless of every other item's state. So this
mixer needs no witness search and no tree enumeration at all: it is
`mixer.py`'s exact same underlying principle (basis exchange on a
matroid), but the matroid's structure here makes the general
construction's expensive step (enumerate the feasible set to find
witnesses) unnecessary -- which is fortunate, since for this gadget the
feasible set is astronomically large (m^k trees) long before it gets
interesting, so `mixer.py`'s enumeration-based construction could not
even be attempted here (see `docs/hard-instance-case-study.md`).

## Circuit

One-hot encoding: qubit `gadget.qubit(i, j)` = 1 iff item i is assigned
to bucket j. Per QAOA layer:
  - cost layer: RZ per qubit (linear_coeff) + RZZ per item-pair per
    matching bucket index (quad_coeff) -- the Ising expansion of
    `PartitionGadget`'s true successors-squared cost, verified exactly in
    `verify_partition_mixer.py`.
  - mixer layer: RXX+RYY ("exchange_swap", the same Givens-rotation
    primitive `mixer.py`'s `_swap_block` uses) between every pair of
    bucket positions within one item's register, UNCONDITIONED -- no
    controlled gates anywhere in this mixer, unlike the general one.
Initial state: every item defaults to bucket 0 (a trivially valid one-hot
assignment), same role as `mixer.py`'s "start from a fixed spanning tree"
convention.
"""

from __future__ import annotations

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import SparsePauliOp

from partition_gadget import PartitionGadget


def _swap_block(theta) -> QuantumCircuit:
    qc = QuantumCircuit(2, name="exchange_swap")
    qc.rxx(theta, 0, 1)
    qc.ryy(theta, 0, 1)
    return qc


def cost_layer(gadget: PartitionGadget, gamma) -> QuantumCircuit:
    qc = QuantumCircuit(gadget.n_qubits, name="partition_cost")
    k, m = gadget.k, gadget.m
    for i in range(k):
        coeff = gadget.linear_coeff(i)
        if coeff == 0.0:
            continue
        for j in range(m):
            qc.rz(2 * gamma * coeff, gadget.qubit(i, j))
    for i in range(k):
        for i2 in range(i + 1, k):
            coeff = gadget.quad_coeff(i, i2)
            if coeff == 0.0:
                continue
            for j in range(m):
                qc.rzz(2 * gamma * coeff, gadget.qubit(i, j), gadget.qubit(i2, j))
    return qc


def mixer_layer(gadget: PartitionGadget, beta) -> QuantumCircuit:
    qc = QuantumCircuit(gadget.n_qubits, name="partition_mixer")
    block = _swap_block(beta).to_gate()
    k, m = gadget.k, gadget.m
    for i in range(k):
        for j in range(m):
            for j2 in range(j + 1, m):
                qc.append(block, [gadget.qubit(i, j), gadget.qubit(i, j2)])
    return qc


def qaoa_circuit(gadget: PartitionGadget, p: int, measure: bool = True) -> QuantumCircuit:
    """`p` layers of (cost, mixer), with independent `Parameter`s per
    layer (`beta_l`, `gamma_l`) -- bind at execution time via
    `assign_parameters`, same as `mixer.py`'s Parameter-based circuits."""
    n = gadget.n_qubits
    qc = QuantumCircuit(n, n if measure else 0, name="partition_qaoa")
    for i in range(gadget.k):
        qc.x(gadget.qubit(i, 0))

    for layer in range(p):
        gamma = Parameter(f"gamma_{layer}")
        beta = Parameter(f"beta_{layer}")
        qc.compose(cost_layer(gadget, gamma), inplace=True)
        qc.compose(mixer_layer(gadget, beta), inplace=True)

    if measure:
        qc.measure(range(n), range(n))
    return qc


def cost_hamiltonian(gadget: PartitionGadget) -> SparsePauliOp:
    """The RZ/RZZ Hamiltonian `cost_layer` implements, as a `SparsePauliOp`
    -- lets any Aer backend (statevector, matrix_product_state, ...)
    compute its EXACT expectation value internally via
    `QuantumCircuit.save_expectation_value`, without ever materializing a
    2^n_qubits array in Python. Add `hamiltonian_constant_offset(gadget)`
    to the result to get the expected value of the true (`bucket_cost`)
    reconfiguration cost, not just the Hamiltonian's own (shifted) value."""
    n = gadget.n_qubits
    labels: list[str] = []
    coeffs: list[float] = []
    for i in range(gadget.k):
        c = gadget.linear_coeff(i)
        if c == 0.0:
            continue
        for j in range(gadget.m):
            lab = ["I"] * n
            lab[n - 1 - gadget.qubit(i, j)] = "Z"
            labels.append("".join(lab))
            coeffs.append(c)
    for i in range(gadget.k):
        for i2 in range(i + 1, gadget.k):
            c = gadget.quad_coeff(i, i2)
            if c == 0.0:
                continue
            for j in range(gadget.m):
                lab = ["I"] * n
                lab[n - 1 - gadget.qubit(i, j)] = "Z"
                lab[n - 1 - gadget.qubit(i2, j)] = "Z"
                labels.append("".join(lab))
                coeffs.append(c)
    return SparsePauliOp(labels, coeffs)


def hamiltonian_constant_offset(gadget: PartitionGadget) -> float:
    """`bucket_cost(x) - diag_hamiltonian_value(x)` is the same constant
    for every assignment x (checked exactly in
    `verify_partition_mixer.verify_ising_coefficients`) -- computed here
    from one arbitrary reference assignment."""
    ref = tuple([0] * gadget.k)
    return gadget.bucket_cost(ref) - gadget.diag_hamiltonian_value(ref)


def bits_to_assignment(gadget: PartitionGadget, bits: str) -> tuple[int, ...] | None:
    """Decodes a Qiskit measurement bitstring (index 0 = LAST character,
    Qiskit's usual little-endian classical-register convention) into a
    per-item bucket assignment, or None if some item's register is not
    exactly one-hot -- the leakage check: the mixer is claimed to
    preserve one-hot-ness exactly, this is how that claim is checked
    against real sampled/simulated output."""
    n = len(bits)
    values = [int(bits[n - 1 - q]) for q in range(n)]
    assignment = []
    for i in range(gadget.k):
        ones = [j for j in range(gadget.m) if values[gadget.qubit(i, j)] == 1]
        if len(ones) != 1:
            return None
        assignment.append(ones[0])
    return tuple(assignment)
