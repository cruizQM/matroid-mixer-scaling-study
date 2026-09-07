"""E.ON's secondary requirement: show a classical tensor-network (MPS)
simulator also struggles on the hard-instance case study
(docs/hard-instance-case-study.md), not just an exact MILP solver. This
also closes (as far as is honestly possible) the approximation-ratio gap
that document flags as open: the expected cost at fixed, representative
angles (0.37 half-turns, the same convention `measure.py`/
`verify_no_leakage` use elsewhere in this repo for "a concrete,
reproducible value" rather than a claimed optimum) IS the approximation
ratio metric used throughout -- no COBYLA optimization loop is
introduced, consistent with this repo's existing scope (circuit
characterization, not running an actual optimizer).

## Two things measured together, because they turn out to be entangled (literally)

1. **Exact ground truth where possible.** `PartitionGadget.n_qubits =
   3*m^2` (not linear in m!) -- exact statevector simulation is only
   tractable up to m=3 (27 qubits; m=4 is 48 qubits, confirmed to need
   4.3 exabytes and fail outright). At m=2 and m=3, the EXACT expected
   cost (via `Statevector`-equivalent `save_expectation_value`, no
   sampling noise at all) is computed directly, at p=1..5, and IS the
   real approximation-ratio number for those sizes.

2. **MPS bond-dimension convergence.** A calibration run at m=3 (where
   exact ground truth exists) found something specific, not assumed: at
   p=1, MPS is EXACT even at bond dimension 2 -- because the circuit
   starts from a computational basis state, and a diagonal (RZ/RZZ) cost
   layer applied to a basis state creates no entanglement at all,
   regardless of how many qubits its terms connect; the mixer layer only
   acts WITHIN each item's own small register, never across items. Real
   entanglement only appears once a SECOND cost layer acts on the
   superposition the first mixer created -- so difficulty scales with
   QAOA depth p, not just instance size m. At p=2 and p=3, low bond
   dimension (2-8) diverges substantially from the exact value; bond
   dimension ~32-64 was needed to converge to within ~2 at m=3, p=3.

   This script fixes p (deliberately >1, to actually engage the
   phenomenon above) and sweeps m past the exact-feasible boundary,
   at several bond dimensions, to see whether the bond dimension that
   sufficed at m=3 keeps sufficing as m grows into the case study's
   actual hard-instance range -- if it does, MPS remains a viable
   classical shortcut; if the bond dimension needed to stay converged
   keeps climbing, that is direct evidence a classical tensor-network
   method also struggles here.
"""

from __future__ import annotations

import csv
import signal
import time
from dataclasses import dataclass
from pathlib import Path

from qiskit_aer import AerSimulator

from partition_gadget import PartitionGadget
from partition_mixer import cost_hamiltonian, hamiltonian_constant_offset, qaoa_circuit

REPRESENTATIVE_ANGLE = 0.37  # same convention as measure.py / verify_no_leakage elsewhere in this repo
P_LAYERS = 2  # p=1 is exactly solvable at any bond dimension (see module docstring) -- not informative here
BOND_DIMS = (4, 16, 64)
PER_RUN_TIMEOUT_S = 240.0  # enforced via SIGALRM -- a run hitting this IS the data point (MPS became impractical)


class _TimeoutError(Exception):
    pass


def _run_with_timeout(fn, timeout_s: float):
    def _handler(signum, frame):
        raise _TimeoutError()

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(int(timeout_s))
    try:
        return fn()
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)


def instance_for_m(m: int) -> PartitionGadget:
    """Same YES-instance family used throughout the case study (fours=2m,
    sevens=m -- always has a perfect partition), for a consistent series
    across m."""
    items = tuple([4] * (2 * m) + [7] * m)
    return PartitionGadget(items=items, m=m, B=15)


def bound_circuit_with_hamiltonian(gadget: PartitionGadget, p: int):
    H = cost_hamiltonian(gadget)
    offset = hamiltonian_constant_offset(gadget)
    qc = qaoa_circuit(gadget, p=p, measure=False)
    qc = qc.assign_parameters({param: REPRESENTATIVE_ANGLE for param in qc.parameters})
    qc.save_expectation_value(H, range(gadget.n_qubits), label="H")
    return qc, offset


@dataclass
class RunResult:
    m: int
    n_qubits: int
    method: str
    bond_dim: int | None
    expected_cost: float | None
    wall_clock_s: float
    exact_optimum: int
    error: str


# NOT using `transpile()` here at all, deliberately -- two independent problems ruled it out, not one:
# (1) `transpile(qc, sim)` derives a coupling map from the backend object, and AerSimulator's inferred coupling map
#     caps out at 63 qubits regardless of method (confirmed via `sim.configuration().n_qubits == 63`) -- including
#     matrix_product_state, which has no such real limit. Crashed this script outright at m=5 (75 qubits) with
#     `CircuitTooWideForTarget`, uncaught because it happened outside the try/except (a real bug fixed below too:
#     the first run's crash lost all of m=2..4's already-good data because the CSV was only written once, at the
#     very end of main() -- now written incrementally per row).
# (2) The natural-looking fix, `transpile(qc, basis_gates=[...], optimization_level=1)` with NO backend, fails
#     differently and more confusingly: this circuit contains `save_expectation_value`, a special Aer instruction
#     (not a standard gate), and without a real backend/target object to recognize it, transpile's default target
#     inference degrades into an unrelated generic target and can no longer even find equivalence rules for
#     `rzz`/`rxx`/`ryy` ("unable to translate ... to target basis: {if_else, ...}").
# The fix: skip transpile entirely and use `.decompose(reps=4)`, which fully unravels the circuit's own
# `exchange_swap` custom gate and Qiskit's standard rxx/ryy/rzz gates down to {u, cx, save_expval} -- all of which
# every AerSimulator method (including matrix_product_state) executes natively, with no coupling-map or
# target-inference step involved at all.
def _prepare(qc):
    return qc.decompose(reps=4)


def run_exact(gadget: PartitionGadget, p: int) -> RunResult:
    qc, offset = bound_circuit_with_hamiltonian(gadget, p)
    sim = AerSimulator(method="statevector")
    t0 = time.perf_counter()
    try:
        tqc = _prepare(qc)
        result = _run_with_timeout(lambda: sim.run(tqc, shots=1).result(), PER_RUN_TIMEOUT_S)
        cost = result.data(0)["H"].real + offset
        error = ""
    except _TimeoutError:
        cost, error = None, f"TIMEOUT after {PER_RUN_TIMEOUT_S:.0f}s"
    except Exception as e:  # noqa: BLE001 -- reporting failure IS the data point here
        cost, error = None, str(e)[:200]
    elapsed = time.perf_counter() - t0
    return RunResult(gadget.m, gadget.n_qubits, "statevector", None, cost, elapsed, gadget.exact_optimum(), error)


def run_mps(gadget: PartitionGadget, p: int, bond_dim: int) -> RunResult:
    qc, offset = bound_circuit_with_hamiltonian(gadget, p)
    sim = AerSimulator(method="matrix_product_state", matrix_product_state_max_bond_dimension=bond_dim)
    t0 = time.perf_counter()
    try:
        tqc = _prepare(qc)
        result = _run_with_timeout(lambda: sim.run(tqc, shots=1).result(), PER_RUN_TIMEOUT_S)
        cost = result.data(0)["H"].real + offset
        error = ""
    except _TimeoutError:
        cost, error = None, f"TIMEOUT after {PER_RUN_TIMEOUT_S:.0f}s"
    except Exception as e:  # noqa: BLE001
        cost, error = None, str(e)[:200]
    elapsed = time.perf_counter() - t0
    return RunResult(gadget.m, gadget.n_qubits, "mps", bond_dim, cost, elapsed, gadget.exact_optimum(), error)


def main() -> None:
    out_path = Path(__file__).resolve().parent.parent / "results" / "mps_scaling_check.csv"
    fieldnames = ["m", "n_qubits", "method", "bond_dim", "expected_cost", "wall_clock_s", "exact_optimum", "error"]

    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        def emit(r: RunResult) -> None:
            writer.writerow(vars(r))
            f.flush()  # a crash mid-sweep (already happened once -- see TRANSPILE_BASIS comment above) should not
            # lose already-computed rows

        for m in (2, 3, 4, 5, 6, 7):
            gadget = instance_for_m(m)
            print(f"=== m={m} n_qubits={gadget.n_qubits} p={P_LAYERS} exact_optimum={gadget.exact_optimum()} ===",
                  flush=True)

            if gadget.n_qubits <= 27:  # m<=3 -- see module docstring for why m=4 (48 qubits) already fails
                r = run_exact(gadget, P_LAYERS)
                emit(r)
                print(f"  statevector (exact): cost={r.expected_cost}  time={r.wall_clock_s:.1f}s  {r.error}",
                      flush=True)

            for bd in BOND_DIMS:
                r = run_mps(gadget, P_LAYERS, bd)
                emit(r)
                print(f"  mps bond_dim={bd:4d}: cost={r.expected_cost}  time={r.wall_clock_s:7.2f}s  {r.error}",
                      flush=True)
                if r.wall_clock_s > PER_RUN_TIMEOUT_S:
                    print(f"  (bond_dim={bd} already took {r.wall_clock_s:.0f}s -- stopping bond-dim increases "
                          f"for this m)", flush=True)
                    break

    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
