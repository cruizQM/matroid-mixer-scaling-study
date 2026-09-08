"""Does the hard-instance family become hard for EVERY classical route
once the bucket size B grows -- not just for a general-purpose solver
on the natural formulation?

docs/hard-instance-case-study.md's Part 1 sweep has a caveat its own
cross-check exposes: at fixed B=15 with items from {4,5,6,7}, the exact
DP (`PartitionGadget.exact_optimum`) solves every instance in
milliseconds, because its memo key -- the sorted tuple of partial
bucket sums -- collides massively when items take only four values.
Strong NP-hardness of 3-PARTITION needs B to grow with m. With random
large-B items the sums don't collide and the DP hits its true bound,
about m^k/m! states: measured directly (scratch probe before this
script), instant at m=4, >300s at m=5, B=1000. So genuine hardness
needs large B *and* m >= 5.

## Why the weighted formulation, and why B=15 is re-run under it

The Part 1 sweep encodes item weight a_i as (a_i - 1) leaf nodes hung
off each item node (Khodabakhsh's unit-demand reduction, P2). At
B=10,000 that is ~60,000 leaves, every one of them forced -- CP-SAT's
time would then measure presolve on a huge model, not the combinatorial
core. Khodabakhsh's original objective (P1) uses nodal demands
directly; giving item node v_i demand a_i and bucket node u_j demand 1
yields the identical argmin (same objective up to the constant
sum(a_i - 1)) on a graph of only 1 + m + k nodes. The DP and the QAOA
circuit already use weights, so this puts CP-SAT on the same footing.
B=15 is included under this formulation too, so every row is
like-for-like; the Part 1 table remains the leaf-formulation record.

## What is measured per (m, B)

- CP-SAT proof time and status on the weighted MILP (30-minute cap).
- DP time (5-minute cap) -- the structure-aware classical route.
- The QAOA circuit's size, which is 3*m^2 qubits with a gate count
  that depends only on (m, k), never on B: only Hamiltonian
  coefficients change. Stated in the CSV so the point is on record.

Random items (seed 0) drawn uniformly from (B/4, B/2) and walked to
sum exactly m*B. At large B such instances are almost surely
NO-instances, and proving NO is the exhaustive case -- no hand-built
perturbation needed.
"""

from __future__ import annotations

import csv
import random
import signal
import time
from pathlib import Path

from ortools.sat.python import cp_model

from partition_gadget import PartitionGadget

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_PATH = RESULTS_DIR / "large_b_hardness_sweep.csv"
CPSAT_TIME_LIMIT_S = 1800.0
DP_TIME_LIMIT_S = 300.0
WORKERS = 8
M_VALUES = (4, 5, 6)
B_VALUES = (15, 1000, 10000)
SEED = 0


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def random_items(m: int, B: int, seed: int) -> tuple[int, ...]:
    rng = random.Random(seed)
    lo, hi = B // 4 + 1, (B - 1) // 2
    k = 3 * m
    items = [rng.randint(lo, hi) for _ in range(k)]
    target = m * B
    while sum(items) != target:
        i = rng.randrange(k)
        if sum(items) < target and items[i] < hi:
            items[i] += 1
        elif sum(items) > target and items[i] > lo:
            items[i] -= 1
    return tuple(items)


def solve_weighted(gadget: PartitionGadget):
    """Weighted-demand radial-reconfiguration MILP on the leafless gadget:
    root 0, bucket nodes 1..m (demand 1), item nodes m+1..m+k (demand
    a_i). Two flow systems on the same edge set: unit-demand for tree
    feasibility (verbatim pattern from run_hardness_sweep.py), and
    weighted-demand whose per-edge magnitude is the successor weight
    the objective squares."""
    m, k = gadget.m, gadget.k
    root = 0
    buckets = list(range(1, 1 + m))
    items_nodes = list(range(1 + m, 1 + m + k))
    n = 1 + m + k

    edges = [(root, u) for u in buckets] + [(u, v) for v in items_nodes for u in buckets]
    n_edges = len(edges)
    weight = [0] + [1] * m + list(gadget.items)
    total_w = sum(weight)

    model = cp_model.CpModel()
    x = [model.new_bool_var(f"x{i}") for i in range(n_edges)]

    incident: list[list[tuple[int, bool]]] = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        incident[u].append((i, True))
        incident[v].append((i, False))

    def flow_system(cap: int, demand: list[int], root_supply: int, tag: str):
        fwd = [model.new_int_var(0, cap, f"{tag}f{i}") for i in range(n_edges)]
        bwd = [model.new_int_var(0, cap, f"{tag}b{i}") for i in range(n_edges)]
        for i in range(n_edges):
            model.add(fwd[i] <= cap * x[i])
            model.add(bwd[i] <= cap * x[i])
        for node in range(n):
            outflow, inflow = [], []
            for i, is_u in incident[node]:
                if is_u:
                    outflow.append(fwd[i])
                    inflow.append(bwd[i])
                else:
                    outflow.append(bwd[i])
                    inflow.append(fwd[i])
            net_out = sum(outflow) - sum(inflow)
            model.add(net_out == (root_supply if node == root else -demand[node]))
        return fwd, bwd

    flow_system(n - 1, [1] * n, n - 1, "u")
    c_fwd, c_bwd = flow_system(total_w, weight, total_w - weight[root], "w")
    model.add(sum(x) == n - 1)

    sq = []
    for i in range(n_edges):
        mag = model.new_int_var(0, total_w, f"mag{i}")
        model.add(mag == c_fwd[i] + c_bwd[i])
        s = model.new_int_var(0, total_w * total_w, f"sq{i}")
        model.add_multiplication_equality(s, [mag, mag])
        sq.append(s)
    model.minimize(sum(sq))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = CPSAT_TIME_LIMIT_S
    solver.parameters.num_workers = WORKERS
    t0 = time.perf_counter()
    status = solver.solve(model)
    elapsed = time.perf_counter() - t0
    obj = int(round(solver.objective_value)) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return solver.status_name(status), status == cp_model.OPTIMAL, elapsed, obj, n, n_edges


def time_dp(gadget: PartitionGadget):
    signal.alarm(int(DP_TIME_LIMIT_S))
    t0 = time.perf_counter()
    try:
        val = gadget.exact_optimum()
        return time.perf_counter() - t0, val, False
    except _Timeout:
        return time.perf_counter() - t0, None, True
    finally:
        signal.alarm(0)


def main() -> None:
    signal.signal(signal.SIGALRM, _alarm)
    RESULTS_DIR.mkdir(exist_ok=True)
    fieldnames = ["m", "B", "k", "n_qubits", "graph_nodes", "graph_edges", "cpsat_status", "cpsat_proved_optimal",
                  "cpsat_time_s", "cpsat_objective", "dp_time_s", "dp_timed_out", "dp_objective_weighted",
                  "matches", "items"]
    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for m in M_VALUES:
            for B in B_VALUES:
                items = random_items(m, B, SEED)
                g = PartitionGadget(items=items, m=m, B=B)
                print(f"=== m={m} B={B} k={g.k} n_qubits={g.n_qubits} ===", flush=True)

                dp_time, dp_val, dp_to = time_dp(g)
                # leaf-formulation DP value -> weighted objective: subtract the leaf edges' constant
                dp_weighted = None if dp_val is None else dp_val - sum(a - 1 for a in items)
                print(f"  DP: {'TIMEOUT' if dp_to else f'{dp_time:.2f}s'}  objective={dp_weighted}", flush=True)

                status, opt, cp_time, cp_obj, gn, ge = solve_weighted(g)
                matches = None if (cp_obj is None or dp_weighted is None) else (cp_obj == dp_weighted)
                print(f"  CP-SAT: {status:<10} proved={opt}  {cp_time:8.1f}s  objective={cp_obj}  "
                      f"matches_dp={matches}", flush=True)

                writer.writerow({
                    "m": m, "B": B, "k": g.k, "n_qubits": g.n_qubits, "graph_nodes": gn, "graph_edges": ge,
                    "cpsat_status": status, "cpsat_proved_optimal": opt, "cpsat_time_s": round(cp_time, 2),
                    "cpsat_objective": cp_obj, "dp_time_s": round(dp_time, 2), "dp_timed_out": dp_to,
                    "dp_objective_weighted": dp_weighted, "matches": matches, "items": " ".join(map(str, items)),
                })
                f.flush()
    print(f"\nwrote {OUT_PATH}")


if __name__ == "__main__":
    main()
