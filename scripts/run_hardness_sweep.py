"""CP-SAT hardness sweep over the gadget graph (gadget_graph.py),
cross-checked against the independent exact DP
(`PartitionGadget.exact_optimum`) at every size -- never trust a solver
time without a ground-truth check alongside it.

Requires `ortools` (not otherwise a dependency of this repo -- add it if
running this script: `pip install ortools`). This script is the classical
half of the hard-instance case study (see docs/hard-instance-case-study.md):
it measures how CP-SAT's time to PROVE optimality grows as the gadget
scales, independent of and prior to anything about the QAOA circuit.
"""

import csv
import time
from pathlib import Path

from ortools.sat.python import cp_model

from gadget_graph import build_gadget_graph
from partition_gadget import PartitionGadget

TIME_LIMIT_S = 1800.0
WORKERS = 8


def solve(gadget: PartitionGadget):
    g = build_gadget_graph(gadget)
    n = g.n_nodes
    edges = g.edges
    m_edges = len(edges)
    root = g.root

    model = cp_model.CpModel()
    x = [model.new_bool_var(f"x{i}") for i in range(m_edges)]

    f_fwd = [model.new_int_var(0, n - 1, f"f_fwd{i}") for i in range(m_edges)]
    f_bwd = [model.new_int_var(0, n - 1, f"f_bwd{i}") for i in range(m_edges)]
    for i in range(m_edges):
        model.add(f_fwd[i] <= (n - 1) * x[i])
        model.add(f_bwd[i] <= (n - 1) * x[i])

    incident: list[list[tuple[int, bool]]] = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        incident[u].append((i, True))
        incident[v].append((i, False))

    for node in range(n):
        outflow = []
        inflow = []
        for i, is_u in incident[node]:
            if is_u:
                outflow.append(f_fwd[i])
                inflow.append(f_bwd[i])
            else:
                outflow.append(f_bwd[i])
                inflow.append(f_fwd[i])
        net_out = sum(outflow) - sum(inflow)
        if node == root:
            model.add(net_out == n - 1)
        else:
            model.add(net_out == -1)

    model.add(sum(x) == n - 1)

    flow_mag = [model.new_int_var(0, n - 1, f"fm{i}") for i in range(m_edges)]
    for i in range(m_edges):
        model.add(flow_mag[i] == f_fwd[i] + f_bwd[i])
    sq = [model.new_int_var(0, (n - 1) ** 2, f"sq{i}") for i in range(m_edges)]
    for i in range(m_edges):
        model.add_multiplication_equality(sq[i], [flow_mag[i], flow_mag[i]])
    model.minimize(sum(sq))

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = TIME_LIMIT_S
    solver.parameters.num_workers = WORKERS
    t0 = time.perf_counter()
    status = solver.solve(model)
    elapsed = time.perf_counter() - t0
    obj = int(solver.objective_value) if status in (cp_model.OPTIMAL, cp_model.FEASIBLE) else None
    return solver.status_name(status), status == cp_model.OPTIMAL, elapsed, obj, g


def run_case(label: str, gadget: PartitionGadget, writer) -> None:
    exact = gadget.exact_optimum()
    status, opt, elapsed, obj, g = solve(gadget)
    match = obj == exact
    print(f"[{label}] m={gadget.m} k={gadget.k}  nodes={g.n_nodes} edges={len(g.edges)}  "
          f"status={status:<10} proved_optimal={opt}  time={elapsed:8.1f}s  "
          f"{'OK' if match else f'*** MISMATCH: cpsat={obj} exact={exact} ***'}", flush=True)
    writer.writerow({
        "label": label, "m": gadget.m, "k": gadget.k, "n_nodes": g.n_nodes, "n_edges": len(g.edges),
        "status": status, "proved_optimal": opt, "wall_clock_s": round(elapsed, 3),
        "cpsat_objective": obj, "exact_optimum": exact, "matches_exact": match,
    })


def main() -> None:
    out_path = Path(__file__).resolve().parent.parent / "results" / "hard_instance_hardness_sweep.csv"
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "label", "m", "k", "n_nodes", "n_edges", "status", "proved_optimal",
            "wall_clock_s", "cpsat_objective", "exact_optimum", "matches_exact",
        ])
        writer.writeheader()
        # YES: fours=2m, sevens=m -- always partitionable (m buckets of
        # {4,4,7}). NO: fours=2m-1, sixes=3, sevens=m-2 -- found by
        # exhaustive multiset search to admit no valid 3-partition at this
        # size/range.
        for m in (4, 5, 6, 7, 8, 9, 10, 11):
            yes_items = tuple([4] * (2 * m) + [7] * m)
            no_items = tuple([4] * (2 * m - 1) + [6, 6, 6] + [7] * (m - 2))
            run_case(f"m{m}_yes", PartitionGadget(items=yes_items, m=m, B=15), writer)
            f.flush()
            run_case(f"m{m}_no", PartitionGadget(items=no_items, m=m, B=15), writer)
            f.flush()
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
