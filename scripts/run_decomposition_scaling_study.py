"""Tests whether zone decomposition keeps mixer witness size bounded as
the REAL FEEDER SIZE grows, not just at the one fixed 33-bus instance
available. Only one real topology exists to validate against directly
(see `real_feeders.py`); this uses `graphs.py::generate_feeder_graph_long_range_ties`
-- a synthetic family that reproduces the specific failure mode found on
real data (long-range tie switches, unlike the original nearest-tie
synthetic model) -- at controllable, arbitrary network size, with a FIXED
target zone size (not a fixed zone count) so the per-subproblem qubit
budget doesn't grow with the network.

For each network size, reports:
  - "naive" whole-graph outlook: the longest fundamental cycle in the
    whole (undecomposed) graph -- a fast, cheap proxy for the witness
    size a whole-graph mixer construction would need (see
    docs/circuit-validity.md's real-feeder section for why this is the
    right proxy).
  - decomposed: max witness size actually found (and exactly verified
    leak-free) across every zone and the assembly problem.
  - **What drives circuit cost beyond witness size** (added after an
    initial run showed CX count varying up to 14x between same-size
    seeds despite near-identical witness size -- see
    docs/circuit-validity.md): term count and valid-pattern count per
    term, tracked explicitly rather than only the final CX number, plus
    a transpile at `optimization_level=3` alongside the repo's usual `1`
    to check how much of the gap a free transpiler setting closes
    (little, it turns out -- 1-3%). The actual fix,
    `_minimize_patterns` in `mixer.py` (Quine-McCluskey logic
    minimization of witness patterns), is implemented there and applies
    automatically to every run of this script via `build_matroid_mixer`
    -- the numbers this script now produces already include it. See
    `results/decomposition_scaling_results_before_minimization.csv` for
    the pre-fix numbers this same sweep produced, kept for the
    before/after comparison in docs/circuit-validity.md.

Seeds per size raised from an initial 3 to 12 for this reason: 3 points
cannot distinguish "flat with high per-seed variance" from "quietly
growing" -- see docs/circuit-validity.md for why that mattered.

Writes results/decomposition_scaling_results.csv (per size, per seed) and
results/decomposition_scaling_summary.csv (per size, mean/median/max),
plus results/decomposition_scaling_plot.png and
results/decomposition_cx_detail_plot.png.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from qiskit import transpile

from graphs import enumerate_spanning_trees, generate_feeder_graph_long_range_ties
from measure import TRANSPILE_BASIS
from mixer import build_matroid_mixer, fundamental_cycle_structure, mixer_circuit, verify_all_terms_no_leakage
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones_by_size

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "decomposition_scaling_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "decomposition_scaling_summary.csv"
PLOT_PATH = RESULTS_DIR / "decomposition_scaling_plot.png"
CX_DETAIL_PLOT_PATH = RESULTS_DIR / "decomposition_cx_detail_plot.png"

N_NODES_RANGE = [10, 20, 40, 80, 120]
K_TIES = 5  # matches the real IEEE 33-bus feeder's actual tie-switch count
TARGET_ZONE_SIZE = 8
SEEDS_PER_SIZE = 12


def naive_max_cycle_length(graph) -> int:
    """Cheap proxy for the witness size a whole-graph mixer would need:
    the size of the largest fundamental-cycle membership set any single
    edge participates in."""
    cycle_of = fundamental_cycle_structure(graph)
    return max((len(c) for c in cycle_of.values()), default=0)


def measure_subproblem_detailed(graph) -> dict:
    """Like run_zone_decomposition_validation.py's measure_subproblem,
    but also: transpiles at optimization_level=3 alongside the repo's
    usual 1, and reports term count / valid-pattern-count statistics
    directly, not just the resulting CX count -- see module docstring for
    why (CX count alone doesn't explain its own variance)."""
    trees = enumerate_spanning_trees(graph)
    if not trees:
        return dict(
            term_count=0, max_witness_size=0, total_controlled_gates=0, max_pattern_count=0,
            fully_connected=False, no_leakage_verified=True,
            cx_opt1=0, depth_opt1=0, cx_opt3=0, depth_opt3=0,
            num_spanning_trees=0,
        )

    construction = build_matroid_mixer(graph, trees)
    witness_sizes = [t.control_count for t in construction.terms]
    pattern_counts = [len(t.valid_patterns) for t in construction.terms]
    leak_free = verify_all_terms_no_leakage(construction) if construction.terms else True

    if construction.terms:
        qc = mixer_circuit(construction, beta=0.37)
        tqc1 = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        tqc3 = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=3)
        cx1, depth1 = tqc1.count_ops().get("cx", 0), tqc1.depth()
        cx3, depth3 = tqc3.count_ops().get("cx", 0), tqc3.depth()
    else:
        cx1 = depth1 = cx3 = depth3 = 0

    return dict(
        term_count=len(construction.terms),
        max_witness_size=max(witness_sizes, default=0),
        total_controlled_gates=sum(pattern_counts),
        max_pattern_count=max(pattern_counts, default=0),
        fully_connected=construction.fully_connected,
        no_leakage_verified=leak_free,
        cx_opt1=cx1, depth_opt1=depth1, cx_opt3=cx3, depth_opt3=depth3,
        num_spanning_trees=len(trees),
    )


def run_one(n_nodes: int, seed: int) -> dict:
    t0 = time.perf_counter()
    graph = generate_feeder_graph_long_range_ties(n_nodes, K_TIES, seed)
    naive_cycle = naive_max_cycle_length(graph)

    zone_of = partition_zones_by_size(graph, TARGET_ZONE_SIZE)
    n_zones = len(set(zone_of.values()))

    agg = dict(
        max_witness=0, max_cx1=0, max_depth1=0, max_cx3=0, max_depth3=0,
        total_cx1=0, total_cx3=0, total_controlled_gates=0, max_pattern_count=0,
        term_count_total=0, all_leak_free=True, all_connected=True,
    )

    def fold(row: dict) -> None:
        agg["max_witness"] = max(agg["max_witness"], row["max_witness_size"])
        agg["max_cx1"] = max(agg["max_cx1"], row["cx_opt1"])
        agg["max_depth1"] = max(agg["max_depth1"], row["depth_opt1"])
        agg["max_cx3"] = max(agg["max_cx3"], row["cx_opt3"])
        agg["max_depth3"] = max(agg["max_depth3"], row["depth_opt3"])
        agg["total_cx1"] += row["cx_opt1"]
        agg["total_cx3"] += row["cx_opt3"]
        agg["total_controlled_gates"] += row["total_controlled_gates"]
        agg["max_pattern_count"] = max(agg["max_pattern_count"], row["max_pattern_count"])
        agg["term_count_total"] += row["term_count"]
        if not row["no_leakage_verified"]:
            agg["all_leak_free"] = False
        if not row["fully_connected"] and row["num_spanning_trees"] > 1:
            agg["all_connected"] = False

    for zid in sorted(set(zone_of.values())):
        zg = build_zone_subgraph(graph, zone_of, zid)
        if zg.n_nodes < 2:
            continue
        fold(measure_subproblem_detailed(zg))

    assembly_graph, _ = build_assembly_graph(graph, zone_of, n_zones)
    fold(measure_subproblem_detailed(assembly_graph))

    elapsed = time.perf_counter() - t0
    return {
        "n_nodes": n_nodes,
        "seed": seed,
        "n_qubits_whole_graph": graph.n_edges,
        "n_zones": n_zones,
        "naive_max_cycle_length": naive_cycle,
        "decomposed_max_witness_size": agg["max_witness"],
        "term_count_total": agg["term_count_total"],
        "max_pattern_count": agg["max_pattern_count"],
        "total_controlled_gates": agg["total_controlled_gates"],
        "max_cx_opt1": agg["max_cx1"],
        "max_depth_opt1": agg["max_depth1"],
        "max_cx_opt3": agg["max_cx3"],
        "max_depth_opt3": agg["max_depth3"],
        "total_cx_opt1": agg["total_cx1"],
        "total_cx_opt3": agg["total_cx3"],
        "all_leak_free": agg["all_leak_free"],
        "all_connected": agg["all_connected"],
        "elapsed_s": round(elapsed, 2),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for n_nodes in N_NODES_RANGE:
        for seed in range(SEEDS_PER_SIZE):
            row = run_one(n_nodes, seed)
            rows.append(row)
            print(f"  {row}")

    fieldnames = list(rows[0].keys())
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {CSV_PATH}")

    def stats(group, key):
        vals = [r[key] for r in group]
        return {
            "mean": round(statistics.mean(vals), 1),
            "median": round(statistics.median(vals), 1),
            "max": max(vals),
        }

    summary_rows = []
    for n_nodes in N_NODES_RANGE:
        group = [r for r in rows if r["n_nodes"] == n_nodes]
        row = {"n_nodes": n_nodes, "n_seeds": len(group)}
        for key in (
            "naive_max_cycle_length", "decomposed_max_witness_size", "max_pattern_count",
            "total_controlled_gates", "max_cx_opt1", "max_cx_opt3",
        ):
            s = stats(group, key)
            row[f"{key}_mean"] = s["mean"]
            row[f"{key}_median"] = s["median"]
            row[f"{key}_max"] = s["max"]
        row["all_leak_free"] = all(r["all_leak_free"] for r in group)
        row["all_connected"] = all(r["all_connected"] for r in group)
        summary_rows.append(row)

    summary_fields = list(summary_rows[0].keys())
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")

    try:
        import matplotlib.pyplot as plt

        ns = [r["n_nodes"] for r in summary_rows]

        fig, ax1 = plt.subplots(figsize=(8, 5))
        ax1.plot(ns, [r["naive_max_cycle_length_mean"] for r in summary_rows], "o-", color="crimson",
                 label="whole-graph naive: max fundamental cycle length (proxy for required witness size)")
        ax1.plot(ns, [r["decomposed_max_witness_size_mean"] for r in summary_rows], "s-", color="steelblue",
                 label=f"decomposed (target zone size={TARGET_ZONE_SIZE}): max witness size across zones+assembly")
        ax1.set_xlabel("network size (n_nodes)")
        ax1.set_ylabel("qubits")
        ax1.set_title(f"Whole-graph vs. zone-decomposed witness requirement ({SEEDS_PER_SIZE} seeds/size)")
        ax1.legend(loc="upper left", fontsize=8)
        ax1.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(PLOT_PATH, dpi=150)
        print(f"wrote {PLOT_PATH}")

        fig2, (axa, axb) = plt.subplots(1, 2, figsize=(12, 5))
        axa.plot(ns, [r["max_cx_opt1_mean"] for r in summary_rows], "o-", color="darkorange", label="optimization_level=1 (mean)")
        axa.plot(ns, [r["max_cx_opt3_mean"] for r in summary_rows], "o-", color="seagreen", label="optimization_level=3 (mean)")
        axa.plot(ns, [r["max_cx_opt1_max"] for r in summary_rows], "o--", color="darkorange", alpha=0.4, label="opt=1 (max)")
        axa.plot(ns, [r["max_cx_opt3_max"] for r in summary_rows], "o--", color="seagreen", alpha=0.4, label="opt=3 (max)")
        axa.set_xlabel("network size (n_nodes)")
        axa.set_ylabel("max CX count (worst subproblem)")
        axa.set_title("Transpiled CX count: optimization_level 1 vs 3")
        axa.legend(fontsize=7)
        axa.grid(alpha=0.3)

        axb.plot(ns, [r["max_pattern_count_mean"] for r in summary_rows], "o-", color="purple", label="max valid patterns / term (mean)")
        axb.plot(ns, [r["max_pattern_count_max"] for r in summary_rows], "o--", color="purple", alpha=0.4, label="max valid patterns / term (max)")
        axb.set_xlabel("network size (n_nodes)")
        axb.set_ylabel("valid patterns per term")
        axb.set_title("What drives CX variance: valid-pattern count per term")
        axb.legend(fontsize=8)
        axb.grid(alpha=0.3)

        fig2.tight_layout()
        fig2.savefig(CX_DETAIL_PLOT_PATH, dpi=150)
        print(f"wrote {CX_DETAIL_PLOT_PATH}")
    except ImportError:
        print("matplotlib not available, skipped plots")

    all_ok = all(r["all_leak_free"] for r in summary_rows)
    print(f"\nSummary: all subproblems leak-free across every size tested = {all_ok}")
    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
