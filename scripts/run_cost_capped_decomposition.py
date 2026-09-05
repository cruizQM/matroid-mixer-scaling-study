"""Targets a HARD guarantee the earlier decomposition work never aimed
for: every leaf subproblem's CX cost <= CX_THRESHOLD (500). Earlier
decomposition (flat and hierarchical) used a fixed or density-scaled
`target_zone_size` chosen BLIND to actual cost -- most zones end up
cheap, but a handful don't (one zone: 3,550 of a 3,566-CX subtotal), and
nothing in the earlier design re-acts to that after the fact.

## The fix: measure, don't guess

Build each subproblem, TRANSPILE it, and check its real CX cost
directly. If it's already under threshold, stop -- it's a leaf. If not,
partition IT further (same `zone_decomposition` machinery) and recurse on
each of ITS zones + assembly. Falls back to raising `cost_alpha` (trading
safety for cost, `run_cost_aware_scaling_ladder_aggressive.py`'s
mechanism, but bounded here -- never past the point of inertness, see
that script's own findings) only when a subproblem is too small to
partition further (`n_nodes <= MIN_PARTITION_SIZE`) and is STILL over
threshold -- a graph that small has no more structure left to exploit,
only cost pressure.

## Two refinements, found by asking "is this actually the best we can do?"
of the real-network results (not just "does it meet the threshold")

1. **Don't stop at the first `cost_alpha` that passes.** The original
   version tried only 0.01 and accepted it the instant it met threshold.
   On CIGRE MV that gave 274 CX at 0.985 mean feasible mass -- but
   `cost_alpha=0.2` gives 64 CX at a PERFECT 1.0, strictly better in both
   dimensions, simply never tried because 0.01 already "worked."
   `_best_under_threshold` now sweeps all of `FALLBACK_COST_ALPHAS` and
   keeps the cheapest passing result, at every leaf, not just the
   irreducible-fallback ones.
2. **Try more than one granularity, keep whichever total is cheaper.** The
   original version always jumped straight to `target_zone_size // 2` the
   instant a graph failed to meet threshold, meaning the granularity this
   repo's own flat decomposition uses directly (`target_zone_size=8`) was
   never tried once decomposition was needed at all -- only its half. On
   IEEE33 this mattered a lot: halving immediately produced a finer,
   HARDER-to-condition 4-node assembly-of-assembly core (508 CX) than the
   one un-halved `target_zone_size=8` produces (132 CX) -- same size,
   different specific graph, because which zones get contracted together
   differs. A first attempt just swapped the rule to "always try the
   un-halved size first" -- checked against the FULL synthetic ladder, not
   just the two real networks, that made 12 of 30 seeds WORSE (up to 15x),
   because shrinking zones sometimes helps and sometimes hurts depending
   on where a graph's real density sits -- there's no universally-better
   fixed rule. `_try_granularity` now tries BOTH the current
   `target_zone_size` and its half, fully recurses each, and
   `decompose_to_threshold` keeps whichever total is actually cheaper --
   the same "measure, don't guess" discipline this whole script is built
   on, applied to the granularity choice itself.

Writes results/cost_capped_decomposition_results.csv (per condition,
size, seed) and _summary.csv (whether every leaf actually met the
threshold, and how many levels/leaves it took).
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from graphs import enumerate_spanning_trees
from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones_by_size
from run_cost_aware_scaling_ladder import CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "cost_capped_decomposition_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "cost_capped_decomposition_summary.csv"

CX_THRESHOLD = 500
MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 2
BASE_TARGET_ZONE_SIZE = 8
MIN_TARGET_ZONE_SIZE = 2
MIN_PARTITION_SIZE = 4  # below this, splitting further isn't meaningful -- fall back to cost pressure instead
MAX_DEPTH = 8
FALLBACK_COST_ALPHAS = [0.01, 0.05, 0.2, 1.0, 5.0]  # escalating fallback if still over threshold after fallback recursion depth
N_STARTING_TREES_PER_SUBPROBLEM = 30
MAX_EXACT_ENUMERATION_COMBINATIONS = 200_000
WALK_STEPS = 300


def _build_and_measure(graph, seed: int, cost_alpha: float):
    from math import comb
    if graph.n_nodes < 2:
        return None, 0, 0, [1.0]
    target_size = graph.n_nodes - 1
    if graph.n_edges >= target_size and comb(graph.n_edges, target_size) <= MAX_EXACT_ENUMERATION_COMBINATIONS:
        trees = enumerate_spanning_trees(graph)
    else:
        trees = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    if not trees:
        return None, 0, 0, [1.0]
    truncated = build_truncated_witness_mixer(
        graph, trees, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
        adaptive=False, cost_alpha=cost_alpha,
    )
    construction = truncated.construction
    if not construction.terms:
        return construction, 0, 0, [1.0]
    qc = mixer_circuit(construction, beta=0.37)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    cx = tqc.count_ops().get("cx", 0)
    rng = np.random.default_rng(seed)
    n_check = min(N_STARTING_TREES_PER_SUBPROBLEM, len(trees))
    masses = [final_feasible_mass(graph, construction, random_spanning_tree(graph, rng), 0.37, sparse=True)
              for _ in range(n_check)]
    return construction, cx, tqc.depth(), masses


def _best_under_threshold(graph, seed: int):
    """Sweep FALLBACK_COST_ALPHAS and keep the CHEAPEST result that still
    meets CX_THRESHOLD -- not just the first one that does. The original
    version stopped at cost_alpha=0.01 the instant it passed, which left
    real headroom on the table: on CIGRE MV, 0.01 gives 274 CX at 0.985
    mean feasible mass, but 0.2 gives 64 CX at a PERFECT 1.0 -- cheaper
    and safer, simply never tried because 0.01 already "worked." If
    nothing meets the threshold, returns the cheapest attempt anyway (met=False)."""
    best = None
    met_any = False
    for alpha in FALLBACK_COST_ALPHAS:
        construction, cx, circ_depth, masses = _build_and_measure(graph, seed, cost_alpha=alpha)
        if cx <= CX_THRESHOLD:
            met_any = True
            if best is None or cx < best[1]:
                best = (construction, cx, circ_depth, masses)
        elif best is None:
            best = (construction, cx, circ_depth, masses)
    return (*best, met_any)


def decompose_to_threshold(graph, seed: int, target_zone_size: int = BASE_TARGET_ZONE_SIZE, depth: int = 0) -> dict:
    """Returns dict(cx, depth, masses, n_leaves, max_leaf_cx, all_leaves_met_threshold, max_depth_reached)."""
    construction, cx, circ_depth, masses = _build_and_measure(graph, seed, cost_alpha=0.01)

    if cx <= CX_THRESHOLD or graph.n_nodes <= MIN_PARTITION_SIZE or depth >= MAX_DEPTH:
        # Whether we're here because 0.01 already fits, or because we're
        # forced to stop splitting regardless of cost: sweep for the
        # cheapest cost_alpha that meets the threshold rather than
        # settling for whatever the first one gave.
        construction, cx, circ_depth, masses, met = _best_under_threshold(graph, seed)
        return dict(cx=cx, depth=circ_depth, masses=masses, n_leaves=1,
                    max_leaf_cx=cx, all_leaves_met_threshold=met, max_depth_reached=depth)

    # Still over threshold and splittable: which granularity to split at is
    # genuinely graph-dependent, not a fixed rule. A first attempt at
    # "always try the current target_zone_size before halving" looked
    # right on IEEE33 (matching flat decomposition's own 4-zone split,
    # 508 -> 132 CX) but, checked against the full synthetic ladder,
    # made 12 of 30 seeds WORSE (up to 15x) -- shrinking zones sometimes
    # helps and sometimes hurts depending on where the real density sits,
    # exactly as this repo's own hierarchical-decomposition findings
    # already say elsewhere. So: try BOTH the current target_zone_size and
    # its half, fully recurse each, and keep whichever total is actually
    # cheaper -- measure, don't guess, applied to the granularity choice
    # itself, not just to whether a given attempt meets threshold.
    candidates = sorted({target_zone_size, max(MIN_TARGET_ZONE_SIZE, target_zone_size // 2)}, reverse=True)
    results = [r for r in (_try_granularity(graph, seed, size, depth) for size in candidates) if r is not None]

    if not results:
        # No candidate granularity could even produce a real split (graph
        # too small/dense to divide at any tried size) -- irreducible,
        # same cost_alpha-escalation fallback as the MIN_PARTITION_SIZE case.
        construction, cx, circ_depth, masses, met = _best_under_threshold(graph, seed)
        return dict(cx=cx, depth=circ_depth, masses=masses, n_leaves=1,
                    max_leaf_cx=cx, all_leaves_met_threshold=met, max_depth_reached=depth)

    return min(results, key=lambda r: r["cx"])


def _try_granularity(graph, seed: int, size: int, depth: int):
    """Partition `graph` at target_zone_size=`size` and fully recurse on
    every zone + the assembly graph (children get `size // 2` as their OWN
    starting granularity, to try again from there). Returns None if this
    size can't produce a real split (n_zones <= 1) -- not a candidate."""
    zone_of = partition_zones_by_size(graph, size)
    n_zones = len(set(zone_of.values()))
    if n_zones <= 1:
        return None

    next_target = max(MIN_TARGET_ZONE_SIZE, size // 2)
    total_cx, total_depth, all_masses = 0, 0, []
    n_leaves, max_leaf_cx, all_met, max_depth_reached = 0, 0, True, depth

    for zid in sorted(set(zone_of.values())):
        zg = build_zone_subgraph(graph, zone_of, zid)
        sub = decompose_to_threshold(zg, seed, next_target, depth + 1)
        total_cx += sub["cx"]; total_depth += sub["depth"]; all_masses.extend(sub["masses"])
        n_leaves += sub["n_leaves"]; max_leaf_cx = max(max_leaf_cx, sub["max_leaf_cx"])
        all_met = all_met and sub["all_leaves_met_threshold"]
        max_depth_reached = max(max_depth_reached, sub["max_depth_reached"])

    assembly_graph, _ = build_assembly_graph(graph, zone_of, n_zones)
    sub = decompose_to_threshold(assembly_graph, seed, next_target, depth + 1)
    total_cx += sub["cx"]; total_depth += sub["depth"]; all_masses.extend(sub["masses"])
    n_leaves += sub["n_leaves"]; max_leaf_cx = max(max_leaf_cx, sub["max_leaf_cx"])
    all_met = all_met and sub["all_leaves_met_threshold"]
    max_depth_reached = max(max_depth_reached, sub["max_depth_reached"])

    return dict(cx=total_cx, depth=total_depth, masses=all_masses, n_leaves=n_leaves,
                max_leaf_cx=max_leaf_cx, all_leaves_met_threshold=all_met, max_depth_reached=max_depth_reached)


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)
    result = decompose_to_threshold(graph, seed)
    unsafe = sum(1 for m in result["masses"] if m < 1.0 - 1e-6)
    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits_whole_graph": graph.n_edges,
        "total_cx": result["cx"], "total_depth": result["depth"],
        "n_leaves": result["n_leaves"], "max_leaf_cx": result["max_leaf_cx"],
        "all_leaves_met_threshold": result["all_leaves_met_threshold"],
        "max_depth_reached": result["max_depth_reached"],
        "unsafe_rate": round(unsafe / len(result["masses"]), 4) if result["masses"] else 0.0,
        "mean_feasible_mass": round(float(np.mean(result["masses"])), 4) if result["masses"] else 1.0,
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            for seed in range(SEEDS_PER_SIZE):
                row = run_one(condition, n_nodes, seed)
                rows.append(row)
                print(f"  {row}", flush=True)

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    summary_rows = []
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            group = [r for r in rows if r["condition"] == condition and r["n_nodes"] == n_nodes]
            cx_vals = [r["total_cx"] for r in group]
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes,
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_max": max(cx_vals),
                "max_leaf_cx_max": max(r["max_leaf_cx"] for r in group),
                "all_seeds_met_threshold": all(r["all_leaves_met_threshold"] for r in group),
                "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
                "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
            })
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")
    for row in summary_rows:
        print(f"  {row}")


if __name__ == "__main__":
    main()
