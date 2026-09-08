"""Verification gates for `cycle_surrogate.py`, on the real IEEE 33-bus
feeder with pandapower's real loads and line resistances. Each gate is a
claim the module's docstring makes; none is taken on the docstring's
word.

1. LOCALITY: a single exchange changes flow ONLY on its own cycle's
   edges. Checked for every (tie, opened edge) -- every nonzero delta
   must land on the cycle. This is the fact the whole surrogate rests on.
2. SINGLE-EXCHANGE EXACTNESS: surrogate == exact for every single
   exchange (true by construction -- the delta table IS that exact
   flow -- so this is a self-consistency check, and it must be 0).
3. DOUBLE-EXCHANGE ERROR: for every pair of exchanges on two different
   ties that still yields a spanning tree, compare surrogate loss to
   exact loss. Split by whether the two cycles share an edge. Disjoint
   pairs MUST be exact; overlapping pairs are where the first-order
   approximation is actually approximating, and their error is the
   number that decides whether this surrogate is usable.
4. OVERLAP STRUCTURE: which ties' cycles share edges, and how many
   joint choices disconnect the network. This is the measured size of
   the coverage caveat.

Writes results/cycle_surrogate_ieee33.csv (one row per double exchange).
"""

from __future__ import annotations

import csv
from itertools import combinations
from pathlib import Path

from cycle_surrogate import build_surrogate, exchange_tree, load_ieee33_physical, signed_flows

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_PATH = RESULTS_DIR / "cycle_surrogate_ieee33.csv"
TOL = 1e-9


def main() -> None:
    feeder = load_ieee33_physical()
    s = build_surrogate(feeder)
    g = feeder.graph
    print(f"IEEE33: {g.n_nodes} buses, {g.n_edges} lines, {len(feeder.tie_idx)} ties, "
          f"total load {sum(feeder.load_mw):.3f} MW, base loss {s.loss(s.base_flow):.4f}")
    print(f"cycle lengths (edges incl. tie): { {t: len(c) for t, c in s.cycles.items()} }")
    print(f"one-hot qubits if each cycle is a register: {sum(len(c) for c in s.cycles.values())}")

    # --- gate 1: locality
    violations = 0
    for t, cyc in s.cycles.items():
        cyc_set = set(cyc)
        for opened, d in s.delta[t].items():
            off_cycle = [e for e, v in d.items() if e not in cyc_set and abs(v) > TOL]
            if off_cycle:
                violations += 1
                print(f"  LOCALITY VIOLATION tie={t} opened={opened}: nonzero delta off-cycle on {off_cycle}")
    n_single = sum(len(d) for d in s.delta.values())
    print(f"gate 1 (locality): {n_single} single exchanges, {violations} off-cycle violations "
          f"-> {'OK' if violations == 0 else '*** FAILED ***'}")

    # --- gate 2: single-exchange exactness
    max_err = 0.0
    for t, d in s.delta.items():
        for opened in d:
            exact = signed_flows(feeder, exchange_tree(feeder, {t: opened}))
            sur = s.flows({t: opened})
            max_err = max(max_err, max(abs(exact[e] - sur[e]) for e in range(g.n_edges)))
    print(f"gate 2 (single exact): max |exact - surrogate| flow = {max_err:.2e} "
          f"-> {'OK' if max_err < TOL else '*** FAILED ***'}")

    # --- gate 4 first (needed to label gate 3 rows): overlap structure
    overlaps = s.overlapping_pairs()
    overlap_set = {(t1, t2) for t1, t2, _ in overlaps}
    print(f"gate 4 (overlap): {len(overlaps)} of {len(list(combinations(s.cycles, 2)))} tie pairs share edges:")
    for t1, t2, shared in overlaps:
        print(f"    ties {t1},{t2} share {len(shared)} edge(s)")

    # --- gate 3: double exchanges
    rows = []
    disconnected = 0
    for t1, t2 in combinations(sorted(s.cycles), 2):
        for o1 in s.delta[t1]:
            for o2 in s.delta[t2]:
                choices = {t1: o1, t2: o2}
                tree = exchange_tree(feeder, choices)
                try:
                    exact = signed_flows(feeder, tree)
                except ValueError:
                    disconnected += 1
                    continue
                sur = s.flows(choices)
                exact_loss, sur_loss = s.loss(exact), s.loss(sur)
                rows.append({
                    "tie1": t1, "open1": o1, "tie2": t2, "open2": o2,
                    "cycles_overlap": (t1, t2) in overlap_set,
                    "exact_loss": exact_loss, "surrogate_loss": sur_loss,
                    "abs_err": abs(exact_loss - sur_loss),
                    "rel_err": abs(exact_loss - sur_loss) / exact_loss if exact_loss else 0.0,
                    "max_flow_err": max(abs(exact[e] - sur[e]) for e in range(g.n_edges)),
                })

    def summarize(label, sel):
        if not sel:
            print(f"    {label}: none")
            return
        print(f"    {label}: n={len(sel)}  max rel loss err={max(r['rel_err'] for r in sel):.3%}  "
              f"mean={sum(r['rel_err'] for r in sel) / len(sel):.3%}  "
              f"max flow err={max(r['max_flow_err'] for r in sel):.4f} MW")

    print(f"gate 3 (double exchanges): {len(rows)} valid trees, {disconnected} joint choices disconnect the network")
    disjoint = [r for r in rows if not r["cycles_overlap"]]
    overlapping = [r for r in rows if r["cycles_overlap"]]
    summarize("disjoint cycles (must be exact)", disjoint)
    summarize("overlapping cycles (first-order)", overlapping)
    if disjoint and max(r["abs_err"] for r in disjoint) > TOL:
        print("    *** FAILED: disjoint-cycle double exchange is not exact ***")

    RESULTS_DIR.mkdir(exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {OUT_PATH}")


if __name__ == "__main__":
    main()
