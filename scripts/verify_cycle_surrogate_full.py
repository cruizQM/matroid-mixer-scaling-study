"""Measure both loss surrogates against the EXACT loss over the entire
configuration space of the real IEEE 33-bus feeder.

Five ties with fundamental cycles of 10, 7, 15, 21 and 11 edges give
10*7*15*21*11 = 242,550 one-open-edge-per-cycle configurations -- small
enough to compute the exact radial flow and loss for every one. So the
surrogate error is not sampled or estimated; it is the full landscape.

Reports, per number of simultaneous exchanges (0..5):
  - first-order flow surrogate (cycle_surrogate.Surrogate): rel. loss error
  - second-order loss surrogate (SecondOrderLoss): rel. loss error
and, for what a QAOA phase separator actually needs from an objective:
  - Spearman rank correlation with the exact loss over all valid configs
  - whether each surrogate's argmin is the true optimum, and the true
    optimum's rank under each surrogate.

Writes results/cycle_surrogate_ieee33_full.csv (one row per configuration).
"""

from __future__ import annotations

import csv
from itertools import product
from pathlib import Path

import numpy as np

from cycle_surrogate import SecondOrderLoss, build_surrogate, exchange_tree, load_ieee33_physical, signed_flows

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_PATH = RESULTS_DIR / "cycle_surrogate_ieee33_full.csv"


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx = np.argsort(np.argsort(x))
    ry = np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def main() -> None:
    feeder = load_ieee33_physical()
    s = build_surrogate(feeder)
    ties = sorted(s.cycles)
    L0 = s.loss(s.base_flow)
    # infeasible_penalty=0.0, deliberately. A first run used 10*L0 and
    # reported second-order mean errors of 34-245% at 3-5 exchanges; most
    # of that was the penalty, not the expansion: 30% of VALID
    # configurations contain a pair of choices that disconnects the
    # network on its own (a third tie reconnects it), so a pairwise
    # penalty leaks into feasible trees. With b=0 for those pairs the
    # error measured below is the expansion's alone -- the number this
    # script exists to report. Feasibility has to be enforced elsewhere
    # (mixer or decomposition), not by pairwise terms; see the doc.
    second = SecondOrderLoss(s, infeasible_penalty=0.0)
    print(f"second-order surrogate: {sum(len(a) for a in second.a.values())} linear terms, "
          f"{second.n_pairwise_terms()} nonzero pairwise terms, "
          f"{second.n_infeasible_pairs} pairs that disconnect on their own (b set to 0, see comment)")

    rows = []
    n_total = n_disconnected = 0
    for combo in product(*(s.cycles[t] for t in ties)):
        n_total += 1
        choices = dict(zip(ties, combo))
        n_active = sum(1 for t, o in choices.items() if o != t)
        try:
            exact = s.loss(signed_flows(feeder, exchange_tree(feeder, choices)))
        except ValueError:
            n_disconnected += 1
            continue
        rows.append({
            "n_active": n_active,
            "exact": exact,
            "first_order": s.loss(s.flows(choices)),
            "second_order": second.loss(choices),
            "choices": " ".join(f"{t}:{o}" for t, o in choices.items()),
        })
    print(f"{n_total} configurations; {n_disconnected} disconnect the network; {len(rows)} valid trees")

    exact = np.array([r["exact"] for r in rows])
    fo = np.array([r["first_order"] for r in rows])
    so = np.array([r["second_order"] for r in rows])
    act = np.array([r["n_active"] for r in rows])

    print(f"\n{'#exch':>5} {'n':>7} | {'1st-order rel err':^22} | {'2nd-order rel err':^22}")
    print(f"{'':>5} {'':>7} | {'max':>10} {'mean':>10} | {'max':>10} {'mean':>10}")
    for k in range(6):
        m = act == k
        if not m.any():
            continue
        e1 = np.abs(fo[m] - exact[m]) / exact[m]
        e2 = np.abs(so[m] - exact[m]) / exact[m]
        print(f"{k:>5} {m.sum():>7} | {e1.max():>10.3%} {e1.mean():>10.3%} | {e2.max():>10.3%} {e2.mean():>10.3%}")

    i_true = int(np.argmin(exact))
    print(f"\ntrue optimum: loss={exact[i_true]:.4f} (base {L0:.4f}, {exact[i_true] / L0:.1%} of base), "
          f"{act[i_true]} exchanges, choices={rows[i_true]['choices']}")
    for name, arr in (("first-order", fo), ("second-order", so)):
        i_sur = int(np.argmin(arr))
        rank_of_true = int((arr < arr[i_true]).sum()) + 1
        print(f"{name:>13}: spearman={spearman(exact, arr):.4f}  argmin is true optimum: {i_sur == i_true}  "
              f"true optimum ranked #{rank_of_true} of {len(rows)}  "
              f"(surrogate argmin's exact loss = {exact[i_sur]:.4f})")

    RESULTS_DIR.mkdir(exist_ok=True)
    with open(OUT_PATH, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT_PATH}")


if __name__ == "__main__":
    main()
