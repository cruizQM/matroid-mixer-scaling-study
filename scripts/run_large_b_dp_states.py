"""How many partial-sum states does the structure-aware DP have to
explore, as a function of m and B? This is the MECHANISM behind the
large-B result (run_large_b_sweep.py): B does not make the instance
bigger, it makes the DP's memo table bigger, because random large
weights produce partial sums that no longer collide.

Same instances as the sweep (same `random_items`, same seed), DP only,
same 5-minute cap. On a timeout the state count reached so far is
recorded -- a lower bound, and the informative number.

Writes results/large_b_dp_states.csv.
"""

from __future__ import annotations

import csv
import signal
import time
from pathlib import Path

from partition_gadget import PartitionGadget
from run_large_b_sweep import B_VALUES, DP_TIME_LIMIT_S, M_VALUES, SEED, random_items

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
OUT_PATH = RESULTS_DIR / "large_b_dp_states.csv"


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def main() -> None:
    signal.signal(signal.SIGALRM, _alarm)
    with open(OUT_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["m", "B", "k", "n_qubits", "dp_time_s", "dp_timed_out", "dp_states"])
        writer.writeheader()
        for m in M_VALUES:
            for B in B_VALUES:
                g = PartitionGadget(items=random_items(m, B, SEED), m=m, B=B)
                stats: dict = {}
                signal.alarm(int(DP_TIME_LIMIT_S))
                t0 = time.perf_counter()
                try:
                    g.exact_optimum(stats=stats)
                    timed_out = False
                except _Timeout:
                    timed_out = True
                finally:
                    signal.alarm(0)
                elapsed = time.perf_counter() - t0
                states = stats.get("states", 0)
                print(f"m={m} B={B:>6}  {'TIMEOUT' if timed_out else f'{elapsed:7.2f}s'}  states={states:,}"
                      f"{' (lower bound)' if timed_out else ''}", flush=True)
                writer.writerow({"m": m, "B": B, "k": g.k, "n_qubits": g.n_qubits, "dp_time_s": round(elapsed, 2),
                                 "dp_timed_out": timed_out, "dp_states": states})
                f.flush()
    print(f"\nwrote {OUT_PATH}")


if __name__ == "__main__":
    main()
