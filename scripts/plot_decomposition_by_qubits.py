"""Re-expresses `results/decomposition_scaling_summary.csv` (produced by
`run_decomposition_scaling_study.py`) with qubit count on the x-axis
instead of node count, so it's directly comparable to headline result 1's
plot (`results/scaling_plot.png`), which is also in qubits. Reads the
already-computed summary rather than re-running the sweep -- this is a
relabeling of existing results, not a new measurement.

`n_qubits = n_nodes - 1 + k_ties`, and `k_ties=5` is fixed throughout
`run_decomposition_scaling_study.py`'s sweep, so this is an exact,
deterministic conversion (`n_nodes + 4`), not an approximation.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
SUMMARY_CSV_PATH = RESULTS_DIR / "decomposition_scaling_summary.csv"
PLOT_PATH = RESULTS_DIR / "decomposition_scaling_by_qubits_plot.png"

K_TIES = 5  # must match run_decomposition_scaling_study.py
HEADLINE_1_MAX_QUBITS = 35  # results/scaling_plot.png's tested range ends here


def main() -> None:
    with open(SUMMARY_CSV_PATH, newline="") as f:
        rows = list(csv.DictReader(f))

    qubits = [int(r["n_nodes"]) - 1 + K_TIES for r in rows]
    naive = [float(r["naive_max_cycle_length_mean"]) for r in rows]
    decomposed = [float(r["decomposed_max_witness_size_mean"]) for r in rows]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.axvspan(0, HEADLINE_1_MAX_QUBITS, color="#4472C4", alpha=0.08)
    ax.axvline(HEADLINE_1_MAX_QUBITS, color="#4472C4", ls="dotted", lw=1.5)
    ax.text(
        HEADLINE_1_MAX_QUBITS - 2, max(naive) * 0.55, "headline result 1\ntested up to here",
        ha="right", va="center", fontsize=8, color="#4472C4",
    )
    ax.plot(qubits, naive, "o-", color="crimson",
            label="whole-graph naive: max fundamental cycle length (proxy for required witness size)")
    ax.plot(qubits, decomposed, "s-", color="steelblue",
            label="decomposed (target zone size=8): max witness size across zones+assembly")
    ax.set_xlabel("whole-graph qubit count (n_qubits = n_nodes - 1 + k_ties)")
    ax.set_ylabel("qubits")
    ax.set_title("Whole-graph vs. zone-decomposed witness requirement, by qubit count")
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(PLOT_PATH, dpi=150)
    print(f"wrote {PLOT_PATH}")


if __name__ == "__main__":
    main()
