"""Reads results/scaling_summary.csv (mean/min/max over several random
seeds per size) and produces the headline plot. Callable standalone (to
re-plot without rerunning the measurement sweep) or from
run_scaling_study.py."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"


def load_rows(csv_path: Path):
    with open(csv_path, newline="") as f:
        return list(csv.DictReader(f))


def make_plot(summary_csv_path: Path, out_path: Path) -> None:
    rows = load_rows(summary_csv_path)
    n_qubits = [int(r["n_qubits"]) for r in rows]
    cx_mean = [float(r["cx_mean"]) for r in rows]
    cx_min = [int(r["cx_min"]) for r in rows]
    cx_max = [int(r["cx_max"]) for r in rows]
    depth_mean = [float(r["depth_mean"]) for r in rows]
    depth_min = [int(r["depth_min"]) for r in rows]
    depth_max = [int(r["depth_max"]) for r in rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1.plot(n_qubits, cx_mean, marker="o", label="mean")
    ax1.fill_between(n_qubits, cx_min, cx_max, alpha=0.2, label="min-max across seeds")
    ax1.set_xlabel("qubits (= feeder edges)")
    ax1.set_ylabel("transpiled CX count")
    ax1.set_title("Mixer circuit CX count vs. instance size")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    ax2.plot(n_qubits, depth_mean, marker="o", color="darkorange", label="mean")
    ax2.fill_between(n_qubits, depth_min, depth_max, alpha=0.2, color="darkorange", label="min-max across seeds")
    ax2.set_xlabel("qubits (= feeder edges)")
    ax2.set_ylabel("transpiled circuit depth")
    ax2.set_title("Mixer circuit depth vs. instance size")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    fig.suptitle(
        f"Matroid basis-exchange mixer scaling, sparse feeder graphs "
        f"(n_qubits={n_qubits[0]}-{n_qubits[-1]}, several random seeds per size)"
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    csv_path = RESULTS_DIR / "scaling_summary.csv"
    if not csv_path.exists():
        print(f"{csv_path} not found -- run run_scaling_study.py first", file=sys.stderr)
        sys.exit(1)
    make_plot(csv_path, RESULTS_DIR / "scaling_plot.png")
    print(f"wrote {RESULTS_DIR / 'scaling_plot.png'}")
