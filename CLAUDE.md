# CLAUDE.md — matroid-mixer-scaling-study

## Purpose
A feasibility-preserving QAOA mixer for grid reconfiguration — every state it explores is a radial (spanning-tree) network, via matroid basis-exchange moves conditioned on fundamental-cycle "witness" qubits — plus a provably hard family of the same problem. Qiskit reference implementation; `qaoa-guppy-feeder-mixer` ports it to Guppy/Selene. **This repo is PUBLIC on GitHub** (`cruizQM/matroid-mixer-scaling-study`) — never commit client data or secrets.

## Current status and open questions
- README restructured 2026-09-08 as two parallel questions; hard-instance case study merged 2026-09-07 (branch `hard-instance-partition-gadget` still exists locally, no upstream).
- Not shown (README "What is and isn't shown"): no working physical objective, no QAOA run for question 1, no boundary-coupling loop; question 2 instances are purpose-built and solution quality at hard sizes is unmeasured.
- Scoped, not built: Tier 2 with a physical objective inside each zone.

## How to run
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt       # this repo predates the uv convention
python scripts/verify_correctness.py && python scripts/verify_leakage_trace.py
python scripts/run_real_networks_hierarchical.py
python scripts/verify_cycle_surrogate.py && python scripts/verify_cycle_surrogate_full.py
python scripts/verify_partition_mixer.py
```
`methodology.md` defines the precise scope boundary.

## Key findings
- Tier 1 (whole graph, capped witness, measured leakage) for fault-tolerant hardware; Tier 2 meets a 500-CX budget on every instance tested; both validated on real published grids, 5 seeds each.
- Real tie switches span 33-45% of network diameter — the expensive long-range case by design.
- Cycle-basis surrogate objective: exact for disjoint cycles, quantified failure on IEEE33 (`docs/cycle-surrogate.md`).
- Hard family (3-PARTITION reduction): at large B every classical route tried fails on a measured instance (m=6, B=10,000 unproved at 30 min), with a 75-qubit circuit no bigger than the easy cases.

## Conventions
- Currently **pip + requirements.txt**, not uv (deviation from the standard convention). <!-- TODO: migrate to uv -->
- Figures and CSVs under `results/` back the README and stay tracked; since 2026-09-10 `.gitignore` excludes new `results/`, `*.csv`, etc. Large sweeps to **ClearML** or `matroid-mixer-scaling-study_artifacts/`. `results/cycle_surrogate_ieee33_full.csv` (8 MB) is tracked.
- Update `repo.yaml` when status or findings change.
