# Repository map

Full script-by-script index, grouped by what each one validates. Each
script's own docstring explains what it measures and why; the `docs/*.md`
file covering that topic (linked from the main README) has the full
narrative.

- **Exact construction**: `graphs.py`, `mixer.py`, `measure.py`,
  `run_scaling_study.py`, `plot.py`, `verify_correctness.py`.
- **Real topology & decomposition**: `real_feeders.py`,
  `run_real_feeder_validation.py`, `investigate_fundamental_cycles.py`,
  `zone_decomposition.py`, `run_zone_decomposition_validation.py`,
  `run_decomposition_scaling_study.py`.
- **Bounded-witness mixer**: `random_trees.py`, `truncated_mixer.py`,
  `leakage_trace.py`, `verify_leakage_trace.py`, `tie_density_sweep.py`,
  `run_bounded_witness_safety_survey.py`, `measure_truncated_mixer.py`,
  `truncated_witness_cap_sweep.py`, `truncated_witness_cap_sweep_longrange.py`,
  `truncated_mixer_search_refinement.py`.
- **Escalating ladder & decomposition**: `run_scaling_study_log_ties.py`,
  `run_cost_aware_scaling_ladder.py`, `run_cost_aware_scaling_ladder_aggressive.py`,
  `run_fixed_alpha_ladder.py`, `run_decomposed_cost_aware_ladder.py`,
  `run_best_of_both_ladder.py`, `run_hierarchical_decomposed_ladder.py`,
  `run_cost_capped_decomposition.py`, `run_real_networks_hierarchical.py`,
  `exact_construction_ladder_check.py`.
- **Figures**: `plot_illustrations.py` (explanatory diagrams, not
  measurements), `plot_results_figures.py` (the README's result figures,
  from already-committed CSVs, no re-measurement).

`results/` holds one CSV/plot pair per script above, all generated, none
hand-edited. `*_before_minimization.*` files are pre-fix numbers, kept for
the before/after comparison in `docs/circuit-validity.md`.
`illustration_*.png` are explanatory diagrams, not measurements.
