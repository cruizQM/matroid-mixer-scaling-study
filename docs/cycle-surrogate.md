# A physical cost surrogate for the general feeder case: exact where cycles are disjoint, and a measured failure where they are not

## The question this document answers

Question 1's mixer (README) has never had a phase separator to pair with. The physical objective for reconfiguration — line loss `Σ_e r_e f_e²`, or a congestion surrogate `Σ_e (f_e/cap_e)²` — depends on the flow `f_e`, and on a radial network the flow on a line is the total demand *downstream of it in the current tree*. That is a global property of the configuration. No low-order polynomial in the switch variables expresses it, so there has been no problem-informed cost Hamiltonian for real feeders — only the mixer.

This document asks whether a change of variables rescues that: parameterize a configuration not by switch bits but by *which edge of each tie's fundamental cycle is open*, linearize the flow along cycles, and see whether the result is a usable local objective on real topology. The answer is measured on the IEEE 33-bus feeder with pandapower's real loads and line resistances, over its *entire* configuration space, not a sample.

`cycle_surrogate.py`, `verify_cycle_surrogate.py`, `verify_cycle_surrogate_full.py`. Raw results: `results/cycle_surrogate_ieee33.csv` (every double exchange), `results/cycle_surrogate_ieee33_full.csv` (every configuration).

## The surrogate

A real feeder sits at a base tree `T₀` (backbone closed, ties open). Each tie `t` has a fundamental cycle `C_t` — the `T₀` path between its endpoints, plus `t`. A single exchange on that cycle (close `t`, open some `e ∈ C_t`) reroutes the subtree cut off by `e` through `t`.

**The fact everything rests on:** a single exchange changes the flow *only on the edges of its own cycle*. Every edge off `C_t` keeps its downstream set; every edge above the cycle's junction carries the same total. So `Δf_e(t, e_open)` is a fixed number, computed once from `T₀` and the real loads. With one choice per tie — which edge of its cycle is open, choosing `t` itself meaning "leave it" — the flow surrogate is

```
f_e(choices) ≈ f_e⁰ + Σ_t Δf_e(t, choice_t)
```

signed against a fixed per-edge orientation (a rerouted segment genuinely reverses; magnitudes do not superpose). The loss `Σ_e r_e f_e²` is then quadratic in one-hot choice variables: Z and ZZ terms across registers, coupled only where two cycles share an edge. That is the hard-instance Hamiltonian's shape (`docs/hard-instance-case-study.md`, Part 2) — "same bucket" coupling there, "shared cycle edge" coupling here.

Exact for edge-disjoint cycles and for any single exchange. First-order otherwise.

A second variant expands the *loss* rather than the flow to second order in "which exchanges happened", tabulating every single and double exchange exactly — still pairwise in the one-hot variables (`SecondOrderLoss`).

## What was checked, in order

Each is a claim the construction makes; none is taken on the derivation's word.

1. **Locality.** For every (tie, opened edge) — 59 single exchanges on IEEE33 — every nonzero flow change lands on that tie's own cycle. **0 violations.**
2. **Single-exchange exactness.** Surrogate equals exact flow on every single exchange to 4×10⁻¹⁶. (True by construction; a self-consistency check.)
3. **Overlap structure.** IEEE33's five ties have cycles of 10, 7, 15, 21 and 11 edges (64 one-hot qubits). **8 of the 10 tie pairs share edges** — one pair shares 9. Only two pairs are disjoint.
4. **Double exchanges** (1,170 valid trees; 164 joint choices disconnect the network):

| cycles | n | first-order surrogate, relative loss error | flow error |
|---|---|---|---|
| disjoint | 114 | **0.000% max** | 0 |
| overlapping | 1,056 | **346% max, 10.7% mean** | up to 3.26 MW — on a 3.7 MW system |

The disjoint rows are exact, as the theory says. The overlapping rows are not a correction term; they are the linearization breaking. When two cycles share edges, the subtree one exchange reroutes can contain the other cycle's path, and a `Δf` computed against `T₀` assumes a downstream set that no longer exists.

5. **The entire configuration space.** Five ties with cycles of those lengths give 10·7·15·21·11 = 242,550 one-open-edge-per-cycle configurations — small enough to compute the exact radial flow and loss for every one.

**144,932 of them (59.8%) disconnect the network.** The cycle-register parameterization is *mostly infeasible* on this feeder. Of the 97,618 valid trees, **30.1% contain a pair of choices that would disconnect the network on their own** (a third tie reconnects it) — so feasibility cannot be enforced by any pairwise penalty: a penalty on the pair leaks into feasible configurations. This was found the hard way: a first version of the second-order surrogate penalized such pairs at 10× base loss and reported mean errors of 34–245% that were mostly that penalty, not the expansion.

**The true optimum uses 4 simultaneous exchanges** (loss 13.48, 71% of the base 18.98). The good solutions live exactly where a low-order expansion is worst.

| simultaneous exchanges | n | first-order: max / mean | second-order (leak-free): max / mean |
|---|---|---|---|
| 0 | 1 | 0 / 0 | 0 / 0 |
| 1 | 59 | 0 / 0 | 0 / 0 |
| 2 | 1,170 | 346% / 10.7% | **0 / 0** |
| 3 | 9,562 | 540% / 27.0% | 264% / 7.9% |
| 4 | 35,346 | 560% / 42.2% | 321% / 19.4% |
| 5 | 51,480 | 574% / 51.8% | 320% / 31.2% |

**What a phase separator actually needs from an objective** is not small values but the right *landscape* — its low region should overlap the true low region:

| | Spearman rank correlation | surrogate's argmin, exact loss | true optimum's rank |
|---|---|---|---|
| first-order flow | 0.61 | 14.76 (9.5% above optimum) | #1,339 of 97,618 |
| second-order loss | **0.78** | **25.48 — worse than the base configuration** | #2,061 |

The second-order surrogate ranks the bulk of the landscape better and gets the bottom wrong: its additive pairwise corrections sum to a spurious minimum. An objective that steers toward a configuration worse than doing nothing is disqualified, whatever its correlation. The first-order surrogate is the opposite — worse correlation, but its minimum is a genuinely good tree. Neither is what a real-feeder phase separator should be built on.

## Verdict

- **The theory is right where it applies.** Locality holds exactly; disjoint-cycle configurations are exact to machine precision. A feeder whose ties have edge-disjoint fundamental cycles gets an exact, pairwise, problem-informed Hamiltonian and a true partition-matroid mixer from this construction.
- **Real feeders do not live there.** Long-range ties are placed for redundancy and their cycles overlap by design — on IEEE33, 8 of 10 pairs. In that regime the parameterization is mostly infeasible, pairwise penalties cannot repair it, and no low-order expansion of the objective tracks the optimum.
- **This is the same wall as before, reached from a different side.** The switch-bit formulation has no local objective because flow is global; the cycle-choice formulation has one only when cycles don't interact — i.e. when the global structure is absent. Changing variables did not remove the globality, it relocated it into the overlap.

## The way back is a construction this repo already has

Zone decomposition (Tier 2) cuts the network so that most cycles fall *inside* zones and the coupling between zones is carried by a small assembly graph. Inside a zone whose cycles are edge-disjoint — or that contains a single cycle — everything above is exact: the surrogate is the true loss, the one-hot registers form a real partition matroid, and the hard instance's witness-free mixer applies unchanged. The overlap that broke the whole-feeder surrogate becomes the assembly problem, which is small by construction and already what Tier 2 solves.

That is Tier 2 with a physical objective, and it is a design, not a result: zone choice would have to be driven by cycle disjointness rather than by tie-line cut alone, and the assembly-level coupling still needs an objective of its own. It is scoped here because the measurements above point at it specifically, not because any part of it has been built.

## Honest scope of this document

- **One feeder, one objective.** IEEE33 and lossless real-power loss. CIGRE MV (3 ties) was not run; a feeder with fewer, shorter ties would land in a friendlier regime, and the point of using IEEE33 was that its ties are the realistic, overlapping kind.
- **P-only, lossless, no voltages** — the real-power half of LinDistFlow. Reactive flow and voltage limits would add terms, not remove the globality problem.
- **The 60% infeasibility figure is a property of the parameterization on this topology**, not of the mixer: it is what any one-open-edge-per-cycle scheme inherits when cycles overlap. It is reported because it is the measured size of a caveat that would otherwise be stated as "some".
- **The exhaustive evaluation is possible only because IEEE33 is classically trivial** — 97,618 valid trees is nothing for a classical enumeration. That is consistent with everything in `docs/hard-instance-case-study.md`: this document is about giving question 1 an objective, not about hardness.
