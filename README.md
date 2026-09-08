# A feasibility-preserving QAOA mixer for grid reconfiguration — and a provably hard family of the same problem

## Two questions, one problem

Quantum optimization has a feasibility problem. On constrained problems,
QAOA's difficulty usually isn't the objective — it's that most of the
states the algorithm explores aren't even allowed. The usual fix,
penalty terms, spends precious circuit depth punishing infeasible states
instead of searching among feasible ones, and still can't promise the
answer it returns is valid.

The cleaner idea is a **mixer that only ever moves between feasible
states**, so infeasibility never enters the search at all. The
difficulty is building one cheaply enough to run.

This repo does that for one real, practically important constraint:
keeping an electrical distribution network **radial** — every customer
connected, no loops — as switches open and close. And it asks two
questions about it, in order:

1. **Can the mixer be made cheap enough for real hardware?** Yes — at
   two tiers, one for fault-tolerant machines and one that meets a
   500-CX budget on every instance tested, validated on real published
   grids.
2. **Is this problem actually hard enough to need quantum?** For real
   grids, honestly, no — a classical solver walks them. But there is a
   provably hard *family* of the same problem. This repo builds
   instances from it, measures a general-purpose solver giving up on
   them, and shows the circuit for them is small — with one caveat,
   stated where it arises, that sets up the next experiment.

The second question is the one most quantum-optimization work skips.

## The problem, and the move that keeps it feasible

A distribution feeder has more switches than it strictly needs: a
normally-closed backbone, plus a handful of normally-open **tie**
switches, closed only to reroute power around a fault. Whatever the
switches do, the network must stay radial.

Each switch is a qubit (1 = closed). A configuration is feasible when
the closed switches form a **spanning tree** — every node reached, no
loops.

![A feasible configuration (spanning tree) vs. an infeasible one (a loop)](results/illustration_feeder_problem.png)

The mixer's move is a **basis exchange**: close one switch, open
another, and land on a different tree. Matroid theory guarantees a
chain of such moves can reach any feasible configuration from any other
— so the whole mixer is just one small circuit term per possible
exchange.

![One basis-exchange move: before, the move itself, and after](results/illustration_basis_exchange_move.png)

**Here is the catch.** *Which* switch has to open depends on the rest
of the current configuration. So each exchange term must be conditioned
on other qubits — specifically, on the loop the closing switch would
create, its **fundamental cycle**. A tie's *range* is how far apart its
two ends sit on the backbone, and that is exactly the length of this
loop. A **short-range** tie closes a small loop and needs a cheap
condition. A **long-range** tie closes a loop that winds across much of
the network, and needs a gate controlled on many qubits — expensive.

![A tie edge's fundamental cycle, short vs long](results/illustration_fundamental_cycle.png)

Real tie switches are long-range on purpose: they exist to reroute
power around a fault, so they deliberately link *distant* parts of a
feeder. Measured on published topologies, real ties span 33–45% of the
network's diameter. Real topology is the expensive case, and everything
in the next section is about paying that cost down.

## Two ways to make it cheap

**Tier 1 — whole graph, for fault-tolerant hardware.** Build the mixer
on the entire network, but cap how many qubits any single term may
condition on. Where the cap bites, the term becomes slightly imprecise:
a small, *measured* fraction of probability can leak into infeasible
states. The search that chooses each term's condition has to weigh
circuit cost against that leakage explicitly — letting it minimize
leakage alone drives cost up by orders of magnitude. Done right, this
gives a complete, fully connected mixer at a cost fault-tolerant
hardware can absorb.

**Tier 2 — zones, for today's hardware.** Cut the network into small
zones along its tie lines, build a Tier 1 mixer inside each zone, and
stitch them together with one small mixer over the contracted "zone
graph." The mathematics of spanning trees guarantees the stitched
result is a valid tree, so the decomposition itself cannot leak. Then
enforce a cost budget: transpile each piece, check its real gate count,
and split or retune any piece that comes in over the line.

Two honest limits. Tier 2 searches only trees that stay connected
*inside* every zone — a subset of all feasible configurations, so an
optimum that doesn't respect the zone boundaries is out of reach. And
leakage inside a zone is prevented in practice by keeping zones small,
not by construction. The 500-CX budget is a target the algorithm
recurses toward; the finding is that on every instance tested it got
there without leaking and without ever invoking its cost-for-safety
fallback.

![A graph partitioned into zones, plus the contracted assembly problem](results/illustration_decomposition.png)

*(A third, simpler construction — exact conditioning with no cap — is
the baseline both tiers grow from. It is not a tier: when no small
exact condition exists it silently drops the exchange, and on real
topology that leaves the mixer disconnected. `docs/circuit-validity.md`
traces that failure.)*

## Results

Where the tiers land is judged against a simple yardstick. With
published two-qubit error rates, a circuit's chance of running cleanly
is roughly `(1 − p)^N_CX`; 500 two-qubit gates is about where today's
hardware stops being viable:

| CX gates | trapped-ion (p=0.001) | superconducting (p=0.005) |
|---|---|---|
| 100 | 90% | 61% |
| 500 | 61% | 8% |
| 5,000 | 0.7% | ~0 |

**Why two kinds of experiment.** Published feeder topologies are scarce
— this repo has two, at 15 and 33 buses. They can show the tiers work
on real structure, but not whether they *keep* working as networks
grow toward realistic sizes. A synthetic ladder can grow, but it only
means something if it is built to resemble real feeders in the ways
that actually drive cost. So the ladder answers *does it scale?*, the
real feeders answer *does it transfer?*, and each is convincing only
because of the other.

**The synthetic ladder** runs from 10 to 150 nodes, three random
feeders per size, matched to real feeders in the two ways that matter:
every tie is long-range, and the number of ties grows only
logarithmically with network size. That second fact comes from real
data — thin data, three published feeders being all that exist at this
scope, but consistent: from 15 to 179 buses the ratio of ties to buses
falls about 6×, which is what logarithmic growth looks like. Larger real grids
don't get proportionally more redundancy, just a little more.

![Synthetic feeders with real-topology tie statistics: the two tiers (mean over 3 seeds, bars = min–max)](results/construction_progression_plot.png)

At the largest size (150 nodes, 156 qubits), Tier 1 costs 10,712 CX at
depth ~18,600; Tier 2 costs 337 CX at depth ~650. That is under the
500-CX line, but by the yardstick above it is viable on trapped-ion
hardware and marginal on superconducting. Cost alone isn't the whole
question, though: Tier 1 gets its cost down by capping conditions, so
it has to be asked how much exactness the cap gave up. The answer is
mean leakage of up to 9%, with roughly a third of starting
configurations leaking at all. Tier 2, asked the same question, gave up
nothing on any instance tested. Cheaper *and* exact in practice, not a
tradeoff.

**The real feeders** are the transfer test: topology nobody designed to
be convenient. Both tiers were built on the CIGRE MV benchmark (15
buses, 3 ties) and the IEEE 33-bus feeder (33 buses, 5 ties), five
random seeds each:

![Real networks: the two deployment tiers](results/real_network_comparison_plot.png)

![Where each real network's two tiers land relative to NISQ feasibility](results/real_nisq_feasibility_plot.png)

Tier 2 lands comfortably inside today's hardware budget on both
networks — 64 CX at depth 147 on 17 qubits, 132 CX at depth 272 on 37
qubits — and its number is identical across every seed: predictable,
not just cheap. Tier 1 stays complete and fully connected on both, at a
cost only fault-tolerant hardware could absorb. But fault tolerance
removes gate error, not leakage, and Tier 1's leakage is permanent: on
IEEE33, 43% of starting configurations leak. There is, in other words,
no exact *and* complete whole-graph construction that works on real
topology — the exact one comes out disconnected, and Tier 1 leaks.

The full account — including the intermediate constructions, the
short-range control condition, and the safety measurements not shown
here — is in `docs/scaling-ladder-and-decomposition.md` and
`docs/bounded-witness-mixer.md`.

## The second question: does it matter?

So the mixer runs, cheaply, on real grids. The obvious next question is
whether a quantum computer was ever needed. For the real grids above —
no. A feeder has only a handful of tie switches, so its set of feasible
configurations is small, and a classical solver simply enumerates it.

The cheapest way to make a real grid harder would be to keep its
topology and make the objective nastier. That was tried first, on the
IEEE 33-bus feeder, with a random quadratic cost and then with a
physically motivated congestion cost. Both were solved to proven
optimality in under a third of a second: with so few feasible
configurations, the objective barely matters.

Hardness, then, has to come from the combinatorics, not the topology.
So this repo takes a published NP-hardness proof for exactly this
problem — radial reconfiguration under a loss-minimizing objective
(Khodabakhsh et al., arXiv:1711.03517) — and builds the hard instance
it describes. Three things were measured, not assumed, each answering a
different objection:

- **"A classical solver would just handle it."** The honest measure of
  classical hardness is not how fast a solver *finds* a good tree —
  heuristics do that quickly on almost anything — but how fast it can
  *prove* the tree is optimal. CP-SAT's proof time grows from a few
  seconds on a 65-node instance to a 30-minute timeout, proof
  unfinished, at 177 nodes. Every point is cross-checked against an
  independent exact calculation — and that calculation is the caveat.
  It is a dynamic program that exploits the very structure the
  hardness proof exposes, and at the fixed bucket size used here it
  runs in milliseconds. So what fails is a general-purpose solver on
  the natural formulation, not every classical route. Genuine hardness
  needs the bucket size to grow with the instance — which inflates
  every classical route but leaves the quantum circuit's size untouched.
  That is the next experiment.
- **"Then the quantum circuit must be enormous."** Tier 1 and Tier 2
  can't even be started on this instance: both begin by enumerating the
  feasible set, and here that set is astronomically large — which is
  precisely *why* the instance is hard. The way in is the hardness
  proof itself. It reveals that almost every edge of the instance is
  forced, and the few free choices are independent of one another.
  Independent choices need no conditioning, so the mixer needs no
  witness search at all — and 48 qubits with about 800 two-qubit gates
  cover the 65-node instance where the classical solver already takes
  seconds.
- **"A classical computer could just simulate that circuit."** A small
  circuit only matters if it can't be shortcut classically, and
  tensor-network (MPS) methods are the strongest classical tool for
  circuits like this one. MPS misses the exact answer by about 8% on
  the smallest instance it can be checked against, and its cost grows
  roughly 70× over a range where the circuit itself grows about 14× —
  the excess is entanglement, not size. This was one simulator at
  modest bond dimension, though; a stronger tensor-network attempt is
  the obvious next test.

Notice what just happened: the hard instance needed a *simpler* mixer
than the real grids did. That is not a contradiction — it is the key to
how the two halves of this repo fit together. Real grids are
*structurally* hard to stay feasible on — long loops, expensive
conditions — but *combinatorially* easy to solve. The
hard family is the mirror image: structurally trivial, combinatorially
hard in general. Each construction handles the axis its problem
actually has.
Whether an instance exists that is hard on both axes at once is an open
question this repo does not answer.

Full account: `docs/hard-instance-case-study.md`.

## What is and isn't shown

Question 1 shows the mixer is correct and cheap enough for both
hardware tiers. It does not show a quantum algorithm beating a
classical one: there is no objective, no QAOA run, and no test of the
boundary-coupling loop a real decomposed optimization would need.

Question 2 adds an objective, a real solver comparison, and a
simulation comparison — but on purpose-built instances that a
structure-aware classical algorithm still solves quickly at the bucket
size used, and it stops short of demonstrating advantage. Solution
quality at the genuinely hard sizes is unmeasured, because both
classical ways of checking it hit walls first. The case study says so,
with numbers.

Nothing broader is claimed for either. This is a measurement study,
not a production mixer-compilation library; `methodology.md` has the
precise boundary.

## How to reproduce

```
python -m venv .venv
source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt
```

**Question 1 — the mixer:**

```
python scripts/verify_correctness.py && python scripts/verify_leakage_trace.py
# the exact construction is leak-free where it applies; the Tier 1
# leakage tooling is verified against an exact reference.

python scripts/run_real_networks_hierarchical.py
# the real-network figures above: both tiers, both networks, 5 seeds each.

python scripts/plot_results_figures.py && python scripts/plot_illustrations.py
# every figure in this README, from already-committed data.
```

**Question 2 — the hard instance:**

```
python scripts/verify_partition_mixer.py
# exact-unitary proof that the case study's mixer preserves feasibility,
# plus an exact check of its cost Hamiltonian against the true objective.

python scripts/measure_partition_mixer.py
# transpiled gate counts and depth for the 48-qubit hard instance.

python scripts/run_hardness_sweep.py       # CP-SAT hardness curve (slow: hours)
python scripts/run_mps_scaling_check.py    # the MPS comparison (~15 min)
```

The circuit-construction scripts are deterministic (fixed seeds) and
reproduce the committed `results/` files, modulo library versions. The
CP-SAT sweep runs 8 parallel workers and its *timings* vary run to run
by up to ~2×; its objective values do not. Further investigations use
the scripts indexed in `docs/repository-map.md`.

## Repository layout

- `methodology.md` — graph generation, move generation, and exactly what
  "gate count" and "depth" mean.
- `docs/mixer-construction.md` — matroid theory, exact circuit
  derivation, correctness arguments.
- `docs/circuit-validity.md` — why the exact construction fails on real
  topology, and the decomposition fix.
- `docs/bounded-witness-mixer.md` — Tier 1 in full: the witness cap,
  the leakage measurement, and the cost-aware search.
- `docs/scaling-ladder-and-decomposition.md` — the synthetic ladder,
  Tier 2 in full, safety measurements, and the real-network validation.
  Everything in Results traces back to it.
- `docs/hard-instance-case-study.md` — question 2 in full: the
  instance, the solver curve, the small circuit, and the MPS check.
- `scripts/` and `results/` — one script per measurement, one CSV/plot
  pair per script, all generated, none hand-edited. Full index:
  `docs/repository-map.md`.

## License

Apache License 2.0 — see `LICENSE`.
