"""3-PARTITION -> radial network reconfiguration hardness gadget.

Based on Khodabakhsh, Sharma & Naeem, "Distribution System Reconfiguration
for Optimal Loss Minimization is NP-hard" (arXiv:1711.03517): the unit-demand
special case of their loss-minimization objective is

    cost(T) = sum_{e=(parent,child) in spanning tree T} (successors(e))^2

where successors(e) is the number of nodes in the subtree hanging below e.
They prove this strongly NP-hard by reduction from 3-PARTITION: given k=3m
positive integers a_1..a_k, each strictly between B/4 and B/2, summing to
exactly m*B, build a tree-shaped candidate graph --

    root r
    m "bucket" nodes u_0..u_{m-1}, each with an edge to r
    k "item" nodes v_0..v_{k-1}, each with an edge to EVERY u_j
    each v_i has (a_i - 1) private leaf nodes, each with an edge to v_i only

-- whose minimum cost equals m*(1+B)^2 + const iff the a_i admit a perfect
3-partition into m groups of 3 summing to B each.

## Why the QAOA circuit (partition_mixer.py) only needs k*m qubits

An exchange argument shows that with this uniform, unweighted cost, EVERY
root-u edge and EVERY v-leaf edge is forced into every optimal tree: a
u_j left disconnected from r would have to route through some v_i,
strictly worsening that v_i's own successor count more than it could
save, and leaves have degree 1 so their edge is trivially forced
regardless. The only real freedom is: each v_i picks exactly ONE u_j to
attach to -- a size-k, size-m one-hot choice per item, independent across
items. This is a PARTITION matroid (exactly one edge selected per item,
out of a bipartite candidate set of size m), not the general graphic
matroid this repo's core `mixer.py` handles -- which is exactly why it
needs its own, much simpler, construction (see partition_mixer.py's
module docstring): no witness search, no tree enumeration, because
swapping which bucket one item is assigned to never depends on any other
item's state.

`gadget_graph.py` builds the FULL graph above for classical (CP-SAT)
verification only -- it is not used by the QAOA circuit.

## Objective in terms of the one-hot qubits

Writing x_{i,j} in {0,1} (1 iff item i is assigned to bucket j), the true
minimum-cost value is m*(1+B)^2+const iff a perfect partition exists, and
in general

    cost(x) = m + 2*sum(a_i) + sum(a_i^2) + 2*sum_{i<i'} a_i*a_i' * S(i,i')

where S(i,i') = sum_j x_{i,j}*x_{i',j} (1 iff i and i' land in the same
bucket, 0 otherwise -- correct precisely because x is one-hot). Only the
last (pairwise) term depends on the assignment; the constants can be
dropped from any QAOA phase circuit (they contribute only an unobservable
global phase).
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PartitionGadget:
    items: tuple[int, ...]  # a_1..a_k, the 3-PARTITION instance
    m: int  # number of buckets
    B: int  # target per-bucket sum

    def __post_init__(self) -> None:
        k = len(self.items)
        if k != 3 * self.m:
            raise ValueError(f"need k=3m items, got k={k} m={self.m}")
        if sum(self.items) != self.m * self.B:
            raise ValueError(f"sum(items)={sum(self.items)} != m*B={self.m * self.B}")
        for a in self.items:
            if not (self.B / 4 < a < self.B / 2):
                raise ValueError(f"item {a} outside required range ({self.B / 4}, {self.B / 2})")

    @property
    def k(self) -> int:
        return len(self.items)

    @property
    def n_qubits(self) -> int:
        return self.k * self.m

    def qubit(self, item: int, bucket: int) -> int:
        return item * self.m + bucket

    def linear_coeff(self, item: int) -> float:
        """Coefficient of the per-item bias term in the cost Hamiltonian
        (see module docstring's Ising expansion): applied identically to
        every one-hot qubit of this item's register."""
        a_i = self.items[item]
        s_total = sum(self.items)
        return -a_i * (s_total - a_i) / 2.0

    def quad_coeff(self, item_a: int, item_b: int) -> float:
        """Coefficient of Z_{item_a,j}*Z_{item_b,j} for every matching
        bucket index j (same value for every j, by symmetry of the
        derivation)."""
        return self.items[item_a] * self.items[item_b] / 2.0

    def diag_hamiltonian_value(self, assignment: tuple[int, ...]) -> float:
        """Diagonal value (up to an assignment-independent additive
        constant) of the RZ/RZZ cost Hamiltonian actually built into the
        circuit (partition_mixer.py), evaluated classically -- no quantum
        involved. Used to verify `linear_coeff`/`quad_coeff` reproduce the
        TRUE `bucket_cost` (up to a constant) for every assignment, i.e.
        that the cost oracle the circuit implements is the one this
        gadget claims to implement, before trusting any sampled output."""
        k, m = self.k, self.m

        def z(i: int, j: int) -> int:
            return -1 if assignment[i] == j else 1

        total = 0.0
        for i in range(k):
            coeff = self.linear_coeff(i)
            for j in range(m):
                total += coeff * z(i, j)
        for i in range(k):
            for i2 in range(i + 1, k):
                coeff = self.quad_coeff(i, i2)
                for j in range(m):
                    total += coeff * z(i, j) * z(i2, j)
        return total

    def bucket_cost(self, assignment: tuple[int, ...]) -> int:
        """Ground-truth cost for one assignment (assignment[i] = which
        bucket item i is in) via the ORIGINAL successors-squared formula,
        not the Ising expansion -- an independent reference used to check
        the Ising coefficients above are correct, not just internally
        consistent."""
        bucket_sum = [0] * self.m
        for i, b in enumerate(assignment):
            bucket_sum[b] += self.items[i]
        cost = sum((1 + s) ** 2 for s in bucket_sum)
        cost += sum(a * a for a in self.items)
        cost += sum(a - 1 for a in self.items)
        return cost

    def exact_optimum(self, stats: dict | None = None) -> int:
        """Exact DP over sorted bucket-sum states (buckets with equal
        current sum are interchangeable) -- exponentially cheaper than
        raw m^k enumeration, exact, not a heuristic.

        Its cost is the number of DISTINCT partial-sum states it meets,
        which is what makes B matter: at B=15 with items from {4,5,6,7}
        the sums collide into a few hundred states at any m; random
        large-B items don't collide and the count approaches m^k/m!.
        Pass a dict as `stats` to have `stats["states"]` count cache
        misses as they happen -- it keeps its value even if the call is
        interrupted by a timeout, which is exactly when the count is
        most informative."""
        import functools

        items_sorted = sorted(self.items, reverse=True)
        m, n = self.m, len(items_sorted)

        @functools.lru_cache(maxsize=None)
        def rec(idx: int, sums: tuple[int, ...]) -> int:
            if stats is not None:
                stats["states"] = stats.get("states", 0) + 1
            if idx == n:
                return sum((1 + s) ** 2 for s in sums)
            item = items_sorted[idx]
            best = None
            seen: set[int] = set()
            for j in range(m):
                val = sums[j] + item
                if val in seen:
                    continue
                seen.add(val)
                new_sums = tuple(sorted(sums[:j] + (val,) + sums[j + 1 :]))
                cand = rec(idx + 1, new_sums)
                if best is None or cand < best:
                    best = cand
            return best

        bucket_min = rec(0, tuple([0] * m))
        rec.cache_clear()
        return bucket_min + sum(a * a for a in items_sorted) + sum(a - 1 for a in items_sorted)
