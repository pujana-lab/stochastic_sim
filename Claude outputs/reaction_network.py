"""Stoichiometry of the Gillespie model, read from the clone classes.

Four pieces, in the order a leap uses them:

1. type_delta(clone, kind)    what ONE event does to the type totals N_<type>  (= s_j)
2. ReactionNetwork(config)    the cell types, crowding pools and observables, read once
3. network.to_arrays(...)     one row per event of the RateMatrix: a_j, s_j, n_i, ...
4. apply_counts(sim, ...)     fire event j R_j times on the real clones

Notation as in docs/PLAN_tau_leaping.md: a_j propensity, R_j firing count,
s_j state-change vector over the type totals.
"""
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from src.gillespie.clone_factory import CloneFactory
from src.gillespie.event import Event
from src.gillespie.event_type import EventType

MUTANT_SEED_SIZE = 2  # cells in a new mutant clone (same N=2 as _introduce_mutation)


# ── 1. What one event does ──────────────────────────────────────────────────

def type_delta(clone, kind: EventType) -> dict[str, int]:
    """Change of each type total when `clone` fires `kind` once.

    Mirrors TumorSimulation._apply_event. Every event except a birth removes
    one cell of the clone itself.
    """
    T = clone.get_type()
    if kind == EventType.BIRTH:
        return {T: +1}
    if kind == EventType.DEATH:
        return {T: -1}
    if kind == EventType.MUTATION:
        if not clone.next_mutation:
            raise ValueError(f"{T} mutates but has no next_mutation")
        return _add({T: -1}, {clone.next_mutation: +MUTANT_SEED_SIZE})
    if kind == EventType.EXHAUSTION:
        if not clone.converts_to:
            raise ValueError(f"{T} converts but has no converts_to")
        return _add({T: -1}, {clone.converts_to: +1})
    raise ValueError(f"No stoichiometry for {kind}")


def _add(a: dict, b: dict) -> dict:
    """Sum two deltas, so a type that is both source and target nets out."""
    out = dict(a)
    for t, v in b.items():
        out[t] = out.get(t, 0) + v
    return out


# ── 2. The network: types, crowding pools, observables ──────────────────────

class ReactionNetwork:
    def __init__(self, config):
        self.types = list(config.cell_parameters)          # column order of S
        self.column = {t: i for i, t in enumerate(self.types)}

        # Y_X: ask each type's crowding_numerator() with one cell of each type
        factory = CloneFactory(config)
        self.pools = {}
        for t in self.types:
            clone = factory.create_clone(clone_type=t, N=0)
            self.pools[t] = [u for u in self.types
                             if clone.crowding_numerator(SimpleNamespace(pop_map={u: 1})) > 0]

        # Observables for choosing tau: every N_<type>, plus every pool of 2+ types
        rows, self.observable_names = [], []
        for t in self.types:
            rows.append([1 if u == t else 0 for u in self.types])
            self.observable_names.append(f"N_{t}")
        for t, pool in self.pools.items():
            if len(pool) > 1:
                rows.append([1 if u in pool else 0 for u in self.types])
                self.observable_names.append(f"Y_{t}")
        self.L = np.array(rows)

    # ── 3. One row per event ────────────────────────────────────────────────

    def to_arrays(self, rate_matrix, tissue_state) -> "LeapArrays":
        events = rate_matrix.events
        clone_ids = list(dict.fromkeys(e.clone_id for e in events))  # unique, in order
        position = {cid: i for i, cid in enumerate(clone_ids)}

        S = np.zeros((len(events), len(self.types)), dtype=int)
        for j, e in enumerate(events):
            for t, v in type_delta(tissue_state.clones[e.clone_id], e.kind).items():
                S[j, self.column[t]] = v

        return LeapArrays(
            events=events,
            rates=np.array([e.rate for e in events]),
            S=S,
            consumes=np.array([e.kind != EventType.BIRTH for e in events]),
            clone=np.array([position[e.clone_id] for e in events]),
            clone_ids=clone_ids,
            n=np.array([tissue_state.clones[cid].N for cid in clone_ids]),
        )

    def drift_variance(self, arrays: "LeapArrays"):
        """CGP drift d_k and variance v_k of each observable (rows of L)."""
        G = arrays.S @ self.L.T                  # G[j, k] = how much event j moves observable k
        return G.T @ arrays.rates, (G ** 2).T @ arrays.rates

    def show(self, arrays: "LeapArrays") -> None:
        print(f"{'event':<28}" + "".join(f"{t:>11}" for t in self.types))
        for e, s in zip(arrays.events, arrays.S):
            print(f"{e.kind.name + ' ' + str(e.clone_id):<28}" + "".join(f"{v:>11}" for v in s))


@dataclass
class LeapArrays:
    """Channel j = events[j]. Clone i = clone_ids[i]."""
    events: list[Event]
    rates: np.ndarray      # a_j
    S: np.ndarray          # s_j as rows; columns are network.types
    consumes: np.ndarray   # True if event j removes a cell of its own clone
    clone: np.ndarray      # i of the clone that fires event j
    clone_ids: list
    n: np.ndarray          # n_i, cells per clone at time t


# ── 4. Fire the counts on the real clones ───────────────────────────────────

def apply_counts(sim, arrays: LeapArrays, R) -> None:
    """Fire event j R[j] times.

    All R_j were sampled from the state at t and every change is an
    increment, so the order of the loop does not matter.
    """
    clones = sim.tissue_state.clones
    root_id = {t: cfg.default_id for t, cfg in sim.config.cell_parameters.items()}

    for event, k in zip(arrays.events, R):
        k = int(k)
        if k == 0:
            continue
        clone = clones[event.clone_id]
        if event.kind == EventType.BIRTH:
            clone.N += k
        elif event.kind == EventType.DEATH:
            clone.N -= k
        elif event.kind == EventType.MUTATION:
            clone.N -= k
            for _ in range(k):  # k mutations = k new clones of MUTANT_SEED_SIZE, not one of 2k
                sim.create_clone(clone_id=clone.next_child_id(), clone_type=clone.next_mutation,
                                 N=MUTANT_SEED_SIZE, parent=clone.clone_id)
        elif event.kind == EventType.EXHAUSTION:
            clone.N -= k
            clones[root_id[clone.converts_to]].N += k

    if any(c.N < 0 for c in clones.values()):
        raise ValueError("A clone went negative: exits must be capped at n_i")
    sim.tissue_state.update_pop_map()
    sim._rate_matrix = None
