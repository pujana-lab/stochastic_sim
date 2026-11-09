"""How the simulation advances: one interchangeable step method.

`TumorSimulation` owns the state (clones, rng, recorded history) and the shared
helpers (`_current_rate_matrix`, `_build_rate_matrix`, `_apply_event`,
`_advance_all_instability`, `_record_rates`, `_record`). A `Stepper` decides how long one step is (τ) and how much happens in it:
`SSAStepper` fires one exact Gillespie event; the tau-leaping stepper of
`docs/PLAN_tau_leaping.md` (Phase 3) will fire many at once.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from src.gillespie.tumor_simulation import TumorSimulation


class Stepper(ABC):
    """One step of the simulation. Returns False when the simulation should stop."""

    @abstractmethod
    def step(self, sim: "TumorSimulation") -> bool:
        ...


class SSAStepper(Stepper):
    """Exact Gillespie SSA: one event per step, rates rebuilt after every event.

    Tags in the comments of step():
      [DYNAMICS]   changes the state or what comes next; leaving it out changes the trajectory
      [REPORTING]  only records for the analysis; leaving it out gives the same trajectory
      [ONLY IF …]  runs only under that condition
    The order of the two rng draws (τ, then the event) fixes seeded trajectories.
    """

    def step(self, sim: "TumorSimulation") -> bool:
        # ── 1. Read the current state ────────────────────────────────────────
        # [REPORTING] step counter: only used for "a snapshot every save_interval steps"
        sim.step_count += 1
        # [DYNAMICS] rates of the current state, built at the end of the previous step
        rate_matrix = sim._current_rate_matrix()
        total_rate = rate_matrix.get_total_rate()                 # a₀
        # [REPORTING] [ONLY IF save_interval > 0] rate history (rate_his.parquet)
        sim._record_rates(rate_matrix)

        # ── 2. Stop: nothing can happen ──────────────────────────────────────
        # [DYNAMICS] [ONLY IF a₀ = 0] e.g. every population extinct or at capacity
        if total_rate <= 0 or not rate_matrix.events:
            return False

        # ── 3. When does the next event happen ───────────────────────────────
        # [DYNAMICS] τ ~ Exp(a₀): first rng draw of the step
        tau = self._sample_waiting_time(sim, total_rate)
        new_t = sim.t + tau

        # ── 4. Stop: the next event falls after T_max ────────────────────────
        # [ONLY IF t + τ > T_max] advance to T_max without any event, record the final state
        if new_t > sim.config.T_max:
            sim._advance_all_instability(sim.config.T_max - sim.t)   # [DYNAMICS]
            sim.t = sim.config.T_max                                  # [DYNAMICS]
            sim._record(final=True)   # [DYNAMICS] + [REPORTING], see TumorSimulation._record
            return False

        # ── 5. Advance: move time, fire one event ────────────────────────────
        sim._advance_all_instability(tau)                    # [DYNAMICS] instability grows over τ
        sim.t = new_t                                        # [DYNAMICS]
        event = rate_matrix.choose_event(sim.rng.random())   # [DYNAMICS] second rng draw: which event
        sim._apply_event(event)                              # [DYNAMICS] birth / death / mutation / conversion

        # ── 6. Rebuild the rates and record ──────────────────────────────────
        sim._record(event)            # [DYNAMICS] + [REPORTING], see TumorSimulation._record
        return True

    @staticmethod
    def _sample_waiting_time(sim: "TumorSimulation", total_rate: float) -> float:
        """τ ~ Exp(a₀), the time to the next event, drawn from the simulation's rng.

        Choosing τ is what differs between steppers: the tau-leaping stepper
        computes it (select_tau, PLAN_tau_leaping.md §5) instead of drawing it.
        """
        if total_rate <= 0:
            raise ValueError("Total rate must be positive.")
        return -np.log(sim.rng.random()) / total_rate

    def run_burst(self, sim: "TumorSimulation", n: int) -> bool:
        """Up to n exact steps; the tau-leaping fallback when a leap would buy too few events.

        Returns False as soon as a step says the simulation should stop (no
        events left or T_max reached), True if all n steps ran.
        """
        for _ in range(n):
            if not self.step(sim):
                return False
        return True


# Steppers by the name used in SimulationConfig.stepper, the JSON and --stepper.
STEPPERS: dict[str, type[Stepper]] = {
    "ssa": SSAStepper,
}


def make_stepper(name: str) -> Stepper:
    """The stepper called `name` (see STEPPERS)."""
    try:
        return STEPPERS[name]()
    except KeyError:
        raise ValueError(f"Unknown stepper {name!r}; known: {sorted(STEPPERS)}") from None
