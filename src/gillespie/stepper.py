"""How the simulation advances: one interchangeable step method.

`TumorSimulation` owns the state (clones, rng, recorded history) and the shared
helpers (`_current_rate_matrix`, `_build_rate_matrix`, `_apply_event`,
`_advance_all_instability`). A `Stepper` decides how much happens in one step:
`SSAStepper` fires one exact Gillespie event; the tau-leaping stepper of
`docs/PLAN_tau_leaping.md` (Phase 3) will fire many at once.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from src.gillespie.tumor_simulation import TumorSimulation


class Stepper(ABC):
    """One step of the simulation. Returns False when the simulation should stop."""

    @abstractmethod
    def step(self, sim: "TumorSimulation") -> bool:
        ...


class SSAStepper(Stepper):
    """Exact Gillespie SSA: one event per step, rates rebuilt after every event."""

    def step(self, sim: "TumorSimulation") -> bool:
        sim.step_count += 1
        rate_matrix = sim._current_rate_matrix()
        total_rate = rate_matrix.get_total_rate()

        # Only record rates if saving history
        if sim.save_interval > 0:
            step_rates = [
                {
                    "time": sim.t,
                    "kind": e.kind.name,
                    "clone_id": e.clone_id,
                    "clone_type": e.clone_type,
                    "rate": e.rate,
                } for e in rate_matrix.events
            ]
            sim.rate_history.append(step_rates)

        if total_rate <= 0 or not rate_matrix.events:
            return False

        tau = sim._sample_waiting_time(total_rate)
        new_t = sim.t + tau

        if new_t > sim.config.T_max:
            tau = sim.config.T_max - sim.t
            sim._advance_all_instability(tau)
            sim.t = sim.config.T_max
            sim.times.append(sim.t)
            sim.tissue_state.t = sim.t
            sim._rate_matrix = sim._build_rate_matrix()
            # Always save final state
            sim.history.append(sim.tissue_state.snapshot())
            print("FINAL RATES AND STATE")
            sim.tissue_state.print_pop_map()
            return False

        sim._advance_all_instability(tau)
        sim.t = new_t

        event = rate_matrix.choose_event(sim.rng.random())
        sim._apply_event(event)

        sim.times.append(sim.t)
        sim.tissue_state.t = sim.t
        # Rates of the post-event state: used by the snapshot below and by the next step
        sim._rate_matrix = sim._build_rate_matrix()

        # Save snapshots conditionally
        if sim.save_interval > 0 and sim.step_count % sim.save_interval == 0:
            sim.history.append(sim.tissue_state.snapshot())

        # Save events if recording all steps
        if sim.save_all_steps:
            sim.events.append(event)

        return True

    def run_burst(self, sim: "TumorSimulation", n: int) -> bool:
        """Up to n exact steps; the tau-leaping fallback when a leap would buy too few events.

        Returns False as soon as a step says the simulation should stop (no
        events left or T_max reached), True if all n steps ran.
        """
        for _ in range(n):
            if not self.step(sim):
                return False
        return True
