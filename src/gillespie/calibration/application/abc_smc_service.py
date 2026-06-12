from __future__ import annotations

import concurrent.futures
import logging
import time
from typing import Dict, List, Optional, Tuple

import numpy as np

from src.gillespie.calibration.domain.abc_result import AbcResult
from src.gillespie.calibration.domain.calibration_config import (
    CalibrationConfig,
    Prior,
    PriorDict,
)
from src.gillespie.calibration.domain.distance_metric import weighted_sse
from src.gillespie.calibration.domain.summary_statistic import (
    DEFAULT_STAT_NAMES,
    compute_summary_stats,
)
from src.gillespie.calibration.ports.calibration_ports import (
    OutputPort,
    ReferenceData,
    ReferenceDataPort,
    SimulatorPort,
)

logger = logging.getLogger(__name__)


class AbcSmcService:
    def __init__(
        self,
        config: CalibrationConfig,
        priors: PriorDict,
        simulator: SimulatorPort,
        reference: ReferenceDataPort,
        output: OutputPort,
        rng: Optional[np.random.Generator] = None,
    ):
        self.config = config
        self.priors = priors
        self.simulator = simulator
        self.output = output
        self.rng = rng if rng is not None else np.random.default_rng(config.seed)

        ref_data = reference.load()
        self.ref_stats: np.ndarray = ref_data.ref_stats
        self.precision_weights: np.ndarray = ref_data.precision_weights
        self.stat_names: List[str] = ref_data.stat_names
        self.cell_types: List[str] = ref_data.cell_types

        self.param_names: List[str] = list(priors.keys())

    # ── Proposal helpers ──────────────────────────────────────────────────────

    def _sample_prior(self) -> Optional[Dict[str, float]]:
        theta: Dict[str, float] = {}
        for pname in self.param_names:
            prior = self.priors[pname]
            theta[pname] = prior.sample(self.rng)

        for pname in self.param_names:
            prior = self.priors[pname]
            if theta[pname] < prior.lo or theta[pname] > prior.hi:
                return None

        for pname in self.param_names:
            prior = self.priors[pname]
            if prior.constraint is not None and not prior.constraint(theta):
                return None

        return theta

    def _sample_from_previous_gen(
        self,
        prev_particles: List[Dict[str, float]],
        prev_weights: np.ndarray,
    ) -> Optional[Dict[str, float]]:
        idx = self.rng.choice(len(prev_particles), p=prev_weights)
        parent = prev_particles[idx]

        sigma = self._compute_perturbation_sigma(prev_particles, prev_weights)
        theta: Dict[str, float] = {}
        for pname in self.param_names:
            theta[pname] = float(self.rng.normal(parent[pname], sigma[pname]))

        for pname in self.param_names:
            prior = self.priors[pname]
            if theta[pname] < prior.lo or theta[pname] > prior.hi:
                return None

        for pname in self.param_names:
            prior = self.priors[pname]
            if prior.constraint is not None and not prior.constraint(theta):
                return None

        return theta

    def _compute_perturbation_sigma(
        self,
        particles: List[Dict[str, float]],
        weights: np.ndarray,
    ) -> Dict[str, float]:
        sigma: Dict[str, float] = {}
        for pname in self.param_names:
            vals = np.array([p[pname] for p in particles])
            mean = np.average(vals, weights=weights)
            var = np.average((vals - mean) ** 2, weights=weights)
            sigma[pname] = 2.0 * float(np.sqrt(max(var, 1e-12)))
        return sigma

    # ── Importance weight computation ────────────────────────────────────────

    def _compute_weights(
        self,
        thetas: List[Dict[str, float]],
        prev_particles: List[Dict[str, float]],
        prev_weights: np.ndarray,
    ) -> np.ndarray:
        n_new = len(thetas)
        n_prev = len(prev_particles)
        sigma = self._compute_perturbation_sigma(prev_particles, prev_weights)

        log_weights = np.zeros(n_new)
        for i, theta in enumerate(thetas):
            log_prior_val = 0.0
            for pname in self.param_names:
                prior = self.priors[pname]
                lp = np.log(max(prior.pdf(theta[pname]), 1e-300))
                log_prior_val += lp

            log_denom_parts = np.zeros(n_prev)
            for j in range(n_prev):
                log_k = 0.0
                for pname in self.param_names:
                    s = sigma[pname]
                    if s > 0:
                        diff = (theta[pname] - prev_particles[j][pname]) / s
                        log_k += -0.5 * diff * diff - np.log(s) - 0.5 * np.log(2 * np.pi)
                    else:
                        log_k += 0.0 if diff == 0 else -np.inf
                log_denom_parts[j] = np.log(max(prev_weights[j], 1e-300)) + log_k

            log_denom = np.logaddexp.reduce(log_denom_parts)

            if np.isfinite(log_prior_val) and np.isfinite(log_denom):
                log_weights[i] = log_prior_val - log_denom
            else:
                log_weights[i] = -np.inf

        max_log_w = np.max(log_weights)
        if not np.isfinite(max_log_w):
            return np.ones(n_new) / n_new

        weights = np.exp(log_weights - max_log_w)
        w_sum = np.sum(weights)
        if w_sum > 0:
            weights /= w_sum
        else:
            weights = np.ones(n_new) / n_new

        return weights

    # ── Worker simulation ────────────────────────────────────────────────────

    def _worker(self, theta: Dict[str, float]) -> Optional[float]:
        try:
            dists = []
            base_seed = int(self.rng.integers(1, 2_000_000_000))
            for rep in range(self.config.n_reps):
                params = dict(theta)
                params["seed"] = base_seed + rep
                times, pop = self.simulator.run(params)
                sim_stats = compute_summary_stats(
                    pop, times, self.stat_names, self.cell_types
                )
                dist = weighted_sse(sim_stats, self.ref_stats, self.precision_weights)
                dists.append(dist)
            return float(np.mean(dists))
        except Exception as e:
            logger.debug("Worker simulation failed: %s", e)
            return None

    # ── Run a single generation ──────────────────────────────────────────────

    def run_generation(
        self,
        gen: int,
        epsilon: float,
        prev_particles: Optional[List[Dict[str, float]]] = None,
        prev_weights: Optional[np.ndarray] = None,
    ) -> Tuple[List[Dict[str, float]], np.ndarray, np.ndarray, Dict]:
        accepted_thetas: List[Dict[str, float]] = []
        accepted_dists: List[float] = []
        attempts = 0
        prior_rejects = 0
        sim_fails = 0
        start = time.time()

        with concurrent.futures.ThreadPoolExecutor(
            max_workers=self.config.n_workers
        ) as pool:
            while len(accepted_thetas) < self.config.n_particles:
                needed = self.config.n_particles - len(accepted_thetas)
                batch_size = min(needed * 3, self.config.max_attempts - attempts)

                proposals: List[Dict[str, float]] = []
                for _ in range(batch_size):
                    if gen == 0:
                        theta = self._sample_prior()
                    else:
                        theta = self._sample_from_previous_gen(
                            prev_particles, prev_weights
                        )
                    attempts += 1
                    if theta is None:
                        prior_rejects += 1
                        continue
                    proposals.append(theta)

                if not proposals:
                    if attempts >= self.config.max_attempts:
                        break
                    continue

                future_map = {
                    pool.submit(self._worker, t): t for t in proposals
                }
                for fut in concurrent.futures.as_completed(future_map):
                    try:
                        dist = fut.result()
                    except Exception:
                        sim_fails += 1
                        continue

                    if dist is not None and dist <= epsilon:
                        accepted_thetas.append(future_map[fut])
                        accepted_dists.append(dist)
                        if len(accepted_thetas) >= self.config.n_particles:
                            break

                if attempts >= self.config.max_attempts:
                    break

        if len(accepted_thetas) < self.config.n_particles:
            raise RuntimeError(
                f"Generation {gen}: only accepted {len(accepted_thetas)}/"
                f"{self.config.n_particles} particles after {attempts} attempts "
                f"(prior_rejects={prior_rejects}, sim_fails={sim_fails}). "
                f"Try wider priors, higher epsilon_0, or more max_attempts."
            )

        accepted_arr = np.array(accepted_dists)
        if gen == 0:
            weights = np.ones(len(accepted_thetas)) / len(accepted_thetas)
        else:
            weights = self._compute_weights(
                accepted_thetas, prev_particles, prev_weights
            )

        ess = 1.0 / float(np.sum(weights ** 2))

        diagnostics = {
            "generation": gen,
            "epsilon": epsilon,
            "attempts": attempts,
            "accepted": len(accepted_thetas),
            "prior_reject": prior_rejects,
            "sim_fail": sim_fails,
            "acceptance_rate": len(accepted_thetas) / max(attempts, 1),
            "ess": ess,
            "ess_ratio": ess / len(accepted_thetas),
            "mean_distance": float(np.mean(accepted_arr)),
            "min_distance": float(np.min(accepted_arr)),
            "elapsed_s": time.time() - start,
        }

        logger.info(
            "Gen %d | eps=%.4f | acc=%d/%d (%.1f%%) | ESS=%.1f | prior_rej=%d | sim_fail=%d | %.1fs",
            gen,
            epsilon,
            len(accepted_thetas),
            attempts,
            100.0 * diagnostics["acceptance_rate"],
            ess,
            prior_rejects,
            sim_fails,
            diagnostics["elapsed_s"],
        )

        return accepted_thetas, weights, accepted_arr, diagnostics

    # ── Full ABC-SMC loop ────────────────────────────────────────────────────

    def run_abc_smc(self) -> AbcResult:
        self.output.save_manifest(
            calibration_config={
                f.name: getattr(self.config, f.name)
                for f in self.config.__class__.__dataclass_fields__.values()
            },
            priors=self.priors,
        )

        result = AbcResult()
        epsilon = self.config.epsilon_0

        for gen in range(self.config.n_generations):
            if gen == 0:
                particles, weights, dists, diag = self.run_generation(
                    gen=0, epsilon=epsilon
                )
            else:
                prev_particles = result.particles[-1]
                prev_weights = result.weights[-1]
                dists_prev = result.distances[-1]
                epsilon = float(np.quantile(dists_prev, self.config.alpha))
                epsilon = max(epsilon, self.config.epsilon_final)
                particles, weights, dists, diag = self.run_generation(
                    gen=gen,
                    epsilon=epsilon,
                    prev_particles=prev_particles,
                    prev_weights=prev_weights,
                )

            result.particles.append(particles)
            result.weights.append(weights)
            result.distances.append(dists)
            result.diagnostics.append(diag)
            result.epsilon_history.append(epsilon)

            self.output.save_generation(gen, particles, weights, dists)
            self.output.save_generation_diagnostics(result.diagnostics)
            self.output.save_epsilon_schedule(result.epsilon_history)

        return result

    # ── Resume ───────────────────────────────────────────────────────────────

    def resume(self) -> AbcResult:
        last_gen = self.output.find_last_generation()
        if last_gen is None:
            logger.info("No saved generations found. Starting fresh.")
            return self.run_abc_smc()

        logger.info("Resuming from generation %d", last_gen)
        result = AbcResult()

        for g in range(last_gen + 1):
            particles, weights, dists = self.output.load_generation(g)
            result.particles.append(particles)
            result.weights.append(weights)
            result.distances.append(dists)

        diagnostics_df = None
        eps_sched_df = None
        try:
            import pandas as pd

            diag_path = (
                f"{self.config.output_dir}/generation_diagnostics.csv"
            )
            eps_path = f"{self.config.output_dir}/epsilon_schedule.csv"

            diagnostics_df = pd.read_csv(diag_path)
            eps_sched_df = pd.read_csv(eps_path)
        except Exception:
            pass

        if diagnostics_df is not None:
            result.diagnostics = diagnostics_df.to_dict(orient="records")
        if eps_sched_df is not None:
            result.epsilon_history = eps_sched_df["epsilon"].tolist()

        for g in range(last_gen + 1, self.config.n_generations):
            prev_particles = result.particles[-1]
            prev_weights = result.weights[-1]
            dists_prev = result.distances[-1]
            epsilon = float(np.quantile(dists_prev, self.config.alpha))
            epsilon = max(epsilon, self.config.epsilon_final)

            particles, weights, dists, diag = self.run_generation(
                gen=g,
                epsilon=epsilon,
                prev_particles=prev_particles,
                prev_weights=prev_weights,
            )

            result.particles.append(particles)
            result.weights.append(weights)
            result.distances.append(dists)
            result.diagnostics.append(diag)
            result.epsilon_history.append(epsilon)

            self.output.save_generation(g, particles, weights, dists)
            self.output.save_generation_diagnostics(result.diagnostics)
            self.output.save_epsilon_schedule(result.epsilon_history)

        return result

    # ── Prior sensitivity analysis ──────────────────────────────────────────

    def prior_sensitivity(
        self,
        width_factors: Optional[List[float]] = None,
        n_samples: int = 1000,
    ) -> Dict[float, Dict]:
        if width_factors is None:
            width_factors = [0.5, 0.8, 1.0, 1.2, 1.5]

        results: Dict[float, Dict] = {}
        eps_ref = self.config.epsilon_0

        for factor in width_factors:
            thetas: List[Dict[str, float]] = []
            prior_rejects = 0
            for _ in range(n_samples):
                theta: Dict[str, float] = {}
                valid = True
                for pname in self.param_names:
                    prior = self.priors[pname]
                    mid = (prior.lo + prior.hi) / 2
                    half = (prior.hi - prior.lo) / 2 * factor
                    lo = mid - half
                    hi = mid + half
                    val = self.rng.uniform(lo, hi)
                    if val < prior.lo or val > prior.hi:
                        valid = False
                        break
                    theta[pname] = val
                if valid:
                    thetas.append(theta)
                else:
                    prior_rejects += 1

            dists: List[float] = []
            sim_fails = 0
            for theta in thetas[:n_samples]:
                try:
                    times, pop = self.simulator.run(theta)
                    sim_stats = compute_summary_stats(
                        pop, times, self.stat_names, self.cell_types
                    )
                    dist = weighted_sse(
                        sim_stats, self.ref_stats, self.precision_weights
                    )
                    dists.append(dist)
                except Exception:
                    sim_fails += 1

            n_accepted = sum(1 for d in dists if d <= eps_ref) if dists else 0

            results[factor] = {
                "n_samples": len(dists),
                "n_accepted": n_accepted,
                "acceptance_rate": n_accepted / max(len(dists), 1),
                "mean_distance": float(np.mean(dists)) if dists else None,
                "prior_rejects": prior_rejects,
                "sim_fails": sim_fails,
            }

        return results
