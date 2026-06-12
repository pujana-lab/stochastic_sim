from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, Optional

import numpy as np
from scipy.stats import truncnorm

PriorConstraint = Callable[[Dict[str, float]], bool]

PriorDict = Dict[str, "Prior"]


@dataclass(frozen=True)
class CalibrationConfig:
    n_particles: int = 200
    n_generations: int = 5
    alpha: float = 0.6
    epsilon_0: float = 1e6
    epsilon_final: float = 1.0
    max_attempts: int = 100000
    n_workers: int = 4
    n_reps: int = 1
    seed: int = 42
    output_dir: str = "calibration_results"


@dataclass(frozen=True)
class Prior:
    param_name: str
    lo: float
    hi: float
    distribution: str = "uniform"
    constraint: Optional[PriorConstraint] = None

    def __post_init__(self) -> None:
        valid = {"uniform", "lognorm", "truncnorm"}
        if self.distribution not in valid:
            raise ValueError(
                f"Invalid distribution '{self.distribution}'. Must be one of {valid}"
            )
        if self.lo >= self.hi:
            raise ValueError(
                f"Prior '{self.param_name}': lo ({self.lo}) must be < hi ({self.hi})"
            )

    def sample(self, rng: np.random.Generator) -> float:
        if self.distribution == "uniform":
            return rng.uniform(self.lo, self.hi)
        if self.distribution == "lognorm":
            log_lo = np.log(self.lo)
            log_hi = np.log(self.hi)
            return np.exp(rng.uniform(log_lo, log_hi))
        if self.distribution == "truncnorm":
            mu = (self.lo + self.hi) / 2
            sigma = (self.hi - self.lo) / 6
            a = (self.lo - mu) / sigma
            b = (self.hi - mu) / sigma
            return float(truncnorm.rvs(a, b, loc=mu, scale=sigma, random_state=rng))
        raise ValueError(f"Unknown distribution: {self.distribution}")

    def pdf(self, theta: float) -> float:
        if theta < self.lo or theta > self.hi:
            return 0.0
        if self.distribution == "uniform":
            return 1.0 / (self.hi - self.lo)
        if self.distribution == "lognorm":
            if theta <= 0:
                return 0.0
            return 1.0 / (theta * (np.log(self.hi) - np.log(self.lo)))
        if self.distribution == "truncnorm":
            mu = (self.lo + self.hi) / 2
            sigma = (self.hi - self.lo) / 6
            a = (self.lo - mu) / sigma
            b = (self.hi - mu) / sigma
            return float(truncnorm.pdf(theta, a, b, loc=mu, scale=sigma))
        raise ValueError(f"Unknown distribution: {self.distribution}")
