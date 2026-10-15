"""Conversions between dimensional, dimensionless, and Gillespie inputs.

The formulas follow sections 2.1-2.2 and 4.3 of ``main.pdf``.  The
Gillespie representation uses absolute cell counts and microscopic
interaction propensities, hence bilinear rates are divided by ``OMEGA``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


def _positive(value: float, name: str) -> float:
    if value <= 0:
        raise ValueError(f"{name} must be positive")
    return value


@dataclass(frozen=True)
class PhysicalParameters:
    lambda_C: float
    lambda_I: float
    mu_C: float
    mu_E: float
    theta: float
    beta: float
    omega: float
    K_C: float
    K_I: float
    C0: float = 0.0
    I0: float = 0.0
    E0: float = 0.0


@dataclass(frozen=True)
class DimensionlessParameters:
    gamma_C: float
    sigma: float
    rho: float
    beta: float
    omega: float
    mu_E: float
    C0: float = 0.0
    I0: float = 0.0
    E0: float = 0.0


@dataclass(frozen=True)
class GillespieParameters:
    omega_system: int
    cell_parameters: dict[str, dict[str, Any]]
    theta_I: float
    beta: float
    source: str = "physical"

    def as_json(self) -> dict[str, Any]:
        return {
            "OMEGA": self.omega_system,
            "cell_parameters": self.cell_parameters,
            "theta_I": self.theta_I,
            "beta": self.beta,
            "scale": True,
        }


def dimensional_to_dimensionless(
    physical: PhysicalParameters,
) -> DimensionlessParameters:
    """Apply equation (4) from the report."""
    lambda_c = _positive(physical.lambda_C, "lambda_C")
    K_i = _positive(physical.K_I, "K_I")
    return DimensionlessParameters(
        gamma_C=physical.mu_C / lambda_c,
        sigma=physical.theta * K_i / lambda_c,
        rho=physical.lambda_I / lambda_c,
        beta=physical.beta * physical.K_C / lambda_c,
        omega=physical.omega * physical.K_C / lambda_c,
        mu_E=physical.mu_E / lambda_c,
        C0=physical.C0 / _positive(physical.K_C, "K_C"),
        I0=physical.I0 / K_i,
        E0=physical.E0 / K_i,
    )


def dimensionless_to_dimensional(
    dimensionless: DimensionlessParameters,
    *,
    lambda_C: float,
    K_C: float,
    K_I: float,
) -> PhysicalParameters:
    """Invert equation (4), using the supplied physical reference scales."""
    lambda_c = _positive(lambda_C, "lambda_C")
    k_c = _positive(K_C, "K_C")
    k_i = _positive(K_I, "K_I")
    return PhysicalParameters(
        lambda_C=lambda_c,
        lambda_I=dimensionless.rho * lambda_c,
        mu_C=dimensionless.gamma_C * lambda_c,
        mu_E=dimensionless.mu_E * lambda_c,
        theta=dimensionless.sigma * lambda_c / k_i,
        beta=dimensionless.beta * lambda_c / k_c,
        omega=dimensionless.omega * lambda_c / k_c,
        K_C=k_c,
        K_I=k_i,
        C0=dimensionless.C0 * k_c,
        I0=dimensionless.I0 * k_i,
        E0=dimensionless.E0 * k_i,
    )


def to_gillespie(
    physical: PhysicalParameters,
    *,
    omega_system: int,
    base0: int = 0,
    mutation_rate: float = 0.0,
    fitness_gain: float = 0.0,
) -> GillespieParameters:
    """Convert physical rates and populations to the current JSON contract.

    ``C`` is represented by the simulator's ``mutated`` population.  The
    optional ``base0`` keeps the wild-type compartment available.
    """
    if omega_system <= 0:
        raise ValueError("ommu_ega_system must be positive")
    k = omega_system
    user_k_c = physical.K_C / k
    user_k_i = physical.K_I / k
    return GillespieParameters(
        omega_system=omega_system,
        # These are configuration-level rates. SimulationConfig performs the
        # single division by OMEGA when scale=True.
        theta_I=physical.theta,
        beta=physical.beta,
        cell_parameters={
            "base": {
                "N": int(round(base0)),
                "K": user_k_c,
                "nu": mutation_rate,
                "next_mutation": "mutated",
            },
            "mutated": {
                "N": int(round(physical.C0)),
                "K": user_k_c,
                "lambda0": physical.lambda_C,
                "mu": physical.mu_C,
                "fitness_gain": fitness_gain,
            },
            "immune": {
                "N": int(round(physical.I0)),
                "K": user_k_i,
                "lambda0": physical.lambda_I,
                "omega_exhaust": physical.omega,
            },
            "exhausted": {
                "N": int(round(physical.E0)),
                "K": user_k_i,
                "mu": physical.mu_E,
            },
        },
        source="physical",
    )


def write_json(parameters: GillespieParameters, path: Path) -> None:
    """Write a calculator result that can be passed through the JSON adapter."""
    path.write_text(json.dumps(parameters.as_json(), indent=2) + "\n")


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert model parameters for Gillespie.")
    parser.add_argument("--output", type=Path, required=True, help="Output JSON path.")
    parser.add_argument("--lambda-c", type=float, required=True)
    parser.add_argument("--lambda-i", type=float, required=True)
    parser.add_argument("--mu-c", type=float, required=True)
    parser.add_argument("--mu-e", type=float, required=True)
    parser.add_argument("--theta", type=float, required=True)
    parser.add_argument("--beta", type=float, required=True)
    parser.add_argument("--omega", type=float, required=True)
    parser.add_argument("--k-c", type=float, required=True)
    parser.add_argument("--k-i", type=float, required=True)
    parser.add_argument("--c0", type=float, default=0.0)
    parser.add_argument("--i0", type=float, default=0.0)
    parser.add_argument("--e0", type=float, default=0.0)
    parser.add_argument("--base0", type=int, default=0)
    parser.add_argument("--system-size", type=int, default=4000)
    parser.add_argument("--mutation-rate", type=float, default=0.0)
    parser.add_argument("--fitness-gain", type=float, default=0.0)
    return parser


def main(args: list[str] | None = None) -> None:
    options = _parser().parse_args(args)
    physical = PhysicalParameters(
        lambda_C=options.lambda_c,
        lambda_I=options.lambda_i,
        mu_C=options.mu_c,
        mu_E=options.mu_e,
        theta=options.theta,
        beta=options.beta,
        omega=options.omega,
        K_C=options.k_c,
        K_I=options.k_i,
        C0=options.c0,
        I0=options.i0,
        E0=options.e0,
    )
    result = to_gillespie(
        physical,
        omega_system=options.system_size,
        base0=options.base0,
        mutation_rate=options.mutation_rate,
        fitness_gain=options.fitness_gain,
    )
    payload = result.as_json()
    payload["_calculator"] = {
        "physical": asdict(physical),
        "dimensionless": asdict(dimensional_to_dimensionless(physical)),
        "notes": [
            "SimulationConfig divides theta_I, beta, and omega_exhaust by OMEGA once.",
            "K is the dimensionless capacity; SimulationConfig scales it to OMEGA*K.",
            "AdaptedCrowding then inflates actual_K to balance natural death.",
        ],
    }
    options.output.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
