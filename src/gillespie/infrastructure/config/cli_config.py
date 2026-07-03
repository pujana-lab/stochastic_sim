from __future__ import annotations

import argparse
from dataclasses import fields as dataclass_fields, MISSING
from pathlib import Path

from src.gillespie.simulation_config import SimulationConfig


def _config_field_names() -> set[str]:
    return {f.name for f in dataclass_fields(SimulationConfig)}


def _get_config_defaults() -> dict:
    """Extract default values directly from SimulationConfig dataclass fields."""
    defaults = {}
    for field in dataclass_fields(SimulationConfig):
        if field.default is not MISSING:
            defaults[field.name] = field.default
        elif field.default_factory is not MISSING:
            defaults[field.name] = field.default_factory()
    return defaults


def build_parser() -> argparse.ArgumentParser:
    defaults = _get_config_defaults()
    
    p = argparse.ArgumentParser(description="Run tumor clone Gillespie simulation.")
    p.add_argument("--config", type=Path, default=None, help="JSON config file (CLI args override it).")

    p.add_argument("--N0", type=int, default=defaults.get("N0"))
    p.add_argument("--lambda0", type=float, default=defaults.get("lambda0"))
    p.add_argument("--mu0", type=float, default=defaults.get("mu0"))
    p.add_argument("--nu0", type=float, default=defaults.get("nu0"))

    p.add_argument("--d1-0", dest="d1_0", type=float, default=defaults.get("d1_0"))
    p.add_argument("--d2-0", dest="d2_0", type=float, default=defaults.get("d2_0"))

    p.add_argument("--instability-0", dest="instability_0", type=float, default=defaults.get("instability_0"))
    p.add_argument("--buildup-0", dest="buildup_0", type=float, default=defaults.get("buildup_0"))
    p.add_argument("--base-instability-buildup", type=float, default=defaults.get("base_instability_buildup"))

    p.add_argument("--mutation-instability-jump", type=float, default=defaults.get("mutation_instability_jump"))
    p.add_argument("--mutation-buildup-gain", type=float, default=defaults.get("mutation_buildup_gain"))

    p.add_argument("--N-immune", dest="N_immune", type=int, default=defaults.get("N_immune"))
    p.add_argument("--N-exhausted", dest="N_exhausted", type=int, default=defaults.get("N_exhausted"))
    p.add_argument("--N-mutant", dest="N_mutant", type=int, default=defaults.get("N_mutant"))

    p.add_argument("--lambda_0-Immune", dest="lambda_Immune", type=float, default=defaults.get("lambda_Immune"))

    p.add_argument("--mu-Immune", dest="omega_exhaust", type=float, default=defaults.get("omega_exhaust"))
    p.add_argument("--mu-Exhausted", dest="mu_Exhausted", type=float, default=defaults.get("mu_Exhausted"))

    p.add_argument("--T-max", dest="T_max", type=float, default=defaults.get("T_max"))
    p.add_argument("--seed", type=int, default=defaults.get("seed"))

    p.add_argument("--use-logistic", dest="use_logistic", action=argparse.BooleanOptionalAction, default=defaults.get("use_logistic"))
    p.add_argument("--use-logistic-adapted", dest="use_logistic_adapted", action=argparse.BooleanOptionalAction, default=defaults.get("use_logistic_adapted"))

    p.add_argument("--K0", type=float, default=defaults.get("K0"))
    p.add_argument("--K-immune", dest="K_immune", type=float, default=defaults.get("K_immune"))
    p.add_argument("--K-mutant", dest="K_mutant", type=float, default=defaults.get("K_mutant"))
    p.add_argument("--decline", type=float, default=defaults.get("decline"))
    p.add_argument("--Kmin", type=float, default=defaults.get("Kmin"))

    p.add_argument("--fitness-gain", dest="fitness_gain", type=float, default=defaults.get("fitness_gain"))
    p.add_argument("--theta-I", dest="theta_I", type=float, default=defaults.get("theta_I"))
    p.add_argument("--beta", type=float, default=defaults.get("beta"))

    p.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=defaults.get("verbose"))
    p.add_argument("--scale", action=argparse.BooleanOptionalAction, default=defaults.get("scale"))
    p.add_argument("--decay", action=argparse.BooleanOptionalAction, default=defaults.get("decay"))

    p.add_argument("--omega", dest="OMEGA", type=int, default=defaults.get("OMEGA"))

    p.add_argument("--save-history", type=Path, default=Path("./history.csv"), help="Write long-format history CSV.")
    p.add_argument("--save-clones", type=Path, default=None, help="Write final clone states CSV.")
    p.add_argument("--save-debug", type=Path, default=Path("./debug.csv"), help="Write full debug trace CSV (all fields + event).")
    p.add_argument("--top", type=int, default=10, help="How many largest final clones to print.")

    return p


def get_explicit_cli_args(namespace: argparse.Namespace, parser: argparse.ArgumentParser) -> dict:
    defaults = vars(parser.parse_args([]))
    config_fields = _config_field_names()
    result = {}
    for k, v in vars(namespace).items():
        if k in config_fields:
            if k in defaults:
                if v != defaults[k]:
                    result[k] = v
            else:
                result[k] = v
    return result


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    return parser.parse_args(args)
