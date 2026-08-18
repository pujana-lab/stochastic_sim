from __future__ import annotations

import argparse
from dataclasses import fields as dataclass_fields
from pathlib import Path

from src.gillespie.simulation_config import SimulationConfig


def _config_field_names() -> set[str]:
    return {f.name for f in dataclass_fields(SimulationConfig)}


def _config_defaults() -> dict:
    cfg = SimulationConfig()
    return {
        "T_max": cfg.T_max,
        "seed": cfg.seed,
        "decline": cfg.decline,
        "Kmin": cfg.Kmin,
        "theta_I": cfg.theta_I,
        "beta": cfg.beta,
        "d1_0": cfg.d1_0,
        "d2_0": cfg.d2_0,
        "instability_0": cfg.instability_0,
        "buildup_0": cfg.buildup_0,
        "base_instability_buildup": cfg.base_instability_buildup,
        "mutation_instability_jump": cfg.mutation_instability_jump,
        "mutation_buildup_gain": cfg.mutation_buildup_gain,
        "verbose": cfg.verbose,
        "scale": cfg.scale,
        "decay": cfg.decay,
        "use_logistic": cfg.use_logistic,
        "use_logistic_adapted": cfg.use_logistic_adapted,
        "save_all_steps": cfg.save_all_steps,
        "save_interval": cfg.save_interval,
        "OMEGA": cfg.OMEGA,
        "cell_parameters": cfg.cell_parameters,
    }


def build_parser() -> argparse.ArgumentParser:
    defaults = _config_defaults()
    p = argparse.ArgumentParser(description="Run tumor clone Gillespie simulation.")
    p.add_argument("--config", type=Path, default=None, help="JSON config file (CLI args override it).")
    p.add_argument("--T-max", dest="T_max", type=float, default=defaults["T_max"])
    p.add_argument("--seed", type=int, default=defaults["seed"])
    p.add_argument("--decline", type=float, default=defaults["decline"])
    p.add_argument("--Kmin", type=float, default=defaults["Kmin"])
    p.add_argument("--theta-I", dest="theta_I", type=float, default=defaults["theta_I"])
    p.add_argument("--beta", type=float, default=defaults["beta"])
    p.add_argument("--d1-0", dest="d1_0", type=float, default=defaults["d1_0"])
    p.add_argument("--d2-0", dest="d2_0", type=float, default=defaults["d2_0"])
    p.add_argument("--instability-0", dest="instability_0", type=float, default=defaults["instability_0"])
    p.add_argument("--buildup-0", dest="buildup_0", type=float, default=defaults["buildup_0"])
    p.add_argument("--base-instability-buildup", type=float, default=defaults["base_instability_buildup"])
    p.add_argument("--mutation-instability-jump", type=float, default=defaults["mutation_instability_jump"])
    p.add_argument("--mutation-buildup-gain", type=float, default=defaults["mutation_buildup_gain"])
    p.add_argument("--use-logistic", dest="use_logistic", action=argparse.BooleanOptionalAction, default=defaults["use_logistic"])
    p.add_argument("--use-logistic-adapted", dest="use_logistic_adapted", action=argparse.BooleanOptionalAction, default=defaults["use_logistic_adapted"])
    p.add_argument("--verbose", action=argparse.BooleanOptionalAction, default=defaults["verbose"])
    p.add_argument("--scale", action=argparse.BooleanOptionalAction, default=defaults["scale"])
    p.add_argument("--decay", action=argparse.BooleanOptionalAction, default=defaults["decay"])
    p.add_argument("--save-all-steps", dest="save_all_steps", action=argparse.BooleanOptionalAction, default=defaults["save_all_steps"])
    p.add_argument("--save-interval", dest="save_interval", type=int, default=defaults["save_interval"])
    p.add_argument("--omega", dest="OMEGA", type=int, default=defaults["OMEGA"])
    return p


def get_explicit_cli_args(namespace: argparse.Namespace, parser: argparse.ArgumentParser) -> dict:
    defaults = vars(parser.parse_args([]))
    config_fields = _config_field_names()
    result = {}
    for k, v in vars(namespace).items():
        if k in config_fields:
            if k in defaults and v != defaults[k]:
                result[k] = v
            elif k not in defaults:
                result[k] = v
    return result


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    return parser.parse_args(args)
