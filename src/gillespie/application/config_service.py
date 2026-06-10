from __future__ import annotations

import argparse
from pathlib import Path

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.infrastructure.config.json_config import load_config_json
from src.gillespie.infrastructure.config.cli_config import build_parser, get_explicit_cli_args, parse_args


def merge_and_build(json_data: dict | None = None, cli_overrides: dict | None = None) -> SimulationConfig:
    merged = dict(json_data or {})
    merged.update(cli_overrides or {})
    return SimulationConfig(**merged)


def build_config_from_cli(cli_args: list[str] | None = None) -> tuple[SimulationConfig, argparse.Namespace]:
    parser = build_parser()
    namespace = parse_args(cli_args) if cli_args is not None else parser.parse_args()
    base: dict = {}
    if namespace.config is not None:
        base = load_config_json(namespace.config)
    cli_overrides = get_explicit_cli_args(namespace, parser)
    base.update(cli_overrides)
    return SimulationConfig(**base), namespace
