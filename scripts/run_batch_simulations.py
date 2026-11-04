from __future__ import annotations

import argparse
from pathlib import Path

from src.gillespie.infrastructure.csv_output import (
    save_clones_parquet,
    save_history_parquet,
    save_rates_history_parquet,
)
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation
from src.gillespie.infrastructure.config.json_config import load_config_json
from src.gillespie.infrastructure.run_config_output import save_run_config


def build_base_params(config_path: Path | None) -> dict:
    base: dict = {}
    if config_path is None:
        return base

    raw = load_config_json(config_path)

    for key, value in raw.items():
        if key == "cell_parameters":
            base[key] = value
        elif key in SimulationConfig.__dataclass_fields__:
            base[key] = value

    base.pop("seed", None)
    return base


def run_single_simulation(seed: int, base_params: dict, output_dir: Path, save_config: bool = True) -> Path:
    params = dict(base_params)
    params["seed"] = seed
    config = SimulationConfig(**params)

    sim = TumorSimulation(config=config)
    times, history, tissue_state, rates_history = sim.run()

    run_dir = output_dir / f"seed_{seed:04d}"
    run_dir.mkdir(parents=True, exist_ok=True)

    if save_config:
        save_run_config(run_dir / "config", config)   # config.txt (table) + config.json (reproducible)

    save_history_parquet(run_dir / "history.parquet", times, history)
    save_clones_parquet(run_dir / "clones.parquet", tissue_state.clones)
    save_rates_history_parquet(run_dir / "rates_history.parquet", rates_history)

    return run_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run the same Gillespie simulation multiple times with different RNG seeds."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/gillespie_defaults.json"),
        help="JSON configuration file used as base parameters for every simulation.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/multi_seed_runs"),
        help="Directory where each seed run is saved.",
    )
    parser.add_argument(
        "--runs",
        type=int,
        default=10,
        help="How many simulations to launch.",
    )
    parser.add_argument(
        "--start-seed",
        type=int,
        default=1,
        help="First RNG seed to use.",
    )
    parser.add_argument(
        "--save-config",
        action="store_true",
        default=True,
        help="Save each run's conditions (input and used values) to config.txt and config.json inside the seed directory.",
    )
   
    parser.add_argument(
        "--no-save-config",
        dest="save_config",
        action="store_false"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base_params = build_base_params(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Base parameters loaded from: {args.config}")
    print(f"Running {args.runs} simulations starting from seed {args.start_seed}")
    print(f"Output directory: {args.output_dir}")
    print(f"Save config snapshot: {args.save_config}")
    print("Starting...")
    for i in range(args.runs):
        seed = args.start_seed + i
        run_dir = run_single_simulation(
            seed=seed,
            base_params=base_params,
            output_dir=args.output_dir,
            save_config=args.save_config,
        )
        print(f"[{i + 1}/{args.runs}] seed={seed} -> {run_dir}")


if __name__ == "__main__":
    main()
