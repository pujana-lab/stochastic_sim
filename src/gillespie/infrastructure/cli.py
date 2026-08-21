from __future__ import annotations

from pathlib import Path

from src.gillespie.infrastructure.config.cli_config import build_parser, get_explicit_cli_args, parse_args
from src.gillespie.infrastructure.config.json_config import flatten_cell_types, load_config_json
from src.gillespie.infrastructure.csv_output import (
    save_clones_parquet,
    save_debug_history_parquet,
    save_history_parquet,
    save_rates_history_parquet,
)
from src.gillespie.infrastructure.display.summary_printer import print_summary
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


def main() -> None:
    parser = build_parser()
    args = parse_args()

    base: dict = {}
    if args.config is not None:
        base = flatten_cell_types(load_config_json(args.config))
        print(base)
    cli_overrides = get_explicit_cli_args(args, parser)
    base.update(cli_overrides)
    config = SimulationConfig(**base)

    sim = TumorSimulation(config=config)
    input("press any key to continue")
    print("Starting...")
    #TODO: mover los writers a dentro del simulador
    times, history, tissue_state,rates_history = sim.run()
# las mierdas que tengo que hacer para guardar un csv son locas. 
    print_summary(times, tissue_state.clones, args.top)
    rates_path = save_rates_history_parquet(Path("./rate_his.parquet"), rates_history)
    print(f"Saved rates history to {rates_path}")

    if getattr(args, "save_history", None) is not None:
        history_path = save_history_parquet(args.save_history, times, history)
        print(f"\nSaved history to {history_path}")

    if getattr(args, "save_clones", None) is not None:
        clones_path = save_clones_parquet(args.save_clones, tissue_state.clones)
        print(f"Saved clones to {clones_path}")

    if getattr(args, "save_debug", None) is not None:
        debug_path = save_debug_history_parquet(args.save_debug, times, history, sim.events)
        print(f"Saved debug trace to {debug_path}")


if __name__ == "__main__":
    main()
