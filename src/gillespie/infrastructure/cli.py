from __future__ import annotations

import sys

from src.gillespie.infrastructure.config.cli_config import build_parser, parse_args
from src.gillespie.application.config_service import merge_and_build
from src.gillespie.infrastructure.config.json_config import load_config_json, flatten_cell_types
from src.gillespie.infrastructure.config.cli_config import get_explicit_cli_args
from src.gillespie.infrastructure.display.summary_printer import print_summary
from src.gillespie.infrastructure.csv_output import save_history_csv, save_clones_csv, save_debug_history_csv
from src.gillespie.tumor_simulation import TumorSimulation
from src.gillespie.simulation_config import SimulationConfig


def main() -> None:
    if len(sys.argv) > 1 and sys.argv[1] == "calibrate":
        from src.gillespie.calibration.application.cli import calibrate_main
        calibrate_main(sys.argv[2:] if len(sys.argv) > 2 else None)
        return

    parser = build_parser()
    args = parse_args()

    base: dict = {}
    if args.config is not None:
        base = flatten_cell_types(load_config_json(args.config))
    cli_overrides = get_explicit_cli_args(args, parser)
    base.update(cli_overrides)
    config = SimulationConfig(**base)

    sim = TumorSimulation(config=config)
    times, history, tissue_state = sim.run()

    print_summary(times, tissue_state.clones, args.top)

    if args.save_history is not None:
        save_history_csv(args.save_history, times, history)
        print(f"\nSaved history to {args.save_history}")

    if args.save_clones is not None:
        save_clones_csv(args.save_clones, tissue_state.clones)
        print(f"Saved clones to {args.save_clones}")

    if args.save_debug is not None:
        save_debug_history_csv(args.save_debug, times, history, sim.events)
        print(f"Saved debug trace to {args.save_debug}")


if __name__ == "__main__":
    main()
