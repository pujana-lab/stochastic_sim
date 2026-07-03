from __future__ import annotations

from src.gillespie.infrastructure.config.cli_config import build_parser, parse_args
from src.gillespie.application.config_service import merge_and_build
from src.gillespie.infrastructure.config.json_config import load_config_json, flatten_cell_types
from src.gillespie.infrastructure.config.cli_config import get_explicit_cli_args
from src.gillespie.infrastructure.display.summary_printer import print_summary
from src.gillespie.infrastructure.csv_output import save_rates_history_csv, save_history_csv, save_clones_csv, save_debug_history_csv
from src.gillespie.tumor_simulation import TumorSimulation
from src.gillespie.simulation_config import SimulationConfig
from pathlib import Path
# EStoy un poco harto del tema del cli la verdad, es complicarse la vida para nada por lo menos antes estaba todo en un archivo ahora es un puto lio

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
    save_rates_history_csv(Path("./rate_his.csv"), rates_history)
    print(f"Saved rates history to {Path("./rate_his.csv")}")
    
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
