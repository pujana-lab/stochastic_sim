from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Dict, Any, List

from src.gillespie.infrastructure.csv_output import (
    save_clones_parquet,
    save_history_parquet,
    save_rates_history_parquet,
)
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


# 5 Combinaciones distintas de (gamma, omega)
PARAM_COMBINATIONS = [
    # {"gamma": 0.10, "omega": 0.50, "name": "azul"},
    # {"gamma": 0.10, "omega": 0.80, "name": "roja"},
    # {"gamma": 0.7, "omega": 0.30, "name": "verde"},
    {"gamma": 0.30, "omega": 0.90,"N_immune": 100, "N_mutated": 100, "name": "roj_verde_1"},
    {"gamma": 0.30, "omega": 0.90,"N_immune": 1200, "N_mutated": 100, "name": "roj_verde_2"},
]


def build_config_for_combination(
    base_params: dict, gamma: float, omega: float, seed: int = 42,N_immune: int = 100, N_mutated: int = 100
) -> SimulationConfig:
    """Crea una instancia de SimulationConfig actualizando gamma y omega_exhaust."""
    params = copy.deepcopy(base_params)
    params["seed"] = seed

    # Asignar gamma a la tasa de muerte tumoral (mu) en base y mutated
    if "base" in params["cell_parameters"]:
        params["cell_parameters"]["base"]["mu"] = gamma
    if "mutated" in params["cell_parameters"]:
        params["cell_parameters"]["mutated"]["mu"] = gamma
        params["cell_parameters"]["mutated"]["N"] = N_mutated  # Valor por defecto si no se especifica

    # Asignar omega a la tasa de exhaustacion en celulas inmunes
    if "immune" in params["cell_parameters"]:
        params["cell_parameters"]["immune"]["omega_exhaust"] = omega
        params["cell_parameters"]["immune"]["N"] = N_immune  # Valor por defecto si no se especifica
    return SimulationConfig(**params)


def save_config_snapshot(config: SimulationConfig, path: Path) -> None:
    """Guarda un resumen de texto claro de la configuración utilizada."""
    lines: list[str] = ["SimulationConfig Snapshot"]
    for field_name in config.__dataclass_fields__:
        if field_name == "crowding_strategy":
            continue
        value = getattr(config, field_name)
        if field_name == "cell_parameters":
            lines.append("cell_parameters:")
            for cell_type, cell_cfg in value.items():
                lines.append(f"  {cell_type}: {cell_cfg}")
        else:
            lines.append(f"{field_name}={value!r}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def load_base_params(config_path: Path | None) -> dict:
    """Carga parametros base desde JSON o usa los defaults de SimulationConfig."""
    if config_path and config_path.exists():
        with config_path.open() as f:
            raw = json.load(f)
        return raw

    # Valores base por defecto si no se especifica JSON
    

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Ejecuta 5 simulaciones con combinaciones de gamma (muerte tumoral) y omega (exhaustion)."
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Ruta opcional a un archivo JSON de configuracion base.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/gamma_omega_sweep"),
        help="Directorio donde se guardaran los resultados de cada combinacion.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Semilla aleatoria para reproducibilidad.",
    )
    args = parser.parse_args()

    base_params = load_base_params(args.config)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    summary_results: List[Dict[str, Any]] = []

    print("=" * 70)
    print(" EJECUTANDO BARRIDO DE PARAMETROS (5 COMBINACIONES GAMMA / OMEGA)")
    print("=" * 70)

    for idx, combo in enumerate(PARAM_COMBINATIONS, 1):
        gamma = combo["gamma"]
        omega = combo["omega"]
        N_immune = combo["N_immune"]
        N_mutated = combo["N_mutated"]
        name = combo["name"]

        print(f"\n[{idx}/5] Ejecutando: {name} | Gamma (muerte)={gamma} | Omega (exhaustion)={omega}")

        # Construir configuración e inicializar la simulación
        config = build_config_for_combination(base_params, gamma=gamma, omega=omega, seed=args.seed, N_immune=N_immune, N_mutated=N_mutated)
        sim = TumorSimulation(config=config)

        # Ejecutar la simulación
        times, history, tissue_state, rates_history = sim.run()

        # Guardar artefactos
        run_dir = args.output_dir / name
        run_dir.mkdir(parents=True, exist_ok=True)

        save_config_snapshot(config, run_dir / "config.txt")
        save_history_parquet(run_dir / "history.parquet", times, history)
        save_clones_parquet(run_dir / "clones.parquet", tissue_state.clones)
        save_rates_history_parquet(run_dir / "rates_history.parquet", rates_history)

        # Extraer poblaciones finales para el resumen
        final_pops = tissue_state.pop_map
        summary_results.append(
            {
                "Combo": name,
                "Gamma": gamma,
                "Omega": omega,
                "Final_Time": f"{sim.t:.2f}",
                "Base": final_pops.get("base", 0),
                "Mutated": final_pops.get("mutated", 0),
                "Immune": final_pops.get("immune", 0),
                "Exhausted": final_pops.get("exhausted", 0),
            }
        )

    # Imprimir resumen en consola
    print("\n" + "=" * 85)
    print(" RESUMEN DE RESULTADOS DE LAS 5 SIMULACIONES")
    print("=" * 85)
    header = f"{'Nombre Combo':<32} | {'Gamma':<6} | {'Omega':<6} | {'T_Final':<8} | {'Base':<5} | {'Mut':<5} | {'Imm':<5} | {'Exh':<5}"
    print(header)
    print("-" * len(header))
    for res in summary_results:
        print(
            f"{res['Combo']:<32} | {res['Gamma']:<6.2f} | {res['Omega']:<6.2f} | "
            f"{res['Final_Time']:<8} | {res['Base']:<5} | {res['Mutated']:<5} | "
            f"{res['Immune']:<5} | {res['Exhausted']:<5}"
        )
    print("-" * len(header))


if __name__ == "__main__":
    main()