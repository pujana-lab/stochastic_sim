"""
main.py - Punto de entrada del simulador del Proceso de Moran.

Uso:
    python main.py --config config.yaml

El fichero YAML define las poblaciones, parámetros de simulación y salida.
"""
import argparse
import logging
import sys

import yaml

from src.evolution_engine_bernoulli import EvolutionEngineBernoulli
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.population import Population, SubPopulation
from src.simulator import Simulator

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def load_config(config_path: str) -> dict:
    """Carga la configuración desde un fichero YAML.

    Args:
        config_path: Ruta al fichero YAML de configuración.

    Returns:
        dict: Configuración parseada.

    Raises:
        FileNotFoundError: Si el fichero no existe.
        ValueError: Si el YAML tiene campos obligatorios ausentes.
    """
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    required_fields = ["populations", "steps", "engine", "output"]
    missing = [field for field in required_fields if field not in config]
    if missing:
        raise ValueError(f"Campos obligatorios ausentes en el YAML: {missing}")

    return config


def build_population(populations_config: dict) -> Population:
    """Construye una Population a partir de la sección 'populations' del YAML.

    Args:
        populations_config: Dict con nombre -> {n, fitness}.

    Returns:
        Population: Población inicializada.
    """
    groups = {
        name: SubPopulation(
            name=name,
            n=int(info["n"]),
            fitness=float(info.get("fitness", 1.0)),
        )
        for name, info in populations_config.items()
    }
    return Population(groups)


def build_engine(engine_config: dict):
    """Construye el motor de evolución según la configuración.

    Args:
        engine_config: Dict con 'type' (deterministic|bernoulli) y 'seed' (opcional).

    Returns:
        EvolutionEngineInterface: Motor instanciado.

    Raises:
        ValueError: Si el tipo de motor es desconocido.
    """
    engine_type = engine_config.get("type", "deterministic").lower()
    seed = engine_config.get("seed", None)

    if engine_type == "deterministic":
        logger.info("Motor: determinista")
        return EvolutionEngineDeterministic()
    elif engine_type == "bernoulli":
        logger.info(f"Motor: Bernoulli (seed={seed})")
        return EvolutionEngineBernoulli(seed=seed)
    else:
        raise ValueError(f"Tipo de motor desconocido: '{engine_type}'. Usa 'deterministic' o 'bernoulli'.")


def run_simulation(config: dict) -> None:
    """Ejecuta la simulación y guarda el resultado en un fichero XLSX.

    Args:
        config: Configuración completa cargada desde el YAML.
    """
    population = build_population(config["populations"])
    engine = build_engine(config["engine"])
    steps = int(config["steps"])
    output_path = config["output"]

    logger.info(f"Población inicial: { {k: v.n for k, v in population.groups.items()} }")
    logger.info(f"Pasos: {steps}")

    simulator = Simulator(population=population, evol=engine)

    for step in range(steps):
        simulator.run()
        logger.debug(f"Paso {step + 1}: { {k: v.n for k, v in simulator.population.groups.items()} }")

    df = simulator.get_tracking_summary_df()
    df.to_excel(output_path, index=False)
    logger.info(f"Resultados guardados en: {output_path}")


def parse_args(argv=None) -> argparse.Namespace:
    """Parsea los argumentos de línea de comandos.

    Args:
        argv: Lista de argumentos (por defecto sys.argv).

    Returns:
        argparse.Namespace: Argumentos parseados.
    """
    parser = argparse.ArgumentParser(
        description="Simulador del Proceso de Moran",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplo de uso:
    python main.py --config config.yaml

Estructura del fichero config.yaml:
    populations:
      dominant:
        n: 60
        fitness: 2.0
      weak:
        n: 40
        fitness: 1.0
    steps: 100
    engine:
      type: bernoulli   # deterministic | bernoulli
      seed: 42
    output: results.xlsx
        """,
    )
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Ruta al fichero YAML de configuración.",
    )
    return parser.parse_args(argv)


def main(argv=None) -> None:
    """Punto de entrada principal."""
    args = parse_args(argv)

    try:
        config = load_config(args.config)
        run_simulation(config)
    except FileNotFoundError as e:
        logger.error(f"Fichero no encontrado: {e}")
        sys.exit(1)
    except ValueError as e:
        logger.error(f"Error de configuración: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error inesperado: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()