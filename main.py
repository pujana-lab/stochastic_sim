"""
Proyecto Moran - Análisis e implementación del índice I de Moran
=========================================================================

Este módulo principal coordina los componentes del análisis espacial
usando el índice I de Moran para detección de autocorrelación espacial.

Estructura del proyecto:
- src/: Módulos de lógica principal
- tests/: Suite de pruebas unitarias
- data/: Datos de entrada
- results/: Salidas y visualizaciones

Autor: Luis
Fecha: 2026-03-13
"""

from pathlib import Path

# Constantes del proyecto
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
SRC_DIR = PROJECT_ROOT / "src"
TESTS_DIR = PROJECT_ROOT / "tests"


def main() -> None:
    """
    Función principal que orquesta el análisis espacial.
    
    El flujo esperado es:
    1. Cargar datos espaciales
    2. Calcular matriz de pesos espaciales (W)
    3. Computar índice I de Moran
    4. Realizar pruebas de significancia
    5. Generar visualizaciones
    """
    print("Iniciando análisis Moran...")
    # TODO: Implementar flujo principal
    

if __name__ == "__main__":
    main()
