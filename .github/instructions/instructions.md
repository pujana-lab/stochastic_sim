---
applyTo: '**'
---

# Instrucciones para Copilot - Proyecto Moran

## 1. Rol y Responsabilidades

Eres un **desarrollador competente** con expertise en modelado de poblaciones y probabilidad. Tu objetivo es:

- Escribir código de **alta calidad**, limpio y mantenible
- Crear código **testeable** con pruebas unitarias exhaustivas siguiendo TDD
- Documentar claramente para que matemáticos y biólogos lo entiendan
- Mantener una estructura de código **modular y escalable**

## 2. Contexto del Proyecto

Este proyecto implementa el **Proceso de Moran**, un modelo estocástico de evolución de poblaciones finitas.

En cada paso del proceso:
1. Un individuo es seleccionado para **reproducirse** (proporcional al fitness)
2. Un individuo es seleccionado para **morir** (reemplazado por la copia del reproductor)
3. El tamaño total de la población se mantiene constante

Referencias:
- https://en.wikipedia.org/wiki/Moran_process

## 3. Principios de Desarrollo

### 3.1 Legibilidad y Claridad

- Usa nombres de variables descriptivos en inglés: `reproductor_group`, `victim_group`, no `rep`, `v`
- Escribe docstrings en formato Google para todas las funciones y clases
- Máximo 120 caracteres por línea
- Agrupa imports: stdlib → third-party → local

### 3.2 Type Hints y Documentación

```python
# ✅ CORRECTO
from src.population import Population

def get_reproductor_group(self, population: Population) -> str:
    """
    Selecciona el grupo que se reproduce en este paso.

    Args:
        population: Estado actual de la población.

    Returns:
        str: Nombre del grupo reproductor.

    Example:
        >>> engine = EvolutionEngineDeterministic()
        >>> engine.get_reproductor_group(population)
        'dominant'
    """
    return max(
        population.groups.values(),
        key=lambda g: g.fitness * g.n
    ).name
```

### 3.3 Estructura de Proyecto

```
moran/
├── src/
│   ├── __init__.py
│   ├── population.py                          # SubPopulation + Population
│   ├── evolution_engine_interface.py          # Interfaz abstracta EvolutionEngineInterface
│   ├── evolution_engine_deterministic.py      # Motor determinista (selección por fitness × n)
│   ├── evolution_engine_bernoulli.py          # Motor probabilístico (selección proporcional al fitness)
│   └── simulator.py                           # Simulator: orquesta pasos y guarda tracking
├── tests/
│   ├── __init__.py
│   ├── test_population_kata_01_minimal_builder.py          # Kata 1: construcción básica de Population
│   ├── test_population_kata_02_basic_fitness.py            # Kata 2: fitness en SubPopulation
│   ├── test_population_kata_03_deterministic_evolution.py  # Kata 3: motor determinista
│   ├── test_population_kata_04_simulation_tracking.py      # Kata 4: tracking con Simulator
│   ├── test_population_kata_05_refactor_determinic_evolution.py  # Kata 5: refactor "tell, don't ask"
│   └── test_population_kata_06_bernoulli_evolution_engine.py     # Kata 6: motor Bernoulli
├── main.py              # Punto de entrada (pendiente de implementar)
├── makefile
├── requirements
└── README.md
```

#### Responsabilidades por módulo

| Módulo | Clase | Responsabilidad |
|--------|-------|-----------------|
| `population.py` | `SubPopulation` | Dataclass que almacena `name`, `n` (tamaño) y `fitness` de un subgrupo. |
| `population.py` | `Population` | Contiene un dict de `SubPopulation`. Expone `individuals`, `fitness`, `n`. Métodos: `append_individual`, `remove_individual`. |
| `evolution_engine_interface.py` | `EvolutionEngineInterface` | Contrato abstracto con `get_reproductor_group(population)` y `get_victim_group(population)`. |
| `evolution_engine_deterministic.py` | `EvolutionEngineDeterministic` | Selecciona reproductor (max `fitness × n`) y víctima (min `fitness × n`) de forma determinista. |
| `evolution_engine_bernoulli.py` | `EvolutionEngineBernoulli` | Selecciona reproductor y víctima de forma probabilística proporcional al fitness. |
| `simulator.py` | `Simulator` | Orquesta pasos con `run()`. Guarda histórico en `tracking`. Expone `get_tracking_dataframe()` y `get_tracking_summary_df()`. |

## 4. Modelos de Datos

### SubPopulation

```python
@dataclass
class SubPopulation:
    name: str       # Identificador del subgrupo
    n: int          # Número de individuos
    fitness: float  # Aptitud (por defecto 1.0)
```

### Population

```python
# Construcción esperada
groups = {
    'dominant': SubPopulation('dominant', n=60, fitness=2.0),
    'weak':     SubPopulation('weak',     n=40, fitness=1.0),
}
population = Population(groups)
```

### Simulator

```python
simulator = Simulator(population, evol=EvolutionEngineDeterministic())
simulator.run()                          # Ejecuta un paso
simulator.get_tracking_summary_df()     # DataFrame con evolución temporal
```

## 5. Testing - Requisito Obligatorio

### 5.1 Estándares de Testing

Toda función debe tener pruebas que cubran:

- **Caso base**: Comportamiento normal esperado
- **Casos límite**: Población de un solo grupo, n=1, fitness=0
- **Casos de error**: Inputs inválidos

### 5.2 Estructura de Tests

Los tests siguen una nomenclatura de **kata incremental**:
`test_population_kata_NN_descripcion.py`, donde `NN` indica la fase.

```python
# ✅ CORRECTO
import pytest
from src.population import Population, SubPopulation
from src.evolution_engine_deterministic import EvolutionEngineDeterministic

class TestDeterministicEvolution:
    """Suite de pruebas para el motor de evolución determinista."""

    @pytest.fixture
    def population(self):
        """Fixture: población con dos grupos y fitness distinto."""
        return Population({
            'dominant': SubPopulation('dominant', n=40, fitness=2.0),
            'weak':     SubPopulation('weak',     n=60, fitness=1.0),
        })

    def test_reproductor_is_highest_fitness_group(self, population):
        """Prueba: el grupo con mayor fitness × n se reproduce."""
        engine = EvolutionEngineDeterministic()
        assert engine.get_reproductor_group(population) == 'dominant'

    def test_victim_is_lowest_fitness_group(self, population):
        """Prueba: el grupo con menor fitness × n pierde un individuo."""
        engine = EvolutionEngineDeterministic()
        assert engine.get_victim_group(population) == 'weak'
```

### 5.3 Comandos de Testing

```bash
# Ejecutar todas las pruebas
pytest

# Con cobertura
pytest --cov=src --cov-report=html

# Tests específicos
pytest tests/test_population_kata_03_deterministic_evolution.py -v

# Con output detallado
pytest -v --tb=short
```

## 6. Estándares de Código Python

### 6.1 Imports y Organización

```python
# ✅ CORRECTO
import logging
from abc import ABC, abstractmethod

import pandas as pd

from src.population import Population, SubPopulation
from src.evolution_engine_interface import EvolutionEngineInterface

# ❌ INCORRECTO
from src.population import *
```

### 6.2 Manejo de Errores

```python
# ✅ CORRECTO
def remove_individual(self, group: str) -> None:
    if group not in self.individuals:
        raise ValueError(f"El grupo '{group}' no tiene individuos para eliminar.")
    self.individuals.remove(group)
    self.groups[group].n -= 1
    self.n -= 1
```

## 7. Logging y Debugging

```python
import logging

logger = logging.getLogger(__name__)

def run(self):
    logger.debug(f"Estado antes del paso: {self.population.groups}")
    reproductor = self.evol.get_reproductor_group(self.population)
    victim = self.evol.get_victim_group(self.population)
    logger.info(f"Paso: reproductor={reproductor}, víctima={victim}")
```

## 8. Checklist Antes de Commitear

- [ ] Código sigue PEP 8 y estándares del proyecto
- [ ] Todas las funciones tienen docstrings completos
- [ ] Type hints en todas las funciones
- [ ] Tests cubren casos base, límite y error
- [ ] Tests pasan: `pytest --cov=src`
- [ ] Cobertura ≥ 80%
- [ ] Nombres descriptivos (no abreviaturas)
- [ ] Sin código comentado (borrar o explicar)
- [ ] Logs apropiados con niveles correctos

## 9. Dependencias del Proyecto

```
pytest>=6.2.0
pytest-cov>=2.12.0
pandas>=1.3.0
```