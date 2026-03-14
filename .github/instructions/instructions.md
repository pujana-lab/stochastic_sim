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

El simulador también soporta **mutaciones**: en paralelo al paso de evolución, un motor de mutaciones puede transformar un individuo de un grupo existente en un individuo de un nuevo grupo con fitness distinto. El evento queda registrado en la columna `events` del tracking.

### Estado actual de mutaciones

- La columna `events` está implementada y registra eventos de mutación.
- El `Simulator` acepta un único `mutation_engine`.
- **Pendiente (kata 09+):** soporte para lista de `mutation_engines` simultáneos y columna `lineage` para reconstruir el árbol genealógico de grupos.

### Objetivo a medio plazo: linaje poligénico básico

El diseño previsto permite construir cadenas de mutación entre grupos (`weak → mutant_A → mutant_B`):

```python
simulator = Simulator(
    population=population,
    evol_engine=engine,
    mutation_engines=[
        MutationEngineBernoulli(p=0.05, victim='weak',     new_group='mutant_A', new_fitness=3.0),
        MutationEngineBernoulli(p=0.02, victim='mutant_A', new_group='mutant_B', new_fitness=6.0),
    ]
)
```

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
│   ├── population.py                          # SubPopulation + Population
│   ├── evolution_engine_interface.py          # Interfaz abstracta EvolutionEngineInterface
│   ├── evolution_engine_deterministic.py      # Motor determinista (selección por fitness × n)
│   ├── evolution_engine_bernoulli.py          # Motor probabilístico (selección proporcional al fitness)
│   ├── mutation_engine_interface.py           # Interfaz abstracta MutationEngineInterface
│   ├── mutation_engine_disabled.py            # Motor sin mutaciones (por defecto)
│   ├── mutation_engine_deterministic.py       # Mutación cada N pasos
│   ├── mutation_engine_bernoulli.py           # Mutación con probabilidad p
│   └── simulator.py                           # Simulator: orquesta evolución + mutación
├── tests/
│   ├── test_population_kata_01_minimal_builder.py
│   ├── test_population_kata_02_basic_fitness.py
│   ├── test_population_kata_03_deterministic_evolution.py
│   ├── test_population_kata_04_simulation_tracking.py
│   ├── test_population_kata_05_refactor_determinic_evolution.py
│   ├── test_population_kata_06_bernoulli_evolution_engine.py
│   ├── test_population_kata_07_deterministic_mutation_engine.py
│   └── test_population_kata_08_bernoulli_mutation_engine.py
├── configs/
│   ├── deterministic.yaml
│   ├── deterministic_mutation.yaml
│   └── bernoulli_mutation.yaml
├── main.py
├── Makefile
└── requirements.txt
```

#### Responsabilidades por módulo

| Módulo | Clase | Responsabilidad |
|--------|-------|-----------------|
| `population.py` | `SubPopulation` | Dataclass con `name`, `n` y `fitness` de un subgrupo. |
| `population.py` | `Population` | Contiene un dict de `SubPopulation`. Gestiona `individuals` y expone `append_individual` / `remove_individual`. |
| `evolution_engine_interface.py` | `EvolutionEngineInterface` | Contrato abstracto: `get_reproductor_group()` y `get_victim_group()`. |
| `evolution_engine_deterministic.py` | `EvolutionEngineDeterministic` | Selecciona reproductor (`max fitness × n`) y víctima (`min fitness × n`). Devuelve `None` si la población está fijada. |
| `evolution_engine_bernoulli.py` | `EvolutionEngineBernoulli` | Selecciona reproductor y víctima de forma probabilística proporcional al fitness. Acepta `seed`. |
| `mutation_engine_interface.py` | `MutationEngineInterface` | Contrato abstracto: `should_mutate(step)` y `get_mutation(population)`. |
| `mutation_engine_disabled.py` | `MutationEngineDisabled` | Implementación nula: nunca muta. Usado por defecto cuando no se configura mutación. |
| `mutation_engine_deterministic.py` | `MutationEngineDeterministic` | Muta cada `every_n_steps` pasos de forma determinista. |
| `mutation_engine_bernoulli.py` | `MutationEngineBernoulli` | Muta con probabilidad `p` en cada paso. Acepta `seed`. |
| `simulator.py` | `Simulator` | Orquesta pasos con `run()`. Aplica evolución y mutación. Guarda histórico en `tracking` con columna `events`. |

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
population = Population({
    'dominant': SubPopulation('dominant', n=60, fitness=2.0),
    'weak':     SubPopulation('weak',     n=40, fitness=1.0),
})
```

### Simulator (estado actual)

```python
simulator = Simulator(
    population=population,
    evol_engine=EvolutionEngineBernoulli(seed=42),
    mutation_engine=MutationEngineDeterministic(
        every_n_steps=10,
        victim_group='weak',
        new_group_name='mutated',
        new_fitness=5.0,
    ),
)
simulator.run()                          # Ejecuta un paso (evolución + mutación)
simulator.get_tracking_summary_df()     # DataFrame con columnas: time, grupos..., events
```

### MutationEngineInterface

```python
class MutationEngineInterface(ABC):
    def should_mutate(self, step: int) -> bool: ...
    def get_mutation(self, population: Population) -> tuple[str, str, float]:
        """Devuelve (victim_group, new_group_name, new_fitness)"""
        ...
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
from src.mutation_engine_deterministic import MutationEngineDeterministic
from src.simulator import Simulator

class TestDeterministicMutation:

    @pytest.fixture
    def population(self):
        return Population({
            'dominant': SubPopulation('dominant', n=80, fitness=2.0),
            'weak':     SubPopulation('weak',     n=20, fitness=1.0),
        })

    @pytest.fixture
    def mutation_engine(self):
        return MutationEngineDeterministic(
            every_n_steps=10,
            victim_group='weak',
            new_group_name='mutated',
            new_fitness=5.0,
        )

    def test_should_mutate_at_correct_steps(self, mutation_engine):
        assert mutation_engine.should_mutate(10) is True
        assert mutation_engine.should_mutate(5) is False

    def test_get_mutation_returns_correct_values(self, mutation_engine, population):
        assert mutation_engine.get_mutation(population) == ('weak', 'mutated', 5.0)
```

### 5.3 Katas


*Nota*: Es posible que, a medida de que se incrementa la complejidad del código, cambie la estructura de los test, dejando de hacer katas.

| Kata | Fichero | Contenido |
|------|---------|-----------|
| 01 | `test_population_kata_01_minimal_builder.py` | Construcción básica de `Population` sin fitness |
| 02 | `test_population_kata_02_basic_fitness.py` | `SubPopulation` con fitness en el constructor |
| 03 | `test_population_kata_03_deterministic_evolution.py` | Motor determinista: reproductor y víctima |
| 04 | `test_population_kata_04_simulation_tracking.py` | `Simulator`: tracking y summary DataFrame |
| 05 | `test_population_kata_05_refactor_determinic_evolution.py` | Refactor "tell, don't ask" |
| 06 | `test_population_kata_06_bernoulli_evolution_engine.py` | Motor Bernoulli probabilístico |
| 07 | `test_population_kata_07_deterministic_mutation_engine.py` | Motor mutación determinista + columna `events` |
| 08 | `test_population_kata_08_bernoulli_mutation_engine.py` | Motor mutación Bernoulli |
| 09 | _pendiente_ | `Simulator` acepta lista de `mutation_engines` |
| 10 | _pendiente_ | `events` registra `origin → new_group` por paso |
| 11 | _pendiente_ | `get_lineage_df()`: árbol completo de grupos |

### 5.4 Comandos de Testing

```bash
make test        # todos los tests
make test-cov    # con cobertura HTML
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
from src.mutation_engine_interface import MutationEngineInterface
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
    reproductor = self.evol_engine.get_reproductor_group(self.population)
    victim = self.evol_engine.get_victim_group(self.population)
    logger.info(f"Paso: reproductor={reproductor}, víctima={victim}")
```

## 8. Checklist Antes de Commitear

- [ ] Código sigue PEP 8 y estándares del proyecto
- [ ] Todas las funciones tienen docstrings completos
- [ ] Type hints en todas las funciones
- [ ] Tests cubren casos base, límite y error
- [ ] Tests pasan: `make test-cov`
- [ ] Cobertura ≥ 80%
- [ ] Nombres descriptivos (no abreviaturas)
- [ ] Sin código comentado (borrar o explicar)
- [ ] Logs apropiados con niveles correctos

## 9. Dependencias del Proyecto

```
pytest>=6.2.0
pytest-cov>=2.12.0
pandas>=1.3.0
pyyaml>=6.0
openpyxl>=3.0.0
```