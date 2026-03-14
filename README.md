# moran
Basic Moran Process Simulator

**Authors:**
- Luis Palomero <lpalomerol@gmail.com>
- Victor Manso <victor.mansov@gmail.com>

---

## ¿Qué es el Proceso de Moran?

El [Proceso de Moran](https://en.wikipedia.org/wiki/Moran_process) es un modelo estocástico de evolución de poblaciones finitas. En cada paso:

1. Un individuo es seleccionado para **reproducirse** (proporcional al fitness)
2. Un individuo es seleccionado para **morir**
3. El tamaño total de la población se mantiene constante

El proceso termina cuando un grupo fija toda la población (`n = 0` para el resto).

### Mutaciones

El simulador soporta un motor de mutaciones que actúa en paralelo al proceso de evolución. En cada paso, el motor de mutaciones puede transformar un individuo de un grupo existente en un individuo de un nuevo grupo con un fitness distinto. El evento queda registrado en la columna `events` del tracking.

### Linaje poligénico ⚠️ _planificado_

El objetivo a medio plazo es soportar un **linaje poligénico básico**: múltiples motores de mutación simultáneos que permitan construir cadenas de origen entre grupos (`weak → mutant_A → mutant_B`). Esto permitirá reconstruir el árbol genealógico de los grupos a partir del tracking.

El diseño previsto es:

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

---

## Estructura del Proyecto

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
│   ├── deterministic.yaml                     # Evolución determinista, sin mutación
│   ├── deterministic_mutation.yaml            # Evolución Bernoulli + mutación determinista
│   └── bernoulli_mutation.yaml                # Evolución Bernoulli + mutación Bernoulli
├── main.py
├── Makefile
└── requirements.txt
```

### Responsabilidades por módulo

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

---

## Instalación

```bash
pip install -r requirements.txt
```

Dependencias:
```
pytest>=6.2.0
pytest-cov>=2.12.0
pandas>=1.3.0
pyyaml>=6.0
openpyxl>=3.0.0
```

---

## Uso

### Desde CLI

```bash
python main.py --config configs/deterministic.yaml
```

O usando el `Makefile`:

```bash
make run-deterministic           # evolución determinista, sin mutación
make run-deterministic-mutation  # evolución Bernoulli + mutación determinista
make run-bernoulli-mutation      # evolución Bernoulli + mutación Bernoulli
make run-all                     # los tres escenarios en secuencia
```

### Estructura del `config.yaml`

```yaml
populations:
  dominant:
    n: 60
    fitness: 2.0
  weak:
    n: 40
    fitness: 1.0

steps: 200

engine:
  type: bernoulli   # deterministic | bernoulli
  seed: 42

mutation:           # opcional — si se omite, no hay mutaciones
  type: bernoulli   # deterministic | bernoulli
  p: 0.1            # solo para bernoulli
  every_n_steps: 10 # solo para deterministic
  seed: 42          # solo para bernoulli
  victim_group: weak
  new_group_name: mutated
  new_fitness: 5.0

output: results.xlsx
```

### Desde Python

```python
from src.population import Population, SubPopulation
from src.evolution_engine_bernoulli import EvolutionEngineBernoulli
from src.mutation_engine_deterministic import MutationEngineDeterministic
from src.simulator import Simulator

population = Population({
    'dominant': SubPopulation('dominant', n=60, fitness=2.0),
    'weak':     SubPopulation('weak',     n=40, fitness=1.0),
})

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

for _ in range(200):
    simulator.run()

df = simulator.get_tracking_summary_df()
df.to_excel('results.xlsx', index=False)
```

---

## Testing

El proyecto sigue **TDD** con una nomenclatura de kata incremental:

```bash
make test        # todos los tests
make test-cov    # con cobertura HTML
```

### Katas

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
| 09 | _pendiente_ | Lista de `mutation_engines` simultáneos |
| 10 | _pendiente_ | Columna `lineage` en tracking |
| 11 | _pendiente_ | `get_lineage_df()`: árbol completo de grupos |

---

## Roadmap

### Próximas katas: linaje poligénico

El objetivo es construir un **linaje poligénico básico** a partir del tracking de mutaciones. El diseño se basa en linaje de **grupos** (no de individuos), que es suficiente para el modelo y manejable en memoria.

| Kata | Objetivo |
|------|----------|
| 09 | `Simulator` acepta `mutation_engines: list` en lugar de uno solo |
| 10 | `events` registra `origin → new_group` por cada mutación ocurrida en el paso |
| 11 | `get_lineage_df()` devuelve el árbol de origen de cada grupo con su paso de aparición |

---

## Desarrollo

Este repositorio se construye siguiendo **TDD**: las estructuras de datos y funciones se desarrollan de manera incremental, asegurando que cada componente pase sus pruebas antes de avanzar al siguiente.

