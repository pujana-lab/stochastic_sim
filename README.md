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

---

## Estructura del Proyecto

```
moran/
├── src/
│   ├── population.py                          # SubPopulation + Population
│   ├── evolution_engine_interface.py          # Interfaz abstracta EvolutionEngineInterface
│   ├── evolution_engine_deterministic.py      # Motor determinista (selección por fitness × n)
│   ├── evolution_engine_bernoulli.py          # Motor probabilístico (selección proporcional al fitness)
│   └── simulator.py                           # Simulator: orquesta pasos y guarda tracking
├── tests/
│   ├── test_population_kata_01_minimal_builder.py
│   ├── test_population_kata_02_basic_fitness.py
│   ├── test_population_kata_03_deterministic_evolution.py
│   ├── test_population_kata_04_simulation_tracking.py
│   ├── test_population_kata_05_refactor_determinic_evolution.py
│   └── test_population_kata_06_bernoulli_evolution_engine.py
├── main.py              # Punto de entrada CLI
├── config.yaml          # Ejemplo de configuración
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
| `simulator.py` | `Simulator` | Orquesta pasos con `run()`. Guarda histórico en `tracking`. Expone `get_tracking_dataframe()` y `get_tracking_summary_df()`. |

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
python main.py --config config.yaml
```

### Ejemplo de `config.yaml`

```yaml
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
```

### Desde Python

```python
from src.population import Population, SubPopulation
from src.evolution_engine_deterministic import EvolutionEngineDeterministic
from src.simulator import Simulator

population = Population({
    'dominant': SubPopulation('dominant', n=60, fitness=2.0),
    'weak':     SubPopulation('weak',     n=40, fitness=1.0),
})

simulator = Simulator(population=population, evol=EvolutionEngineDeterministic())

for _ in range(100):
    simulator.run()

df = simulator.get_tracking_summary_df()
df.to_excel('results.xlsx', index=False)
```

---

## Testing

El proyecto sigue **TDD** con una nomenclatura de kata incremental:

```bash
# Ejecutar todos los tests
pytest

# Con cobertura
pytest --cov=src --cov-report=html

# Test específico
pytest tests/test_population_kata_03_deterministic_evolution.py -v
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

---

## Desarrollo

Este repositorio se construye siguiendo **TDD**: las estructuras de datos y funciones se desarrollan de manera incremental, asegurando que cada componente pase sus pruebas antes de avanzar al siguiente.

