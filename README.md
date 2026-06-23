# Tumor Gillespie Simulator

Simulación de clones tumorales mediante un algoritmo de Gillespie con eventos de nacimiento, muerte, mutación y agotamiento.

## Visión general

Este proyecto modela un tejido con cuatro tipos de clones:

- `base` (células sanas o normales)
- `mutated` (células tumorales)
- `immune` (células inmunes)
- `exhausted` (células agotadas)

Cada clon tiene tasas base de nacimiento, muerte, mutación y agotamiento. Las tasas efectivas se calculan en función del estado global del tejido (`TissueState`) y de interacciones entre los tipos celulares.

## Matemáticas del modelo

### Eventos del modelo

El simulador considera los siguientes eventos de Gillespie para cada clon vivo:

- `BIRTH`: división celular
- `DEATH`: muerte celular
- `MUTATION`: creación de un clon mutado hijo
- `EXHAUSTION`: transformación entre clones inmune/exhaustos

### Tasas efectivas por clon

Para un clon de tipo `i` con población total `N_i` en el estado actual:

- tasa de nacimiento:
  $$ r^B_i = \lambda_i \, N_i \, C_i(t) $$
- tasa de muerte:
  $$ r^D_i = \mu_i \, N_i + I_i(t) $$
- tasa de mutación:
  $$ r^M_i = \nu_i \, N_i \, (1 + \text{instability}_i) $$
- tasa de agotamiento:
  $$ r^E_i = \epsilon_i \, N_i $$

Aquí `C_i(t)` es el factor de `crowding` y `I_i(t)` es un término de interacción extra dependiente de otros tipos celulares.

#### Crowding

Las dos estrategias de crowding implementadas son:

- `SimpleCrowding`
- `AdaptedCrowding`

Ambas usan la forma general:

$$ C_i(t) = \max\left(0, 1 - \frac{N_{crowd}}{K_i(t)} \right) $$

con:

- `N_{crowd}`: población competitiva para ese clon
- `K_i(t)`: capacidad de carga efectiva en tiempo `t`

Para `SimpleCrowding`:

$$ K_i(t) = \max(K_{min}, K_i - decline \cdot t) $$

Para `AdaptedCrowding`:

$$ K_i(t) = \max\left(K_{min}, \frac{K_i}{1 - \mu_i / \lambda_i} - decline \cdot t \right) $$

### Interacciones específicas

#### `MutatedClone`

La muerte efectiva incluye un término inmune:

$$ r^D_{mut} = \mu_{mut} N_{mut} + N_{mut} N_{immune} \theta_I $$

#### `ImmuneClone`

- nacimiento activo:
  $$ r^B_{immune} = \lambda_{immune} N_{immune} C_{immune} + N_{immune} N_{mut} \beta $$
- agotamiento dependiente del tumor:
  $$ r^E_{immune} = \epsilon_{immune} N_{immune} N_{mut} $$

#### `ExhaustedClone`

- no se divide:
  $$ r^B_{exh} = 0 $$

### Evolución de la inestabilidad

Cada clon vivo acumula inestabilidad en el tiempo:

$$ \mathit{instability}_i(t + \Delta t) = \mathit{instability}_i(t) + (\mathit{base\_instability\_buildup} + \mathit{buildup}_i) \cdot \Delta t $$

Esto afecta la tasa de mutación en `MutatedClone` y otros eventos si se extiende el modelo.

### Gillespie paso a paso

En cada paso de la simulación:

1. Actualizar `pop_map` en `TissueState`
2. Construir la matriz de tasas de eventos (`RateMatrix`)
3. Sumar la tasa total
   $$ R = \sum_{j} r_j $$
4. Muestrear el tiempo hasta el próximo evento:
   $$ \tau = -\frac{1}{R} \, \ln(u), \quad u \sim U(0,1) $$
5. Seleccionar un evento proporcional a su tasa
6. Aplicar el evento al clon elegido
7. Actualizar el tiempo y guardar un snapshot del estado

La simulación termina cuando se alcanza `T_max` o no hay más eventos posibles.

## Arquitectura del código

### `src/gillespie/simulation_config.py`

Define los parámetros globales del modelo:

- poblaciones iniciales: `N0`, `N_mutant`, `N_immune`, `N_exhausted`
- tasas base: `lambda0`, `lambda_Immune`, `mu0`, `omega_exhaust`, `mu_Exhausted`, `nu0`
- parámetros de interacción: `theta_I`, `beta`
- capacidades: `K0`, `K_immune`, `K_mutant`, `decline`, `Kmin`
- tiempo máximo: `T_max`
- parámetros de inestabilidad y fitness

### `src/gillespie/clone_factory.py`

Construye los clones iniciales y asigna sus parámetros base según el tipo:

- `base` → `WildTypeClone`
- `mutated` → `MutatedClone`
- `immune` → `ImmuneClone`
- `exhausted` → `ExhaustedClone`

### `src/gillespie/clone.py`

Define la clase base `Clone` y subclases especializadas.
Cada clon implementa métodos de tasa efectiva que reciben el `TissueState` completo y pueden leer interacciones globales desde `self.config`.

### `src/gillespie/tissue_state.py`

Encapsula el estado actual del tejido:

- `clones: Dict[CloneId, Clone]`
- `pop_map: Dict[CloneType, int]`
- `snapshot()` para historial
- consultas por tipo de clon

### `src/gillespie/crowding_strategy.py`

Implementa las estrategias de crowding que modifican la tasa de nacimiento.

### `src/gillespie/tumor_simulation.py`

Controla la simulación completa:

- instancia `CloneFactory` y crea el `TissueState`
- construye la matriz de eventos
- muestrea tiempos con la fórmula de Gillespie
- aplica eventos y guarda el historial

## Uso sencillo

```python
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation

config = SimulationConfig(
    N0=100,
    N_mutant=5,
    N_immune=20,
    N_exhausted=0,
    lambda0=0.5,
    mu0=0.2,
    nu0=0.01,
    theta_I=0.001,
    beta=0.001,
    K0=1000,
    K_mutant=1000,
    K_immune=500,
    use_logistic=True,
    use_logistic_adapted=True,
    T_max=1000,
)

sim = TumorSimulation(config)
results = sim.run()
```

## Reescalado de parámetros

Si quieres adaptar el modelo al tamaño del sistema, prepara la configuración antes de crear `TumorSimulation`.

```python
scaled_config = SimulationConfig(
    N0=int(base_N0 * scale),
    K0=int(base_K0 * scale),
    theta_I=base_theta_I / scale,
    beta=base_beta / scale,
    # ... otros parámetros según tu estrategia de escalado
)
```

## Consideraciones

- `TissueState` contiene el estado del tejido; no debe mezclar lógica de simulación.
- El cálculo de tasas y eventos debe realizarse en `TumorSimulation`.
- Si clones del mismo tipo comparten tasas efectivas, conviene calcularlas una vez por tipo en `_build_rate_matrix()`.

## Testing

```bash
python -m pytest tests/gillespie
```

O:

```bash
make test
```
