# Refactorización: Aislamiento del Tissue State

## Objetivo
Aislar el estado del tejido (conjunto de clones) como un objeto independiente (`TissueState`) para:
- Facilitar el historial de estados del sistema en el tiempo
- Pasar el estado completo a los clones para que calculen tasas efectivas basadas en interacciones
- Desacoplar la lógica del simulador de la lógica de cálculo de tasas

## Cambios Implementados

### 1. Nueva Clase: `TissueState` (`src/gillespie/tissue_state.py`)

```python
@dataclass
class TissueState:
    clones: Dict[CloneId, Clone]
```

**Métodos principales:**
- `total_population()` - Retorna el número total de células
- `population_by_type(cell_type)` - Retorna población total de un tipo específico
- `get_clone(clone_id)` - Obtiene un clon específico
- `get_clones_by_type(cell_type)` - Obtiene todos los clones de un tipo
- `snapshot()` - Crea un snapshot del estado actual para el historial

### 2. Modificaciones en `Clone` (`src/gillespie/clone.py`)

Los métodos de cálculo de tasas efectivas ahora reciben `TissueState`:

```python
# Antes
def birth_rate_effective(self, crowding: float) -> float
def death_rate_effective(self) -> float
def exhaustion_rate_effective(self) -> float

# Ahora
def birth_rate_effective(self, tissue_state: TissueState, crowding: float) -> float
def death_rate_effective(self, tissue_state: TissueState) -> float
def exhaustion_rate_effective(self, tissue_state: TissueState) -> float
```

Además, `Clone` mantiene una referencia a `self.config` para leer parámetros globales de interacción como `theta_I` y `beta`.
### 3. Subclases de Clone Actualizadas

#### `MutatedClone`
- Accede a `tissue_state.population_by_type("immune")` en `death_rate_effective()`.
- Usa el parámetro global `self.config.theta_I` para el término de killing inmune.

#### `ImmuneClone`
- Accede a `tissue_state.population_by_type("mutated")` en:
  - `birth_rate_effective()` - para calcular boost de activación
  - `exhaustion_rate_effective()` - para calcular agotamiento dependiente de cáncer
- Usa el parámetro global `self.config.beta` para activación.

#### `ExhaustedClone`
- Simplificada: siempre retorna 0 en `birth_rate_effective()`

### 4. Actualización de `TumorSimulation` (`src/gillespie/tumor_simulation.py`)

**Cambios clave:**
- Ahora mantiene `self.tissue_state: TissueState` en lugar de `self.clones: Dict`
- `_build_rate_matrix()` pasa `tissue_state` a los clones
- `step()` usa `tissue_state.snapshot()` para el historial
- `run()` retorna `TissueState` en lugar de `Dict[CloneId, Clone]`

### 5. Actualización de `CloneFactory` (`src/gillespie/clone_factory.py`)

- Pasa `config` al constructor de `Clone`.
- Configura los parámetros base de cada clon.
- Los parámetros de interacción global (`theta_I`, `beta`) se leen directamente desde `self.config` dentro de los clones, no se inyectan como atributos individuales.

## Uso

### Acceder al estado del tejido:
```python
# En TumorSimulation
state = self.tissue_state

# Obtener poblaciones por tipo
n_mutated = state.population_by_type("mutated")
n_immune = state.population_by_type("immune")
total = state.total_population()

# Obtener un clon específico
clone = state.get_clone(clone_id)
```

### En los métodos de tasa efectiva:
```python
# Un clon puede acceder a cualquier población del tejido
# y leer parámetros globales desde su config.
def death_rate_effective(self, tissue_state: TissueState) -> float:
    base_death = super().death_rate_effective(tissue_state)
    n_immune = tissue_state.population_by_type("immune")
    immune_killing = self.N * n_immune * self.config.theta_I
    return base_death + immune_killing
```

### Historial de estados:
```python
times, history, final_state = simulation.run()

# history es una lista de snapshots (Dict por CloneId)
# final_state es el TissueState final
for i, snapshot in enumerate(history):
    print(f"Time {times[i]}: {snapshot}")
```

## Ventajas

1. **Encapsulación**: El estado del sistema es explícito y accesible
2. **Flexibilidad**: Los clones pueden acceder a cualquier información de población
3. **Escalabilidad**: Fácil añadir nuevos tipos de clones o interacciones
4. **Historial limpio**: El método `snapshot()` mantiene un registro consistente
5. **Desacoplamiento**: La lógica de tasas está separada del simulador
6. **No hay dependencia circular en ejecución**: `TissueState` contiene los clones, pero los clones solo reciben `tissue_state` como parámetro en los cálculos de tasa.

## Posibles Mejoras Futuras

1. Mover el cálculo de crowding dentro del Clone usando `TissueState`
2. Crear estrategias de interacción como clases separadas
3. Añadir métodos para consultar interacciones específicas en `TissueState`
4. Implementar validaciones de integridad del estado
