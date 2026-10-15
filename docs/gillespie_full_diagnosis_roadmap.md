# Diagnóstico completo y roadmap para corregir la implementación Gillespie / ODE

> **Historical diagnosis.** The claims and proposed fixes below describe an
> earlier investigation. They are not a current specification; use
> [`ARCHITECTURE.md`](ARCHITECTURE.md) and the source/tests as the authority.

## Resumen ejecutivo

Hay un problema estructural en la implementación actual: el modelo estocástico y el modelo determinista no están representando la misma dinámica. La diferencia no es solo ruido estocástico; hay un desajuste conceptual y varios bugs de cálculo que hacen que la simulación Gillespie no sea comparable con el sistema ODE que se usa para visualizarla.

La causa principal es que la matriz de tasas se está construyendo por clon individual, mientras que el modelo agregado del ODE se expresa por tipo de población. Además, hay varios puntos en el código donde la lógica es inconsistente, especialmente en la construcción de tasas, la logística de crowding y la comparación con el sistema determinista.

El objetivo de este documento es:

1. diagnosticar la causa raíz de la diferencia,
2. separar lo que es ruido estocástico de lo que es un bug real,
3. proponer un roadmap claro para arreglar la implementación sin cambiar la intención del modelo.

---

## 1. Diagnóstico técnico

### 1.1. Error principal: la matriz de tasas se calcula por clon y por tipo a la vez

El problema más grave aparece en `src/gillespie/tumor_simulation.py` dentro del método `_build_rate_matrix()`.

El código recorre `self.tissue_state.clones.items()`, y para cada clon vivo hace esto:

- calcula `birth_rate_effective`, `death_rate_effective`, etc.,
- crea un evento por tipo de reacción para ese clon,
- añade cada evento a la matriz de tasas.

Esto significa que, si hay varios clones del mismo tipo, cada uno genera su propia contribución a la tasa de nacimiento, muerte, mutación o agotamiento. Eso está mezclando dos niveles de agregación:

- por clon individual, y
- por tipo de población.

Esto es inconsistente con la idea del Gillespie SSA clásico, donde los eventos se definen sobre reacciones globales y los reactantes se cuentan por población total o por grupos de reacción correctamente definidos.

La propia nota en el código lo indica:

- `#BUG: aqui estamos mezclando dos logicas...`
- `#TODO: hay que volver a poner los rates por TIPO...`

Es decir, el problema ya estaba identificado por el autor, pero no se corrigió.

### Por qué esto afecta a la comparación con ODE

El modelo determinista usa variables agregadas por tipo:

- `C` (células base / WT)
- `T` (inmunes / effectors)
- `E` (exhausted)

Mientras que el Gillespie actual modela cada clon individualmente como una entidad separada con su propio `clone_id`. Eso produce:

- eventos duplicados por poblaciones equivalentes,
- doble conteo cuando existen varios clones del mismo tipo,
- cambios de tasa no equivalentes a la dinámica continua agregada,
- patrones que no corresponden a la misma ecuación diferencial.

### Resultado

Lo que se ve como “diferencia determinista vs. simulado” no es solo un efecto estocástico. Hay un error de especificación del modelo de reacción.

---

### 1.2. El ODE de comparación no usa el mismo modelo que la simulación Gillespie

La figura / script de comparación está en `notebooks/visualization/plotear_populations.py`.

Allí el sistema ODE se define como:

```python
def system_ode(t, x, p):
    C, T, E = x
    dC = C * (1.0 - C) - p["gamma"] * C - p["sigma"] * C * T
    dT = p["rho"] * T * (1.0 - T - E) + (p["beta"] - p["omega"]) * C * T
    dE = p["omega"] * C * T - p["mu"] * E
    return [dC, dT, dE]
```

Esto es un sistema adimensional diferente al que describe la lógica de `Clone` / `TissueState` / `SimulationConfig`:

- las tasas del Gillespie se construyen con `lambda0`, `mu`, `nu`, `omega_exhaust`, `theta_I`, `beta`, crowding, etc.,
- las ecuaciones del ODE usan parámetros con nombres distintos (`rho`, `sigma`, `gamma`, `omega`, `mu`, `beta`),
- la escala del sistema no coincide necesariamente con la de la simulación,
- el esquema de población también cambia: `C/T/E` no están ligadas 1:1 a `base/mutated/immune/exhausted` en la implementación de clones.

En otras palabras: la comparación es entre dos modelos distintos, no entre una simulación y su ODE equivalente.

---

### 1.3. El crowding / logística está mal alineado con la dinámica de reacción

La estrategia de crowding está en `src/gillespie/crowding_strategy.py`.

La función `crowding()` devuelve:

```python
return max(0.0, 1.0 - numerator_N / Kt)
```

Esto significa que la reducción de crecimiento es proporcional a la densidad relativa del conjunto competitivo. Eso puede ser razonable, pero hay varios problemas:

1. `numerator_N` depende del tipo de clone y del estado actual de `tissue_state`,
2. `Kt` se calcula con `clone.actual_K`, que depende de `birth_rate`, `death_rate`, etc.,
3. `clone.crowding_numerator` y `crowding_strategy` usan lógica no siempre coherente con la agregación por tipo,
4. el método `AdaptedCrowding.calculate_K` puede producir capacidades que no están bien acopladas al modelo de reacción en tiempo continuo.

Esto puede crear:

- crecimiento saturado demasiado rápido,
- desbalance entre tipos,
- tasas de natalidad artificialmente reducidas o exageradas,
- variabilidad no comparable con ODE.

---

### 1.4. Hay una posible doble contabilización de eventos por clon

Este punto se ve claramente en `Clone.birth_rate_effective()`:

```python
return self.birth_rate * self.N * self.config.crowding_strategy.crowding(self,tissue_state=tissue_state)
```

y en el comentario:

```python
#BUG: related to bug in tumor_simulation.py ahora mismo el birth rate va por tipo y no por clon. lo cual esta double counting
```

Es decir, el crecimiento de un tipo puede estar siendo contado dos veces: una vez por la lógica del clon y otra por la agregación del tipo o de la interacción. El resultado es un sesgo hacia tasas mayores o hacia una mutación/crecimiento excesivo en la trayectoria simulada.

---

### 1.5. El mecanismo de mutación y agotamiento tiene una semántica débil

En `TumorSimulation._introduce_mutation()` y `_induce_exhaustion()` vemos que la mutación y la agotamiento se modelan como cambios de clon y tamaño poblacional, pero la matriz de tasas se construye sin respetar consistentemente el conjunto de reacciones globales.

Puntos problemáticos:

- `_introduce_mutation()` crea un clon hijo con `N=2` y mata al clon padre,
- el modelo de agotamiento hace `self.tissue_state.clones[(-2,)].divide()`, lo que es una operación de una “clase global” de exhaustos,
- el contador de mutación / agotamiento depende de identidad de clone, no de una reacción por tipo,
- la `RateMatrix` se elige por `u = rng.random()` sobre eventos ya construidos, pero los eventos tienen tasas no agregadas en el mismo espacio de estado.

Esto genera transiciones que pueden producir trayectorias muy distintas a la evolución continua de un ODE agregado.

---

### 1.6. Hay una discrepancia entre parámetros físicos y escala del sistema

`SimulationConfig` hace scaling con `OMEGA` y otras cantidades, pero la comparación ODE usa valores adimensionales en un rango tipo `[0,1]`. El modelo de Gillespie usa pop_map y K de tipo natural, con escalas muy distintas a la ODE.

Esto provoca que, incluso si ambos modelos fueran correctos, no estén intentando resolver el mismo sistema en la misma escala. La comparación directa sin normalización es engañosa.

---

## 2. Qué es ruido estocástico y qué es bug verdadero

### Lo que es ruido estocástico

La diferencia esperada entre Gillespie y ODE puede venir de:

- aleatoriedad de la trayectoria,
- número finito de células,
- eventos discretos en un sistema de pequeña población,
- fluctuaciones alrededor de la media.

Eso es completamente normal.

### Lo que es bug verdadero

Los siguientes elementos son bugs o inconsistencias de modelado:

1. `RateMatrix` construido por clon individual y no por reacción agregada.
2. doble conteo de tasas por clon/tipo.
3. comparación entre un sistema ODE no equivalente y una simulación Gillespie diferente.
4. crowding / logistic no alineado con la agregación del sistema.
5. falta de una definición clara del “estado macro” para comparar con ODE.

---

## 3. Diagnóstico estructural de la implementación

### 3.1. La arquitectura actual mezcla tres niveles diferentes

La implementación actual mezcla en el mismo código:

- nivel de clone individual,
- nivel de tipo de población,
- nivel de ecuación diferencial agregada.

Eso hace que sea imposible comparar correctamente la simulación con el ODE sin una etapa de agregación explícita.

### 3.2. El sistema no tiene un “esquema de reacción centralizado”

En SSA, lo recomendable es definir de forma explícita:

- reacciones por tipo o por clase,
- tasas de reacción globales,
- cambio de estado por evento,
- actualización del estado total.

La implementación actual no separa claramente esos conceptos. La matriz de tasas se construye ad hoc desde los clones vivos, en vez de desde unas reacciones bien definidas.

### 3.3. La documentación del modelo no coincide con la implementación

El README describe una estructura conceptual clara (base, mutated, immune, exhausted), pero la lógica real del código no siempre cumple ese nivel de abstracción. El modelo está más “mecánico” que “específico”, y eso hace más difícil mantener la coherencia entre ODE y Gillespie.

---

## 4. Roadmap de corrección

### Fase 0: definir el modelo objetivo

Antes de tocar el código, hay que fijar una meta clara:

- ¿Queremos un Gillespie SSA para las poblaciones agregadas por tipo?
- ¿Queremos un Gillespie SSA por clon individual?
- ¿Queremos una comparación con una ODE agregada exacta del mismo modelo?

La recomendación es la siguiente:

- usar Gillespie por tipo agregado para facilitar la comparación con ODE,
- o mantener Gillespie por clon individual si se quiere un modelo más detallado, pero entonces la ODE debe derivarse del mismo sistema y no de una versión “aproximada”.

### Fase 1: separar los niveles de abstracción

#### Objetivo

Crear una separación limpia entre:

- `Clone`: estado individual
- `TissueState`: estado agregado del tejido
- `RateMatrix`: reacciones globales
- `SimulationConfig`: parámetros del modelo

#### Tareas

- definir explícitamente qué eventos existen,
- definir qué variables agregadas importan al modelo,
- eliminar la construcción de tasas dentro del bucle por clon si el modelo objetivo es agregado.

---

### Fase 2: corregir la matriz de tasas

#### Objetivo

Asegurar que cada reacción se define una sola vez y sobre el mismo espacio de estado.

#### Reglas a cumplir

- cada reacción debe tener una tasa bien definida por tipo o por clase,
- si hay varios clones del mismo tipo, su tasa debe agregarse antes de crear la reacción,
- la probabilidad de cada evento debe computarse sobre el total de tasa relevante,
- `choose_event()` debe seleccionar sobre un conjunto consistente de eventos.

#### En la práctica

Reemplazar esta lógica:

- iterar sobre clones,
- crear eventos por clone,

por esta otra:

- iterar sobre tipos relevantes,
- calcular la tasa agregada por tipo,
- añadir un solo evento por tipo y reacción.

Ejemplo conceptual:

```python
for type_name in ["base", "mutated", "immune", "exhausted"]:
    N = tissue_state.pop_map.get(type_name, 0)
    birth_rate = ...
    death_rate = ...
    mutation_rate = ...
    exhaustion_rate = ...
    add_event(type_name, reaction, total_rate)
```

---

### Fase 3: corregir la lógica de crowding y logística

#### Objetivo

Hacer que la forma logística coincida con el modelo de reacción.

#### Tareas

- decidir si el crowding se aplica a la tasa de nacimiento global o a la probabilidad de división por individuo,
- mantener una única fórmula para `C(t)` y aplicarla de forma consistente,
- revisar `SimpleCrowding` y `AdaptedCrowding` para que sus resultados sean comparables con la ODE equivalente,
- evitar que `Kt` se calcule con factores inconsistentes.

---

### Fase 4: unificar el modelo determinista con el modelo estocástico

#### Objetivo

Al crear el ODE, derivarlo del mismo conjunto de reacciones que el sistema Gillespie.

#### Reglas

- si el Gillespie usa tasa agregada por tipo, el ODE debe usar exactamente esas mismas tasas promedio,
- las ecuaciones diferenciales deben obtenerse por balance de reacciones,
- no se debe comparar Gillespie por clon con ODE por tipo sin una agregación explícita.

#### Resultado esperado

Diferencia esperada: solo ruido estocástico alrededor de la trayectoria media, no un desplazamiento sistemático por bug.

---

### Fase 5: pruebas de regresión

#### Objetivo

Crear un conjunto mínimo de tests para proteger el modelo.

#### Casos a cubrir

1. una sola población sin interacciones,
2. dos poblaciones con crecimiento y muerte,
3. mutación mutada → nueva clase,
4. agotamiento de inmunes,
5. comparación ODE vs Gillespie en régimen de gran población,
6. estrategia de crowding con saturación.

#### Recomendación

No dejar la validación solo en la visualización. Debe existir un test que verifique que la tasa total y la suma de reacciones coinciden con la dinámica esperada.

---

### Fase 6: limpieza final y frontera del modelo

#### Objetivo

Dejar la implementación mantenible y legible.

#### Tareas

- eliminar comentarios `BUG`/`TODO` que señalan decisiones inconsistentes,
- documentar formalmente las reacciones del modelo,
- separar las funciones matemáticas de la lógica del simulador,
- mantener una única fuente de verdad para los parámetros del sistema.

---

## 5. Prioridad de corrección

### Prioridad alta

1. arreglar la construcción de `RateMatrix`,
2. definir una agregación por tipo consistente,
3. unificar el Gillespie con la ODE real del modelo.

### Prioridad media

4. revisar crowding,
5. limpiar mutación y agotamiento,
6. fortalecer pruebas.

### Prioridad baja

7. refactor de estilo y documentación,
8. limpieza de comentarios antiguos,
9. optimización secundaria.

---

## 6. Recomendación final

La implementación actual tiene un problema peor que un simple desajuste numérico: el modelo estocástico no está definido sobre la misma variable de estado que la aproximación determinista. Por eso la diferencia observada no es una evidencia de que el Gillespie “se comporte mal”, sino de que el sistema está comparándose con un modelo distinto.

La corrección adecuada no es “ajustar parámetros para que se parezcan”, sino:

- definir un único modelo de reacción,
- agotar la lógica por tipo/estado,
- derivar la ODE del mismo conjunto de reacciones,
- validar con pruebas reales.

---

## 7. Checklist de implementación recomendada

- [ ] definir el espacio de estados del modelo,
- [ ] decidir si el SSA es por tipo agregado o por clon individual,
- [ ] corregir `RateMatrix` para evitar doble conteo,
- [ ] revisar `Clone.birth_rate_effective` y `death_rate_effective`,
- [ ] unificar `crowding` con la tasa global,
- [ ] derivar un ODE equivalente al mismo sistema,
- [ ] añadir pruebas de regresión,
- [ ] comparar en régimen grande para validar equivalencia asintótica.

---

## 8. Conclusión corta

La diferencia observada entre determinista y simulado no es solo ruido: hay un bug estructural en la lógica de tasas y una comparación entre modelos no equivalentes. La solución correcta pasa por unificar la representación del sistema y hacer que Gillespie y ODE operen sobre el mismo balance de reacciones.
