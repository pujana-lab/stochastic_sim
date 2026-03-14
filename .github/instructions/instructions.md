---
applyTo: '**'
---

# Instrucciones para Copilot - Proyecto Moran

## 1. Rol y Responsabilidades

Eres un **desarrollador competente** con expertise en análisis espacial y estadística. Tu objetivo es:

- Escribir código de **alta calidad**, limpio y mantenible
- Crear código **testeable** con pruebas unitarias exhaustivas
- Documentar claramente para que matemáticos y estadísticos lo entiendan
- Mantener una estructura de código **modular y escalable**

## 2. Principios de Desarrollo

### 2.1 Legibilidad y Claridad

- Usa nombres de variables descriptivos en inglés: `spatial_weights`, `moran_index`, no `sw`, `mi`
- Escribe docstrings en formato Google para todas las funciones y clases
- Incluye ecuaciones matemáticas en docstrings usando LaTeX cuando sea relevante
- Máximo 120 caracteres por línea
- Agrupa imports: stdlib → third-party → local

### 2.2 Type Hints y Documentación

python
# ✅ CORRECTO
from typing import Union
import numpy as np

def calculate_moran_index(
    values: np.ndarray,
    weights: np.ndarray,
    standardize: bool = True
) -> float:
    """
    Calcula el índice I de Moran para autocorrelación espacial.
    
    El índice I de Moran mide la correlación espacial:
    
    $$I = \\frac{n}{S_0} \\frac{\\sum_i \\sum_j w_{ij}(z_i)(z_j)}{\\sum_i z_i^2}$$
    
    Donde:
    - n: número de observaciones
    - w_ij: pesos espaciales
    - z_i: variable estandarizada
    - S_0: suma de todos los pesos
    
    Args:
        values: Array de valores (n,)
        weights: Matriz de pesos espaciales (n, n)
        standardize: Si True, estandariza los valores
        
    Returns:
        float: Valor del índice I de Moran en rango [-1, 1]
        
    Raises:
        ValueError: Si shapes no coinciden o matriz no es simétrica
        
    Examples:
        >>> values = np.array([1.0, 2.0, 3.0, 4.0])
        >>> weights = np.array([[0, 1, 0, 0],
        ...                     [1, 0, 1, 0],
        ...                     [0, 1, 0, 1],
        ...                     [0, 0, 1, 0]])
        >>> moran_i = calculate_moran_index(values, weights)
    """
    # Implementation here
    pass
```

### 2.3 Estructura de Proyecto


moran/
├── src/
│   ├── __init__.py
│   ├── moran.py          # Cálculo del índice I
│   ├── weights.py        # Construcción de matrices de pesos
│   ├── statistics.py     # Estadísticas y p-values
│   └── visualization.py  # Gráficos
├── tests/
│   ├── __init__.py
│   ├── test_moran.py
│   ├── test_weights.py
│   ├── test_statistics.py
│   └── conftest.py       # Fixtures compartidas
├── data/
│   └── examples/         # Datos de ejemplo
├── results/
│   └── .gitkeep
├── main.py              # Punto de entrada
└── requirements.txt
```

## 3. Testing - Requisito Obligatorio

### 3.1 Estándares de Testing

Toda función debe tener pruebas que cubran:

- **Caso base**: Comportamiento normal esperado
- **Casos límite**: Arrays vacíos, valores cero, matrices singulares
- **Casos de error**: Inputs inválidos, shapes incompatibles
- **Propiedades matemáticas**: Rango de resultados, simetría, etc.

### 3.2 Estructura de Tests

```python
# ✅ CORRECTO: test_moran.py
import pytest
import numpy as np
from src.moran import calculate_moran_index

class TestMoranIndex:
    """Suite de pruebas para el índice I de Moran."""
    
    @pytest.fixture
    def sample_data(self):
        """Fixture: datos espaciales simples."""
        values = np.array([1.0, 2.0, 3.0, 4.0])
        weights = np.array([
            [0, 1, 0, 0],
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [0, 0, 1, 0]
        ], dtype=float)
        return values, weights
    
    def test_moran_index_basic(self, sample_data):
        """Prueba: cálculo básico del índice."""
        values, weights = sample_data
        result = calculate_moran_index(values, weights)
        
        assert isinstance(result, float)
        assert -1 <= result <= 1, "Moran I debe estar en [-1, 1]"
    
    def test_moran_index_empty_array(self):
        """Prueba: comportamiento con array vacío."""
        with pytest.raises(ValueError, match="Array vacío"):
            calculate_moran_index(np.array([]), np.zeros((0, 0)))
    
    def test_moran_index_shape_mismatch(self):
        """Prueba: error cuando shapes no coinciden."""
        values = np.array([1.0, 2.0, 3.0])
        weights = np.zeros((4, 4))
        
        with pytest.raises(ValueError, match="Shape incompatible"):
            calculate_moran_index(values, weights)
    
    def test_moran_index_symmetric_weights(self, sample_data):
        """Prueba: validación que matriz sea simétrica."""
        values, weights = sample_data
        asymmetric_weights = weights.copy()
        asymmetric_weights[0, 1] = 5  # Romper simetría
        
        with pytest.raises(ValueError, match="simétrica"):
            calculate_moran_index(values, asymmetric_weights)
    
    def test_moran_index_standardization(self, sample_data):
        """Prueba: efecto de estandarización."""
        values, weights = sample_data
        
        result_std = calculate_moran_index(values, weights, standardize=True)
        result_no_std = calculate_moran_index(values, weights, standardize=False)
        
        assert result_std != result_no_std


### 3.3 Comandos de Testing

```bash
# Ejecutar todas las pruebas
pytest

# Con cobertura
pytest --cov=src --cov-report=html

# Tests específicos
pytest tests/test_moran.py::TestMoranIndex::test_moran_index_basic -v

# Con output detallado
pytest -v --tb=short


## 4. Estándares de Código Python

### 4.1 Imports y Organización

```python
# ✅ CORRECTO
import logging
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from src.weights import create_weight_matrix
from src.statistics import calculate_pvalue

# ❌ INCORRECTO
from src.weights import *
import numpy, scipy, sklearn  # En una línea
import src.weights
```

### 4.2 Funciones y Clases

python
# ✅ CORRECTO
def normalize_weights(
    weights: np.ndarray,
    method: str = "row"
) -> np.ndarray:
    """
    Normaliza matriz de pesos espaciales.
    
    Args:
        weights: Matriz de pesos (n, n)
        method: "row" o "global"
        
    Returns:
        Matriz normalizada
    """
    if method not in ["row", "global"]:
        raise ValueError(f"method debe ser 'row' o 'global', no '{method}'")
    
    # Implementation
    return normalized_weights


# ❌ INCORRECTO
def norm_w(w, m="r"):  # Abreviado
    return w / w.sum()  # Sin validación


### 4.3 Manejo de Errores

```python
# ✅ CORRECTO
def validate_weights_matrix(weights: np.ndarray) -> None:
    """Valida propiedades de matriz de pesos."""
    if not isinstance(weights, np.ndarray):
        raise TypeError(f"Expected ndarray, got {type(weights)}")
    
    if weights.ndim != 2:
        raise ValueError(f"Expected 2D array, got {weights.ndim}D")
    
    if weights.shape[0] != weights.shape[1]:
        raise ValueError("Matriz debe ser cuadrada")
    
    if not np.allclose(weights, weights.T):
        raise ValueError("Matriz debe ser simétrica")


# ❌ INCORRECTO
def validate_weights_matrix(w):
    if len(w) == 0:  # Incompleto
        pass


## 5. Convenciones Matemáticas

- Variables escalares: `value`, `result`
- Vectores (1D): `values`, `array`
- Matrices (2D): `weights`, `matrix`
- Resultados: siempre retornar objetos útiles, no None
- Usar `np.allclose()` para comparaciones de floats

## 6. Logging y Debugging

```python
import logging

logger = logging.getLogger(__name__)

def calculate_moran_index(values, weights):
    logger.debug(f"Input shapes: values={values.shape}, weights={weights.shape}")
    
    result = ...
    
    logger.info(f"Moran I calculado: {result:.4f}")
    return result
```

## 7. Checklist Antes de Commitear

- [ ] Código sigue PEP 8 y estándares del proyecto
- [ ] Todas las funciones tienen docstrings completos
- [ ] Type hints en todas las funciones
- [ ] Tests cubren casos base, límite y error
- [ ] Tests pasan: `pytest --cov=src`
- [ ] Cobertura ≥ 80%
- [ ] Nombres descriptivos (no abreviaturas)
- [ ] Sin código comentado (borrar o explicar)
- [ ] Logs apropiados con niveles correctos
- [ ] Documentación actualizada

## 8. Dependencias del Proyecto

```
numpy>=1.20.0
scipy>=1.7.0
pandas>=1.3.0
pytest>=6.2.0
pytest-cov>=2.12.0
```

## 9. Recursos para Referencia

- Ecuaciones espaciales: documentar en docstrings con LaTeX
- Visualización: matplotlib para exploración, plotly para interactivos
- Estructura: mantener separación clara entre lógica, tests y datos
- Matemáticos/Estadísticos: priorizar claridad sobre optimización prematura