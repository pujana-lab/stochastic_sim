
from src.gillespie.tissue_state import TissueState

# esta funcion es una prueba de concepto. si queremos que los clones de cada tipo se comporten de la misma forma deberiamos solo calcular los parametros una vez por tipo.

# Ahor mismo estamos calculando de nuevo los mismos valores para cada clon. he empezado a intentar implementarlo 
def create_parameter_map_by_type(tissue_state: TissueState) -> dict[str, dict]:
    """
    Escanea el tejido, encuentra todos los tipos de células únicos actuales
    y genera un diccionario base mapeando cada tipo a sus nuevos parámetros.
    """

    unique_cell_types = {clone.get_type() for clone in tissue_state.clones.values()}


    parameter_map = {
        cell_type: {
            "birth_rate": None,   
            "death_rate": None,
            "custom_multiplier": 1.0
        }
        for cell_type in unique_cell_types
    }
    
    return parameter_map