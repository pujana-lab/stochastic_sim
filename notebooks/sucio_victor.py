from typing import Dict, Type


class Clone:
    # El registro central y único
    _registry: Dict[str, Type["Clone"]] = {}

    def __init__(self, nombre: str):
        self.nombre = nombre
        # Cada vez que nace un clon, se anota a sí mismo (o a su clase) en la central
        Clone._registry[nombre] = type(self)


# 1. Creamos el primer clon
clon_1 = Clone("Alfa")
print(clon_1._registry)  # Resultado: {'Alfa': <class 'Clone'>}

# 2. Creamos un segundo clon más tarde
clon_2 = Clone("Beta")

# 3. ¿Qué ve el primer clon ahora?
print(clon_1._registry)  # Resultado: {'Alfa': <class 'Clone'>, 'Beta': <class 'Clone'>}