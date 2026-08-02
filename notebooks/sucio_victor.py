from typing import Dict, Type


class Clone:
    # The central, single registry
    _registry: Dict[str, Type["Clone"]] = {}

    def __init__(self, nombre: str):
        self.nombre = nombre
        # Every time a clone is created, it registers itself (or its class) in the central registry
        Clone._registry[nombre] = type(self)


# 1. Create the first clone
clon_1 = Clone("Alfa")
print(clon_1._registry)  # Result: {'Alfa': <class 'Clone'>}

# 2. Create a second clone later
clon_2 = Clone("Beta")

# 3. What does the first clone see now?
print(clon_1._registry)  # Result: {'Alfa': <class 'Clone'>, 'Beta': <class 'Clone'>}
