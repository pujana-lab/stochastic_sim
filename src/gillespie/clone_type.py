from enum import Enum
#TODO: Esto es un lio que flipas, nunca se si llamarlo como value como string como CloneType.(tipo). se puede hacer refactor de esto??
class CloneType(Enum):
    IMMUNE = "immune"
    BASE = "base"
    MUTATED = "mutated"
    EXHAUSTED = "exhausted"


    