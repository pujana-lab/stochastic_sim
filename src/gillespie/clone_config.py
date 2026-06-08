from dataclasses import dataclass

@dataclass
class CloneConfig:
    N: int = 0
    birth_rate:float = 0.0
    death_rate:float = 0.0
    mutation_rate:float = 0.0
    exhaustion_rate:float = 0.0

    K: int = 0
    next_mutation: str =  None
