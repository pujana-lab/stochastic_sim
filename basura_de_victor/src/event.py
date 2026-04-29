
from dataclasses import dataclass
from src.cloneId import CloneId

@dataclass
class Event:
    kind: str
    clone_id: CloneId
    rate: float
    reaction_number: int

