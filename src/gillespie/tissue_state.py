from __future__ import annotations

from typing import Dict, Optional
from dataclasses import dataclass, field

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone


@dataclass
class TissueState:
    """Encapsulates the complete state of the tissue at a given time.
    
    Contains all clones and their populations, and provides methods to query
    population counts by clone type.
    """
    ## igual es mejor guardar como array de clones 
    clones: Dict[CloneId, 'Clone']  # type: ignore
    pop_map: Dict[str, int] = field(init=False)

    def __post_init__(self) -> None:
        self.update_pop_map()

    def update_pop_map(self) -> None:
        self.pop_map = self.get_pop_map()

    def print_pop_map(self) -> None:
        """Print the current pop_map as a simple table."""
        self.update_pop_map()
        header_type = "clone_type"
        header_count = "count"
        max_type = max((len(ct) for ct in self.pop_map), default=len(header_type))
        print(f"{header_type:<{max_type}} | {header_count}")
        print(f"{'-' * max_type}-+-{'-' * len(header_count)}")
        for clone_type, count in sorted(self.pop_map.items()):
            print(f"{clone_type:<{max_type}} | {count}")

    def total_population(self) -> int:
        """Get total number of cells across all clones."""
        return sum(clone.N for clone in self.clones.values())
    
    # REPLACED BY GET_POP_MAP
    # --------------------------------------------------------------------------------
    def population_by_type(self, cell_type: str) -> int:
        """Get total population for a specific cell type.
        
        Args:
            cell_type: The clone type as string ('base', 'mutated', 'immune', 'exhausted')
            
        Returns:
            Total number of cells of this type
        """
        return sum(
            clone.N for clone in self.clones.values() 
            if clone.get_type() == cell_type
        )
    # --------------------------------------------------------------------------------
    
    
    def get_pop_map(self, clone_type: str = None) -> Dict[str, int]:
        
        pop_map: Dict[str, int]={}
        
        for clone in self.clones.values():
            clone_type_str = clone.get_type()
            if clone_type is not None and clone_type_str != clone_type:
                continue
            pop_map[clone_type_str] = pop_map.get(clone_type_str,0) + clone.N
        
        return pop_map

    
    def get_clone(self, clone_id: CloneId) -> Optional['Clone']:  # type: ignore
        """Get a specific clone by its ID."""
        return self.clones.get(clone_id)
    
    def get_clones_by_type(self, cell_type: str) -> Dict[CloneId, 'Clone']:  # type: ignore
        """Get all clones of a specific type."""
        return {
            cid: clone for cid, clone in self.clones.items()
            if clone.get_type() == cell_type
        }
    
    def snapshot(self) -> Dict[CloneId, dict]:
        """Create a snapshot of current state for history recording."""
        return {
            cid: {
                "Type": clone.get_type(),
                "N": clone.N,
                "rb": clone.birth_rate,
                "rd": clone.death_rate,
                "rm": clone.mutation_rate,
                "re": clone.exhaustion_rate,
            }
            for cid, clone in self.clones.items()
        }
