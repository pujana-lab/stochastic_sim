from __future__ import annotations

from typing import Dict, Optional
from dataclasses import dataclass, field

from src.gillespie.cloneId import CloneId
from src.gillespie.clone import Clone
from src.gillespie.clone_type import CloneType


@dataclass
class TissueState:
    """Encapsulates the complete state of the tissue at a given time.
    
    Contains all clones and their populations, and provides methods to query
    population counts by clone type.
    """
    ## igual es mejor guardar como array de clones 
    clones: Dict[CloneId, 'Clone']  # type: ignore
    pop_map: Dict[CloneType, int] = field(init=False)

    def __post_init__(self) -> None:
        self.update_pop_map()

    def update_pop_map(self) -> None:
        self.pop_map = self.get_pop_map()

    def total_population(self) -> int:
        """Get total number of cells across all clones."""
        return sum(clone.N for clone in self.clones.values())
    
    # REPLACED BY GET_POP_MAP
    # --------------------------------------------------------------------------------
    def population_by_type(self, cell_type: CloneType) -> int:
        """Get total population for a specific cell type.
        
        Args:
            cell_type: The CloneType as string ('base', 'mutated', 'immune', 'exhausted')
            
        Returns:
            Total number of cells of this type
        """
        return sum(
            clone.N for clone in self.clones.values() 
            if clone.cell_type == cell_type
        )
    # --------------------------------------------------------------------------------
    
    
    def get_pop_map(self, clone_type: CloneType = None) -> Dict[CloneType,int]:
        
        pop_map: Dict[CloneType, int]={}
        
        for clone in self.clones.values():
            if clone_type is not None and clone.cell_type != clone_type:
                continue
            pop_map[clone.cell_type] = pop_map.get(clone.cell_type,0) + clone.N
        
        return pop_map

    
    def get_clone(self, clone_id: CloneId) -> Optional['Clone']:  # type: ignore
        """Get a specific clone by its ID."""
        return self.clones.get(clone_id)
    
    def get_clones_by_type(self, cell_type: CloneType) -> Dict[CloneId, 'Clone']:  # type: ignore
        """Get all clones of a specific type."""
        return {
            cid: clone for cid, clone in self.clones.items()
            if clone.cell_type == cell_type
        }
    
    def snapshot(self) -> Dict[CloneId, dict]:
        """Create a snapshot of current state for history recording."""
        return {
            cid: {
                "Type": clone.cell_type.value,
                "N": clone.N,
                "rb": clone.birth_rate,
                "rd": clone.death_rate,
                "rm": clone.mutation_rate,
                "re": clone.exhaustion_rate,
            }
            for cid, clone in self.clones.items()
        }
