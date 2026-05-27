from __future__ import annotations

from typing import Dict, Optional
from dataclasses import dataclass

from src.gillespie.cloneId import CloneId

from src.gillespie.clone_type import CloneType


@dataclass
class TissueState:
    """Encapsulates the complete state of the tissue at a given time.
    
    Contains all clones and their populations, and provides methods to query
    population counts by clone type.
    """
    
    clones: Dict[CloneId, 'Clone']  # type: ignore
    
    def total_population(self) -> int:
        """Get total number of cells across all clones."""
        return sum(clone.N for clone in self.clones.values())
    
    def population_by_type(self, cell_type: str) -> int:
        """Get total population for a specific cell type.
        
        Args:
            cell_type: The CloneType as string ('wild_type', 'mutated', 'immune', 'exhausted')
            
        Returns:
            Total number of cells of this type
        """
        return sum(
            clone.N for clone in self.clones.values() 
            if str(clone.cell_type) == cell_type
        )
    
    def get_clone(self, clone_id: CloneId) -> Optional['Clone']:  # type: ignore
        """Get a specific clone by its ID."""
        return self.clones.get(clone_id)
    
    def get_clones_by_type(self, cell_type: str) -> Dict[CloneId, 'Clone']:  # type: ignore
        """Get all clones of a specific type."""
        return {
            cid: clone for cid, clone in self.clones.items()
            if str(clone.cell_type) == cell_type
        }
    
    def snapshot(self) -> Dict[CloneId, dict]:
        """Create a snapshot of current state for history recording."""
        return {
            cid: {
                "N": clone.N,
                "rb": clone.birth_rate,
                "rd": clone.death_rate,
                "rm": clone.mutation_rate,
                "re": clone.exhaustion_rate,
            }
            for cid, clone in self.clones.items()
        }
