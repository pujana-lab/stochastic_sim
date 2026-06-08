# CloneType Refactoring: Problems & Solutions

## Executive Summary

The current `CloneType` enum creates unnecessary complexity and maintenance burden. The type information exists in two places (class hierarchy + enum), leading to confusion, redundant conversions, and scattered type checks throughout the codebase.

**Recommendation**: Replace `CloneType` enum with a simple `CloneRegistry` constant class that maps type names to Clone classes.

---

## Current Problems

### 1. **Multiple Representations of the Same Information**

Currently, a clone's type exists in three forms:
- **Class type**: `WildTypeClone`, `MutatedClone`, `ImmuneClone`, `ExhaustedClone`
- **Enum**: `CloneType.MUTATED`, `CloneType.IMMUNE`, etc.
- **String**: `"mutated"`, `"immune"`, etc.

**Example from code:**
```python
# clone_factory.py - Takes string
def create_clone(self, clone_type: str = "base", ...) -> Clone:

# clone.py - Sets enum attribute
self.cell_type: CloneType = CloneType.BASE

# tumor_simulation.py - Uses string to access enum value
clone_type = clone.cell_type
# ... later uses clone_type.value

# clone.py - Stores string for mutation
self.next_mutation: str = ""
```

**Problem**: Developers must convert between formats constantly. Easy to make mistakes like:
- Comparing strings to enums
- Using `.value` inconsistently
- Confusing which format is expected

---

### 2. **Redundant Type Storage**

Each Clone subclass IS a type (polymorphism), yet also explicitly stores its type:

```python
class WildTypeClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.BASE  # Redundant! Already is WildTypeClone

class ImmuneClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.IMMUNE  # Redundant!
```

**Problem**: 
- Adds noise to the code
- Risk of inconsistency (wrong cell_type assigned)
- Could use `type(clone).__name__` or isinstance() instead

---

### 3. **String-Based Type References Are Fragile**

Mutation types are stored as bare strings:

```python
# clone.py
self.next_mutation: str = ""

# clone_factory.py
clone.next_mutation = CloneType.MUTATED.value  # "mutated" as string

# tumor_simulation.py
self.create_clone(clone_id=..., clone_type=clone.next_mutation, ...)
```

**Problems**:
- No compile-time checking (typo = silent bug)
- Hard to search for all uses of a type ("mutated" appears as plain string)
- No IDE autocomplete or refactoring support
- Easy to add invalid type values

---

### 4. **Factory Method Type Parameter Is Loosely Typed**

```python
def create_clone(self, clone_type: str = "base", ...) -> Clone:
```

Accepts any string, no validation until inside the method. Invalid types create runtime errors.

---

### 5. **Scattered Type Checking Logic**

Type comparisons happen in different ways throughout:

```python
# tissue_state.py
pop_map.get(CloneType.BASE, 0)  # Enum as dict key
pop_map.get(CloneType.MUTATED, 0)

# tumor_simulation.py
if clone_type not in type_rates:
    type_rates[clone_type] = (...)  # clone_type is enum here
```

**Problem**: Inconsistent patterns make the code harder to read and maintain.

---

### 6. **TODO Comment in Code**

```python
# clone_type.py
#TODO: Esto es un lio que flipas, nunca se si llamarlo como value como string como CloneType.(tipo). 
# se puede hacer refactor de esto??
```

The team already acknowledges this is confusing!

---

## Proposed Solution: Two Approaches

### **RECOMMENDED: Use `__init_subclass__` for Auto-Registration**

This is the cleanest solution—**each Clone subclass self-registers** without needing a central registry:

```python
# src/gillespie/clone.py

from typing import Dict, Type

class Clone:
    """Base Clone class with auto-registration via __init_subclass__."""
    
    _registry: Dict[str, Type["Clone"]] = {}
    
    def __init_subclass__(cls, clone_type: str = None, **kwargs):
        """Auto-register each subclass with its type name."""
        super().__init_subclass__(**kwargs)
        if clone_type:
            Clone._registry[clone_type] = cls
    
    def __init__(self, clone_id: CloneId, config: SimulationConfig, N: int = 1, parent: Optional[CloneId] = None):
        self.clone_id: CloneId = clone_id
        self.config: SimulationConfig = config
        self.N = N
        self.parent = parent
        # ... rest of init
        self.next_mutation: str = ""  # Now stores string type name
        
    def get_type(self) -> str:
        """Get the type name of this clone."""
        # Find this class in the registry
        for type_name, cls in Clone._registry.items():
            if isinstance(self, cls) and type(self) == cls:
                return type_name
        raise ValueError(f"Clone type not registered: {type(self)}")

# Subclasses declare their type
class WildTypeClone(Clone, clone_type="base"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class MutatedClone(Clone, clone_type="mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class ImmuneClone(Clone, clone_type="immune"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

class ExhaustedClone(Clone, clone_type="exhausted"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
```

### **Key Benefits of `__init_subclass__` Approach**

1. ✅ **DRY (Don't Repeat Yourself)**: Adding a new type requires ONLY creating the subclass
2. ✅ **No central registry file needed**: Type info lives with the class
3. ✅ **Type-safe**: Constants (not strings) via class definition
4. ✅ **Automatic validation**: Registry is populated at import time
5. ✅ **Easy to extend**: New types don't require touching factory or registry
6. ✅ **Cleaner mutations**: `next_mutation = "mutated"` or use a helper method

### **Factory Becomes Simpler**

```python
# src/gillespie/clone_factory.py

def create_clone(self, clone_type: str = "base", N: int = 1, parent: Optional[Clone] = None) -> Clone:
    if clone_type not in Clone._registry:
        raise ValueError(f"Unknown clone type: {clone_type}")
    
    clone_class = Clone._registry[clone_type]
    clone = clone_class(
        clone_id=clone_id,
        config=self.config,
        N=N,
        parent=parent
    )
    
    # Setup type-specific attributes
    if clone_type == "base":
        clone.next_mutation = "mutated"
        clone.mutation_rate = self.config.nu0
        clone.K = self.config.K0
    elif clone_type == "mutated":
        clone.birth_rate = self.config.lambda0 * (1.0 + self.config.fitness_gain)
        clone.K = self.config.K_mutant
    # ... etc
    
    return clone
```

---

### **ALTERNATIVE: Use a Simple `CloneRegistry` Constant Class**

If your team prefers an explicit central location, keep a registry (but it's less elegant):

```python
# src/gillespie/clone_registry.py

class CloneRegistry:
    """Type names for clone types."""
    BASE = "base"
    MUTATED = "mutated"
    IMMUNE = "immune"
    EXHAUSTED = "exhausted"
    
    # Only needed if NOT using __init_subclass__
    ALL_TYPES = [BASE, MUTATED, IMMUNE, EXHAUSTED]
```

**Downside**: Still requires updating when adding new types (if not using auto-registration)

---

## Migration Steps (__init_subclass__ Approach)

### **Phase 1: Update Clone Base Class**
- [ ] Add `_registry: Dict[str, Type[Clone]]` class variable
- [ ] Implement `__init_subclass__(cls, clone_type: str = None)` method
- [ ] Add `get_type()` method to return type name
- [ ] Remove `self.cell_type: CloneType` attribute
- [ ] Remove `from src.gillespie.clone_type import CloneType` import

### **Phase 2: Update Clone Subclasses**
- [ ] Add `clone_type="base"` to `WildTypeClone(Clone, clone_type="base"):`
- [ ] Add `clone_type="mutated"` to `MutatedClone(Clone, clone_type="mutated"):`
- [ ] Add `clone_type="immune"` to `ImmuneClone(Clone, clone_type="immune"):`
- [ ] Add `clone_type="exhausted"` to `ExhaustedClone(Clone, clone_type="exhausted"):`
- [ ] Remove individual `self.cell_type = CloneType.X` assignments from `__init__`

### **Phase 3: Update CloneFactory**
- [ ] Use `Clone._registry[clone_type]` to get class
- [ ] Add validation: `if clone_type not in Clone._registry: raise ValueError(...)`
- [ ] Update all `clone.next_mutation` to use string type names (e.g., `"mutated"`)
- [ ] Remove references to `CloneType` enum

### **Phase 4: Update TissueState**
- [ ] Change `pop_map: Dict[CloneType, int]` → `Dict[str, int]`
- [ ] Replace `clone.cell_type` with `clone.get_type()` in all methods
- [ ] Update method signatures to use `str` instead of `CloneType`
- [ ] Remove import of `CloneType`
- [ ] Update `snapshot()` to use `clone.get_type()` instead of `clone.cell_type.value`

### **Phase 5: Update TumorSimulation & Utilities**
- [ ] Update `_build_rate_matrix()` to use `clone.get_type()` instead of `clone.cell_type`
- [ ] Update any code accessing `Event.clone_type` to work with strings
- [ ] Update `crowding_strategy.py` to work with string types
- [ ] Update any tests or utilities

### **Phase 6: Delete Old Code**
- [ ] Remove `src/gillespie/clone_type.py` (enum file)
- [ ] Remove all remaining `CloneType` imports
- [ ] Search codebase for `.value` conversions and remove unnecessary ones

---

## Impact Analysis (__init_subclass__ Approach)

### **Files Affected**

| File | Changes | Complexity |
|------|---------|-----------|
| `clone.py` | Add __init_subclass__, get_type(), remove cell_type | Low |
| `clone_factory.py` | Use Clone._registry, add validation | Low |
| `clone_type.py` | DELETE | - |
| `tissue_state.py` | Dict[str, int], replace clone.cell_type with clone.get_type() | Low |
| `tumor_simulation.py` | Use clone.get_type() instead of clone.cell_type | Low |
| `crowding_strategy.py` | Use string types instead of CloneType enum | Low |
| `test_*.py` | Update imports/assertions | Low |
| `event.py` | Update clone_type from CloneType to str | Low |

### **Breaking Changes**

- All references to `CloneType` enum removed
- `pop_map` changes from `Dict[CloneType, int]` → `Dict[str, int]`
- Method signatures change from `cell_type: CloneType` → `cell_type: str`

### **Adding New Types is Now Trivial**

No more changes needed across multiple files:
```python
# ONLY change needed: In clone.py
class HighlyMutatedClone(Clone, clone_type="highly_mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # Automatically registered! No other files to change!
```

---

## Example Refactoring (__init_subclass__)

### **Before:**
```python
# clone.py - Redundant type storage
class WildTypeClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.BASE  # Redundant!

class MutatedClone(Clone):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.cell_type = CloneType.MUTATED  # Redundant!

# clone_factory.py - String manipulation and enum values
def create_clone(self, clone_type: str = "base", ...) -> Clone:
    if clone_type == "mutated":
        clone = MutatedClone(...)
        clone.next_mutation = CloneType.MUTATED.value  # Convert enum to string!

# tissue_state.py - Enum as dict key
pop_map: Dict[CloneType, int] = {}
pop_map[clone.cell_type] = 10  # CloneType enum key

# tumor_simulation.py
rate_matrix.add_event(Event(..., clone_type=clone.cell_type))  # enum
if event.kind == EventType.MUTATION:
    self.create_clone(..., clone_type=clone.next_mutation)  # string
```

### **After (__init_subclass__):**
```python
# clone.py - Self-registering classes, no redundancy
class Clone:
    _registry: Dict[str, Type["Clone"]] = {}
    
    def __init_subclass__(cls, clone_type: str = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if clone_type:
            Clone._registry[clone_type] = cls
    
    def get_type(self) -> str:
        for type_name, cls in Clone._registry.items():
            if type(self) == cls:
                return type_name
        raise ValueError(f"Clone type not registered: {type(self)}")

class WildTypeClone(Clone, clone_type="base"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # No self.cell_type! Type is implicit in the class.

class MutatedClone(Clone, clone_type="mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # No self.cell_type!

# clone_factory.py - Simple, automatic validation
def create_clone(self, clone_type: str = "base", ...) -> Clone:
    if clone_type not in Clone._registry:
        raise ValueError(f"Unknown clone type: {clone_type}")
    
    clone_class = Clone._registry[clone_type]
    clone = clone_class(clone_id, self.config, N, parent)
    
    if clone_type == "mutated":
        clone.next_mutation = "mutated"  # Just a string, no enum conversion!
        clone.birth_rate = self.config.lambda0 * (1.0 + self.config.fitness_gain)

# tissue_state.py - String keys, cleaner
pop_map: Dict[str, int] = {}
pop_map[clone.get_type()] = 10  # String key, not enum

# tumor_simulation.py - Consistent string usage
rate_matrix.add_event(Event(..., clone_type=clone.get_type()))
if event.kind == EventType.MUTATION:
    self.create_clone(..., clone_type=clone.next_mutation)
```

### **Adding a New Type: Before vs After**

**Before** (with CloneType enum or CloneRegistry):
1. Create class in `clone.py`
2. Update `clone_factory.py` with creation logic
3. Add constant to registry/enum
4. Update any documentation

**After** (with __init_subclass__):
```python
# clone.py - That's it!
class HighlyMutatedClone(Clone, clone_type="highly_mutated"):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.birth_rate = self.config.lambda0 * 2.0  # Example setup
```
✅ Type is auto-registered. No other files need changes.

rate_matrix.add_event(Event(..., clone_type=clone.get_type()))
if event.kind == EventType.MUTATION:
    self.create_clone(..., clone_type=clone.next_mutation)
```

---

## Discussion Points for Team

### **On the __init_subclass__ Approach (RECOMMENDED)**

1. **Is the auto-registration mechanism clear enough?**
   - Pro: Very Pythonic, DRY, reduces maintenance burden
   - Con: Less explicit than a central registry file
   - Mitigation: Document the pattern well, add comments to `__init_subclass__`

2. **Should `get_type()` method be cached or computed each time?**
   - Current: Computed by iterating registry (minimal overhead)
   - Alternative: Cache in instance variable (tiny speedup, adds complexity)
   - Recommendation: Keep it simple, compute on demand

3. **What about type constants for easier reference?**
   - Current: Use strings directly ("mutated", "base")
   - Alternative: Add constants to Clone class for IDE autocomplete
   ```python
   class Clone:
       TYPE_BASE = "base"
       TYPE_MUTATED = "mutated"
       TYPE_IMMUNE = "immune"
       TYPE_EXHAUSTED = "exhausted"
   ```
   - Would allow: `clone.next_mutation = Clone.TYPE_MUTATED` (type-safe!)

4. **Should we validate that all required methods are implemented in subclasses?**
   - Current: Loose (Python's duck typing)
   - Alternative: Add abstract methods to Clone base class
   - Trade-off: More structure vs. less flexibility

### **Implementation Timeline**

- Phase 1-2 only change `clone.py` (lowest risk)
- Phase 3-5 are mechanical updates (straightforward)
- Phase 6 cleanup is final

Recommendation: **Do all phases in one go** or break into two PRs:
- PR1: clone.py changes (auto-registration logic)
- PR2: Rest of codebase (using the new system)

---

## Checklist for Implementation

- [ ] Team review & approval of approach
- [ ] Create `clone_registry.py`
- [ ] Update `clone.py` 
- [ ] Update `clone_factory.py`
- [ ] Update all imports across codebase
- [ ] Update/add tests
- [ ] Delete `clone_type.py`
- [ ] Run full test suite
- [ ] Update any documentation

---

## References

- Current TODO comment in `clone_type.py`
- Previous refactoring: TissueState encapsulation (in `tissue_state_refactor.md`)
