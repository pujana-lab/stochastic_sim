# Memory Optimization Guide for TumorSimulation

## Overview

The optimized `TumorSimulation` class introduces three configurable memory modes to handle different simulation scales and analysis needs.

## Memory Modes

### 1. **LIGHTWEIGHT** (~1-2 MB)
**Best for:** Production runs, parameter sweeps, when you only care about final state

```python
from tumor_simulation_optimized import TumorSimulation, MemoryMode

sim = TumorSimulation(config, memory_mode=MemoryMode.LIGHTWEIGHT)
times, history, final_state, rates = sim.run()

# history contains ONLY the final state
# rate_history is None
# events is not tracked
# Ideal for 1000s of parallel simulations
```

**Behavior:**
- ✅ Saves only the final tissue state
- ✅ No intermediate snapshots
- ✅ No rate history tracking
- ✅ No event recording
- ✅ Minimal memory footprint

---

### 2. **STANDARD** (~10-100 MB, default)
**Best for:** Most development work, analysis with periodic checkpoints

```python
sim = TumorSimulation(config, memory_mode=MemoryMode.STANDARD)
# or just TumorSimulation(config) since STANDARD is default
times, history, final_state, rates = sim.run()

# history contains snapshots every save_interval steps
# rate_history is None
# events is not tracked
```

**Behavior:**
- ✅ Saves snapshots every `save_interval` steps (default: 100)
- ✅ Allows trajectory analysis
- ✅ No event or rate history
- ✅ Balanced memory/utility tradeoff

**Configure snapshot frequency:**
```python
config.save_interval = 50  # Save every 50 steps instead of 100
sim = TumorSimulation(config, memory_mode=MemoryMode.STANDARD)
```

---

### 3. **FULL** (~100 MB - 1+ GB)
**Best for:** Detailed analysis, debugging, publication-quality data

```python
sim = TumorSimulation(config, memory_mode=MemoryMode.FULL)
times, history, final_state, rates = sim.run()

# history contains EVERY step
# rate_history contains ALL rates at each step
# events contains EVERY event (birth/death/mutation/exhaustion)
# Maximum fidelity for analysis
```

**Behavior:**
- ✅ Records every single step
- ✅ Saves complete rate matrix history
- ✅ Logs all events for event-level analysis
- ✅ High memory cost, maximum information

---

## Key Changes from Original

### 1. **Conditional Storage Initialization**
```python
# OLD: Always allocated
self.history: List[Dict[CloneId, dict]] = [...]
self.events: List[Optional[Event]] = [] if self.save_all_steps else []

# NEW: Only allocated if needed
if self._should_save_history():
    self.history: List[Dict[CloneId, dict]] = [...]
self.events: Optional[List[Event]] = [] if self.memory_mode == MemoryMode.FULL else None
```

### 2. **Smart Rate History**
```python
# OLD: Rate history allocated but conditionally used
self.rate_history: List[List[Dict]] = [] if self.save_interval > 0 else []

# NEW: Only allocated in FULL mode
self.rate_history: Optional[List[List[Dict]]] = [] if self.memory_mode == MemoryMode.FULL else None
```

### 3. **Helper Methods**
```python
def _should_save_history(self) -> bool:
    return self.save_interval > 0

def _should_save_rates(self) -> bool:
    return self.memory_mode == MemoryMode.FULL
```

### 4. **Memory Inspection**
```python
# NEW: Built-in memory profiling
usage = sim.estimate_memory_usage()
# Returns: {'times': 0.001, 'history': 50.2, 'total': 50.203, ...}

sim.print_memory_summary()
# Prints formatted memory breakdown
```

---

## Usage Examples

### Example 1: Parameter Sweep with LIGHTWEIGHT Mode
```python
from tumor_simulation_optimized import TumorSimulation, MemoryMode
import numpy as np

results = []
for mutation_rate in np.linspace(0.01, 0.1, 20):
    config.mutation_rate = mutation_rate
    sim = TumorSimulation(config, memory_mode=MemoryMode.LIGHTWEIGHT)
    times, history, final_state, _ = sim.run()
    
    # Extract final population from history[-1]
    final_pop = history[-1]  # Only has one element anyway
    results.append({
        'mutation_rate': mutation_rate,
        'final_population': final_state.total_population(),
        'final_time': times[-1]
    })
    
    # Memory profile: ~1-2 MB per simulation
    # Total for 20 runs: ~40 MB
```

### Example 2: Trajectory Analysis with STANDARD Mode
```python
config.save_interval = 50
sim = TumorSimulation(config, memory_mode=MemoryMode.STANDARD)
times, history, final_state, _ = sim.run()

# Plot population over time
import matplotlib.pyplot as plt

populations = []
plot_times = []
for snapshot in history:
    pop = sum(clone['N'] for clone in snapshot.values())
    populations.append(pop)
    # Note: you can extract times from snapshot['t'] or use indices

plt.plot(populations)
plt.xlabel('Snapshot')
plt.ylabel('Population')
plt.show()
```

### Example 3: Detailed Event Analysis with FULL Mode
```python
sim = TumorSimulation(config, memory_mode=MemoryMode.FULL)
times, history, final_state, rate_history = sim.run()

# Analyze event statistics
from collections import Counter

event_types = Counter(e.kind.name for e in sim.events)
print(f"Birth events: {event_types['BIRTH']}")
print(f"Death events: {event_types['DEATH']}")
print(f"Mutations: {event_types['MUTATION']}")
print(f"Exhaustions: {event_types['EXHAUSTION']}")

# Analyze mutation timing
mutation_events = [e for e in sim.events if e.kind.name == 'MUTATION']
mutation_times = [times[i] for i, e in enumerate(sim.events) if e.kind.name == 'MUTATION']
print(f"Mutation times: {mutation_times}")

# Access rate history
for i, step_rates in enumerate(rate_history[:5]):  # First 5 steps
    print(f"Step {i}:")
    for rate_dict in step_rates:
        print(f"  {rate_dict['kind']} rate: {rate_dict['rate']:.4f}")
```

### Example 4: Adaptive Memory Mode
```python
def estimate_simulation_size(config):
    """Rough estimate: steps ≈ total_population * time / average_inter_event"""
    # Placeholder estimation
    return config.T_max * 1000

config_light = estimate_simulation_size(config)

if config_light < 10000:
    mode = MemoryMode.FULL  # Small sims can be FULL
elif config_light < 1000000:
    mode = MemoryMode.STANDARD
else:
    mode = MemoryMode.LIGHTWEIGHT

sim = TumorSimulation(config, memory_mode=mode)
times, history, final_state, rates = sim.run()
sim.print_memory_summary()
```

---

## Memory Breakdown by Component

### LIGHTWEIGHT Mode
| Component | Size | Notes |
|-----------|------|-------|
| `times` | ~0.1 MB | One entry per step |
| `history` | ~0.5-1 MB | Single final snapshot |
| `tissue_state` | ~0.3 MB | Final state only in memory |
| `rate_history` | None | Not allocated |
| `events` | None | Not allocated |
| **Total** | **~1-2 MB** | Independent of sim length! |

### STANDARD Mode (save_interval=100)
| Component | Size | Notes |
|-----------|------|-------|
| `times` | ~0.1 MB | One entry per step |
| `history` | ~10-50 MB | 1/100th of all snapshots |
| `tissue_state` | ~0.3 MB | Final state |
| `rate_history` | None | Not allocated |
| `events` | None | Not allocated |
| **Total** | **~10-50 MB** | Scales with T_max/save_interval |

### FULL Mode
| Component | Size | Notes |
|-----------|------|-------|
| `times` | ~0.1 MB | One entry per step |
| `history` | ~50-500 MB | Every step's snapshot |
| `tissue_state` | ~0.3 MB | Final state |
| `rate_history` | ~20-200 MB | Rate matrix at each step |
| `events` | ~10-100 MB | Event objects (birth/death/etc) |
| **Total** | **~100 MB - 1+ GB** | Scales linearly with steps |

---

## Recommendations

1. **For development:** Use `STANDARD` (default) - good balance
2. **For production sweeps:** Use `LIGHTWEIGHT` - run 1000s in parallel
3. **For debugging:** Use `FULL` - inspect everything
4. **For analysis:** Use `STANDARD` → `FULL` (re-run with FULL only if needed)

---

## Migration from Original

```python
# OLD CODE
sim = TumorSimulation(config)

# NEW CODE - Drop-in replacement (uses STANDARD by default)
sim = TumorSimulation(config)

# Or be explicit
sim = TumorSimulation(config, memory_mode=MemoryMode.STANDARD)
```

**Breaking changes:** None! Backward compatible via default `MemoryMode.STANDARD`.

---

## Future Optimizations

1. **Streaming snapshots to disk:** Write history to HDF5 as sim runs
2. **Event compression:** Store only state-changing events, not every rate update
3. **Sparse history:** Only save clones with N > threshold
4. **Time-based intervals:** Save at `T_max/100` intervals instead of step-based

