"""Text table of the run conditions: input values next to the values actually used.

Pure formatting (returns strings, no I/O). Input = defaults + JSON + CLI before OMEGA
scaling (SimulationConfig.input_values); used = after scaling (SimulationConfig.as_dict()).
"""
from src.gillespie.simulation_config import SimulationConfig

CHANGED = "*"   # marks a value changed by the OMEGA scaling (or by fitness_gain)


def _fmt(value) -> str:
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, dict):
        return "{" + ", ".join(f"{k}: {_fmt(v)}" for k, v in value.items()) + "}"
    return str(value)


def config_rows(config: SimulationConfig) -> list[tuple[str, str, str, bool]]:
    """(parameter, input, used, changed) for every parameter; cell types as "type.param"."""
    given, used = config.input_values, config.as_dict()
    rows = [(name, _fmt(given[name]), _fmt(used[name]), given[name] != used[name])
            for name in given if name != "cell_parameters"]
    for cell_type, params in used["cell_parameters"].items():
        for name, value in params.items():
            before = given["cell_parameters"][cell_type][name]
            rows.append((f"{cell_type}.{name}", _fmt(before), _fmt(value), before != value))
    return rows


def format_config(config: SimulationConfig) -> str:
    rows = config_rows(config)
    width = [max(len(r[i]) for r in rows + [("parameter", "input", "used", False)]) for i in range(3)]
    line = lambda p, i, u, mark="": f"{p:<{width[0]}}  {i:<{width[1]}}  {u:<{width[2]}}  {mark}".rstrip()
    out = [line("parameter", "input", "used"), "-" * (sum(width) + 6)]
    out += [line(p, i, u, CHANGED if changed else "") for p, i, u, changed in rows]
    out.append(f"({CHANGED} = changed by the OMEGA scaling, OMEGA = {config.OMEGA}; scale = {config.scale})")
    return "\n".join(out)
