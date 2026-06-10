from io import StringIO

from src.gillespie.clone import Clone
from src.gillespie.cloneId import CloneId
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.infrastructure.display.summary_printer import print_summary


def _make_clone(clone_id: CloneId, N: int, config: SimulationConfig) -> Clone:
    c = Clone(clone_id=clone_id, N=N, config=config)
    c.birth_rate = 0.5
    c.death_rate = 0.2
    return c


def test_print_summary_outputs_final_time():
    config = SimulationConfig()
    clones = {(): _make_clone((), 10, config)}
    buf = StringIO()
    print_summary([1.0, 2.0], clones, top_k=5, file=buf)
    out = buf.getvalue()
    assert "final_time: 2.000000" in out


def test_print_summary_total_clones():
    config = SimulationConfig()
    clones = {
        (): _make_clone((), 10, config),
        (1,): _make_clone((1,), 5, config),
    }
    buf = StringIO()
    print_summary([1.0], clones, top_k=5, file=buf)
    assert "total_clones_created: 2" in buf.getvalue()


def test_print_summary_living_vs_dead():
    config = SimulationConfig()
    clones = {
        (): _make_clone((), 10, config),
        (1,): _make_clone((1,), 0, config),
    }
    buf = StringIO()
    print_summary([1.0], clones, top_k=5, file=buf)
    out = buf.getvalue()
    assert "living_clones: 1" in out
    assert "final_population: 10" in out


def test_print_summary_largest_clone():
    config = SimulationConfig()
    clones = {(): _make_clone((), 20, config)}
    clones[()].instability = 0.123
    buf = StringIO()
    print_summary([1.0], clones, top_k=5, file=buf)
    out = buf.getvalue()
    assert "largest_clone: root" in out
    assert "largest_clone_size: 20" in out


def test_print_summary_top_clones_limited_by_k():
    config = SimulationConfig()
    clones = {(): _make_clone((), 10, config), (1,): _make_clone((1,), 5, config)}
    buf = StringIO()
    print_summary([1.0], clones, top_k=1, file=buf)
    out = buf.getvalue()
    assert "lambda=0.5000" in out


def test_print_summary_empty_clones():
    buf = StringIO()
    print_summary([1.0], {}, top_k=5, file=buf)
    assert "final_population: 0" in buf.getvalue()
    assert "living_clones: 0" in buf.getvalue()
