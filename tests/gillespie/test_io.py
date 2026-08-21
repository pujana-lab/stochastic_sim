import pytest
import pandas as pd
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone_factory import CloneFactory
from src.gillespie.infrastructure.csv_output import clone_id_to_str, save_clones_parquet, save_history_parquet


# ── clone_id_to_str ───────────────────────────────────────────────────────────

def test_root_clone_id():
    assert clone_id_to_str(()) == "root"

def test_single_level_clone_id():
    assert clone_id_to_str((1,)) == "1"

def test_nested_clone_id():
    assert clone_id_to_str((1, 2, 3)) == "1.2.3"


# ── save_history_parquet ─────────────────────────────────────────────────────

def make_history_data():
    times = [0.0, 1.0]
    history = [
        {(): {"Type": "base", "N": 10, "rb": 0.5, "rd": 0.2}},
        {(): {"Type": "base", "N": 12, "rb": 0.48, "rd": 0.2}, (1,): {"Type": "mutated", "N": 1, "rb": 0.55, "rd": 0.2}},
    ]
    return times, history


def test_save_history_parquet_creates_file(tmp_path):
    path = tmp_path / "history.parquet"
    times, history = make_history_data()
    save_history_parquet(path, times, history)
    assert path.exists()


def test_save_history_parquet_header(tmp_path):
    path = tmp_path / "history.parquet"
    times, history = make_history_data()
    save_history_parquet(path, times, history)
    df = pd.read_parquet(path)
    assert list(df.columns) == ["time", "type", "clone_id", "N", "rb", "rd"]


def test_save_history_parquet_row_count(tmp_path):
    path = tmp_path / "history.parquet"
    times, history = make_history_data()
    save_history_parquet(path, times, history)
    df = pd.read_parquet(path)
    assert len(df) == 3


def test_save_history_parquet_clone_id_encoding(tmp_path):
    path = tmp_path / "history.parquet"
    save_history_parquet(path, [0.0], [{(): {"Type": "base", "N": 5, "rb": 0.5, "rd": 0.2}}])
    df = pd.read_parquet(path)
    assert df.iloc[0]["clone_id"] == "root"


# ── save_clones_parquet ──────────────────────────────────────────────────────

def make_clones() -> dict:
    config = SimulationConfig()
    factory = CloneFactory(config)
    root = factory.create_clone(clone_id=(), clone_type="base", N=10)
    child = factory.create_clone(clone_id=(1,), clone_type="base", N=3, parent=())
    return {(): root, (1,): child}


def test_save_clones_parquet_creates_file(tmp_path):
    path = tmp_path / "clones.parquet"
    save_clones_parquet(path, make_clones())
    assert path.exists()


def test_save_clones_parquet_header(tmp_path):
    path = tmp_path / "clones.parquet"
    save_clones_parquet(path, make_clones())
    header = list(pd.read_parquet(path).columns)
    assert "clone_id" in header
    assert "N" in header
    assert "birth_rate" in header


def test_save_clones_parquet_row_count(tmp_path):
    path = tmp_path / "clones.parquet"
    clones = make_clones()
    save_clones_parquet(path, clones)
    df = pd.read_parquet(path)
    assert len(df) == len(clones)


def test_save_clones_parquet_parent_encoding(tmp_path):
    path = tmp_path / "clones.parquet"
    save_clones_parquet(path, make_clones())
    rows = pd.read_parquet(path).to_dict(orient="records")
    root_row = next(r for r in rows if r["clone_id"] == "root")
    child_row = next(r for r in rows if r["clone_id"] == "1")
    assert root_row["parent"] == ""
    assert child_row["parent"] == "root"
