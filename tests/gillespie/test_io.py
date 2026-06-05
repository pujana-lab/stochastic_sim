import csv
import pytest
from pathlib import Path
from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.clone import Clone
from src.gillespie.cloneId import CloneId
from src.gillespie.infrastructure.csv_output import clone_id_to_str, save_history_csv, save_clones_csv


# ── clone_id_to_str ───────────────────────────────────────────────────────────

def test_root_clone_id():
    assert clone_id_to_str(()) == "root"

def test_single_level_clone_id():
    assert clone_id_to_str((1,)) == "1"

def test_nested_clone_id():
    assert clone_id_to_str((1, 2, 3)) == "1.2.3"


# ── save_history_csv ──────────────────────────────────────────────────────────

def make_history_data():
    times = [0.0, 1.0]
    history = [
        {(): {"N": 10, "rb": 0.5, "rd": 0.2}},
        {(): {"N": 12, "rb": 0.48, "rd": 0.2}, (1,): {"N": 1, "rb": 0.55, "rd": 0.2}},
    ]
    return times, history


def test_save_history_csv_creates_file(tmp_path):
    path = tmp_path / "history.csv"
    times, history = make_history_data()
    save_history_csv(path, times, history)
    assert path.exists()


def test_save_history_csv_header(tmp_path):
    path = tmp_path / "history.csv"
    times, history = make_history_data()
    save_history_csv(path, times, history)
    with path.open() as f:
        reader = csv.reader(f)
        header = next(reader)
    assert header == ["time", "clone_id", "N", "rb", "rd"]


def test_save_history_csv_row_count(tmp_path):
    path = tmp_path / "history.csv"
    times, history = make_history_data()
    save_history_csv(path, times, history)
    with path.open() as f:
        rows = list(csv.reader(f))
    # 1 header + 1 clone at t=0 + 2 clones at t=1 = 4
    assert len(rows) == 4


def test_save_history_csv_clone_id_encoding(tmp_path):
    path = tmp_path / "history.csv"
    save_history_csv(path, [0.0], [{(): {"N": 5, "rb": 0.5, "rd": 0.2}}])
    with path.open() as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["clone_id"] == "root"


# ── save_clones_csv ───────────────────────────────────────────────────────────

def make_clones() -> dict:
    config = SimulationConfig()
    root = Clone(clone_id=(), N=10, config=config)
    child = Clone(clone_id=(1,), N=3, config=config, parent=())
    return {(): root, (1,): child}


def test_save_clones_csv_creates_file(tmp_path):
    path = tmp_path / "clones.csv"
    save_clones_csv(path, make_clones())
    assert path.exists()


def test_save_clones_csv_header(tmp_path):
    path = tmp_path / "clones.csv"
    save_clones_csv(path, make_clones())
    with path.open() as f:
        header = next(csv.reader(f))
    assert "clone_id" in header
    assert "N" in header
    assert "birth_rate" in header


def test_save_clones_csv_row_count(tmp_path):
    path = tmp_path / "clones.csv"
    clones = make_clones()
    save_clones_csv(path, clones)
    with path.open() as f:
        rows = list(csv.reader(f))
    assert len(rows) == len(clones) + 1  # +1 header


def test_save_clones_csv_parent_encoding(tmp_path):
    path = tmp_path / "clones.csv"
    save_clones_csv(path, make_clones())
    with path.open() as f:
        rows = list(csv.DictReader(f))
    root_row = next(r for r in rows if r["clone_id"] == "root")
    child_row = next(r for r in rows if r["clone_id"] == "1")
    assert root_row["parent"] == ""
    assert child_row["parent"] == "root"
