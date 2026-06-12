#!/usr/bin/env python3
"""
generate_synthetic_reference.py

Run N replicates of TumorSimulation with known ground-truth params.
Aggregate trajectories → mean_N, std_N per (time, type) → reference CSV + precision weights.

Usage:
    python scripts/generate_synthetic_reference.py --scenario homeostasis --n-replicates 50 --seed 42 --output-dir data/reference
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.gillespie.simulation_config import SimulationConfig
from src.gillespie.tumor_simulation import TumorSimulation


# ============================================================
# Ground truth params per scenario
# ============================================================

SCENARIOS = {
    "homeostasis": {
        "label": "Near-stable WT, no mutations, no immune response",
        "config": SimulationConfig(
            N0=50, N_mutant=0, N_immune=10, N_exhausted=0,
            lambda0=0.30, mu0=0.25, nu0=0.0,
            lambda_Immune=0.01, mu_Immune=0.01, mu_Exhausted=0.01,
            T_max=30, use_logistic=True, use_logistic_adapted=True,
            K0=100, K_immune=100, K_mutant=200,
            theta_I=0.0, beta=0.0, fitness_gain=0.0,
            verbose=False,
        ),
        "calibrate": ["lambda0", "mu0"],
    },
    "tumour-growth": {
        "label": "Net growth + mutation (Makefile-based)",
        "config": SimulationConfig(
            N0=20, N_mutant=0, N_immune=50, N_exhausted=0,
            lambda0=0.35, mu0=0.20, nu0=0.01,
            lambda_Immune=0.01, mu_Immune=0.005, mu_Exhausted=0.005,
            T_max=30, use_logistic=True, use_logistic_adapted=True,
            K0=500, K_immune=500, K_mutant=500,
            theta_I=0.0005, beta=0.0004, fitness_gain=0.02,
            instability_0=0.1,
            base_instability_buildup=0.005,
            mutation_instability_jump=0.05,
            verbose=False,
        ),
        "calibrate": ["lambda0", "mu0", "nu0"],
    },
    "immune-response": {
        "label": "Immune escape with exhaustion",
        "config": SimulationConfig(
            N0=30, N_mutant=5, N_immune=80, N_exhausted=0,
            lambda0=0.35, mu0=0.18, nu0=0.005,
            lambda_Immune=0.02, mu_Immune=0.006, mu_Exhausted=0.008,
            T_max=30, use_logistic=True, use_logistic_adapted=True,
            K0=500, K_immune=500, K_mutant=500,
            theta_I=0.001, beta=0.001, fitness_gain=0.15,
            verbose=False,
        ),
        "calibrate": ["theta_I", "beta", "fitness_gain"],
    },
}

CELL_TYPES = ["base", "mutated", "immune", "exhausted"]
SUMMARY_STAT_NAMES = ["mean", "final"]


# ============================================================
# Core logic
# ============================================================

def run_replicate(config: SimulationConfig, seed: int) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    c = SimulationConfig(**{**config.__dict__, "seed": seed, "verbose": False})
    sim = TumorSimulation(config=c)
    times, history, _ = sim.run()
    t = np.array(times)
    pop = {ct: np.zeros(len(t), dtype=float) for ct in CELL_TYPES}
    for i, snap in enumerate(history):
        for cid, cd in snap.items():
            typ = cd["Type"]
            if typ in pop:
                pop[typ][i] += cd["N"]
    return t, pop


_N_GRID_POINTS = 100


def _make_grid(t_reps: list[np.ndarray]) -> np.ndarray:
    t_min = min(t[0] for t in t_reps)
    t_max = max(t[-1] for t in t_reps)
    return np.linspace(t_min, t_max, _N_GRID_POINTS)


def aggregate(t_reps: list[np.ndarray], pop_reps: list[dict[str, np.ndarray]]) -> dict:
    grid = _make_grid(t_reps)

    interp = {ct: np.zeros((len(pop_reps), len(grid))) for ct in CELL_TYPES}
    for i, (t, pop) in enumerate(zip(t_reps, pop_reps)):
        for ct in CELL_TYPES:
            interp[ct][i] = np.interp(grid, t, pop[ct], left=pop[ct][0], right=pop[ct][-1])

    result = {"time": grid}
    for ct in CELL_TYPES:
        result[f"mean_{ct}"] = interp[ct].mean(axis=0)
        result[f"std_{ct}"] = interp[ct].std(axis=0, ddof=1)
    result["n_replicates"] = len(pop_reps)
    return result


def _active_types(agg: dict) -> list[str]:
    return [ct for ct in CELL_TYPES if agg[f"mean_{ct}"].max() > 5.0]


def compute_weights(agg: dict) -> dict[str, float]:
    w = {}
    active = _active_types(agg)
    for stat in SUMMARY_STAT_NAMES:
        for ct in active:
            key = f"{stat}_{ct}"
            mean_arr = agg[f"mean_{ct}"]
            std_arr = agg[f"std_{ct}"]
            if stat == "mean":
                var = np.mean(std_arr ** 2)
            elif stat == "final":
                var = std_arr[-1] ** 2
            else:
                var = 1.0
            w[key] = 1.0 / max(var, 1e-4)
    return w


def save_reference(agg: dict, weights: dict, out_dir: Path, scenario: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / f"reference_{scenario}.csv"
    n = len(agg["time"])
    with open(csv_path, "w") as f:
        cols = ["time", "type", "mean_N", "std_N", "n_replicates"]
        f.write(",".join(cols) + "\n")
        for i in range(n):
            t = agg["time"][i]
            for ct in CELL_TYPES:
                f.write(f"{t},{ct},{agg[f'mean_{ct}'][i]:.6f},{agg[f'std_{ct}'][i]:.6f},{agg['n_replicates']}\n")
    print(f"  saved {csv_path}")

    w_path = out_dir / f"reference_{scenario}_weights.json"
    with open(w_path, "w") as f:
        json.dump(weights, f, indent=2)
    print(f"  saved {w_path}")


def save_ground_truth(out_dir: Path, scenario: str, config: SimulationConfig) -> None:
    gt = {"scenario": scenario, "calibrate": SCENARIOS[scenario]["calibrate"]}
    for k in gt["calibrate"]:
        gt[k] = getattr(config, k)
    path = out_dir / f"reference_{scenario}_ground_truth.json"
    with open(path, "w") as f:
        json.dump(gt, f, indent=2)
    print(f"  saved {path}")


def plot_diagnostic(agg: dict, out_dir: Path, scenario: str, label: str) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        plt.set_loglevel("WARNING")
    except ImportError:
        print("  [skip] matplotlib not available — no diagnostic plot")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    t = agg["time"]
    colors = {"base": "#2196F3", "mutated": "#F44336", "immune": "#4CAF50", "exhausted": "#9E9E9E"}

    for ct in CELL_TYPES:
        mean = agg[f"mean_{ct}"]
        std = agg[f"std_{ct}"]
        if mean.max() < 0.5:
            continue
        ax.plot(t, mean, color=colors[ct], label=ct, linewidth=1.5)
        ax.fill_between(t, mean - std, mean + std, color=colors[ct], alpha=0.15)

    ax.set_xlabel("Time")
    ax.set_ylabel("Population")
    ax.set_title(f"{label} — {agg['n_replicates']} replicates")
    ax.legend()
    ax.grid(True, alpha=0.3)

    path = out_dir / f"reference_{scenario}.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {path}")


# ============================================================
# Main
# ============================================================

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenario", choices=list(SCENARIOS), required=True)
    parser.add_argument("--n-replicates", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("data/reference"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    scenario = args.scenario
    info = SCENARIOS[scenario]
    config = info["config"]
    rng = np.random.default_rng(args.seed)

    print(f"Scenario: {scenario} — {info['label']}")
    print(f"  replicates: {args.n_replicates}")
    print(f"  ground truth: {info['calibrate']}")
    print()

    t_reps = []
    pop_reps = []
    for rep in range(args.n_replicates):
        seed = int(rng.integers(1, 2_000_000_000))
        t, pop = run_replicate(config, seed)
        t_reps.append(t)
        pop_reps.append(pop)
        if (rep + 1) % 25 == 0:
            print(f"  replicate {rep + 1}/{args.n_replicates}")

    agg = aggregate(t_reps, pop_reps)
    weights = compute_weights(agg)
    out_dir = args.output_dir

    save_reference(agg, weights, out_dir, scenario)
    save_ground_truth(out_dir, scenario, config)

    if not args.no_plot:
        plot_diagnostic(agg, out_dir, scenario, info["label"])

    print(f"\nDone. Results in {out_dir.resolve()}")


if __name__ == "__main__":
    main()
