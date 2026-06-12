from __future__ import annotations

import json
import platform
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


class CsvOutputAdapter:
    def __init__(self, output_dir: str):
        self.out = Path(output_dir)
        self.out.mkdir(parents=True, exist_ok=True)

    def save_generation(
        self,
        gen: int,
        particles: List[Dict[str, float]],
        weights: np.ndarray,
        distances: np.ndarray,
    ) -> None:
        if not particles:
            return
        df = pd.DataFrame(particles)
        df["weight"] = weights
        df["distance"] = distances
        path = self.out / f"gen_{gen:02d}.csv"
        df.to_csv(path, index=False)

    def save_epsilon_schedule(self, eps_history: List[float]) -> None:
        path = self.out / "epsilon_schedule.csv"
        df = pd.DataFrame({"generation": range(len(eps_history)), "epsilon": eps_history})
        df.to_csv(path, index=False)

    def save_generation_diagnostics(self, diagnostics: List[Dict]) -> None:
        path = self.out / "generation_diagnostics.csv"
        df = pd.DataFrame(diagnostics)
        df.to_csv(path, index=False)

    def save_manifest(
        self,
        calibration_config: dict,
        priors: dict,
    ) -> None:
        try:
            git_commit = (
                subprocess.run(
                    ["git", "rev-parse", "HEAD"],
                    capture_output=True,
                    text=True,
                    timeout=5,
                )
                .stdout.strip()
            )
        except Exception:
            git_commit = "unknown"

        manifest = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "git_commit": git_commit,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "calibration_config": {k: v for k, v in calibration_config.items()},
            "priors": {
                name: {
                    "param_name": p.param_name if hasattr(p, "param_name") else name,
                    "lo": p.lo,
                    "hi": p.hi,
                    "distribution": p.distribution,
                }
                for name, p in priors.items()
            },
        }

        path = self.out / "run_manifest.json"
        with open(path, "w") as f:
            json.dump(manifest, f, indent=2)

    def find_last_generation(self) -> Optional[int]:
        gens = []
        for f in self.out.glob("gen_*.csv"):
            try:
                gen = int(f.stem.split("_")[1])
                gens.append(gen)
            except (IndexError, ValueError):
                continue
        return max(gens) if gens else None

    def load_generation(
        self, gen: int
    ) -> Tuple[List[Dict[str, float]], np.ndarray, np.ndarray]:
        path = self.out / f"gen_{gen:02d}.csv"
        df = pd.read_csv(path)
        particles = df.drop(columns=["weight", "distance"]).to_dict(orient="records")
        weights = df["weight"].values
        distances = df["distance"].values
        return particles, weights, distances
