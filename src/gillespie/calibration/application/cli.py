from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from src.gillespie.calibration.adapters.csv_output import CsvOutputAdapter
from src.gillespie.calibration.adapters.csv_reference import CsvReferenceAdapter
from src.gillespie.calibration.adapters.gillespie_adapter import (
    GillespieSimulatorAdapter,
)
from src.gillespie.calibration.application.abc_smc_service import AbcSmcService
from src.gillespie.calibration.domain.calibration_config import (
    CalibrationConfig,
    Prior,
    PriorDict,
)
from src.gillespie.simulation_config import SimulationConfig

logger = logging.getLogger(__name__)


def _parse_priors(priors_dict: dict) -> PriorDict:
    priors: PriorDict = {}
    for pname, spec in priors_dict.items():
        if isinstance(spec, dict):
            priors[pname] = Prior(
                param_name=pname,
                lo=spec.get("lo", 0.0),
                hi=spec.get("hi", 1.0),
                distribution=spec.get("distribution", "uniform"),
            )
        elif isinstance(spec, (list, tuple)) and len(spec) == 2:
            priors[pname] = Prior(
                param_name=pname,
                lo=float(spec[0]),
                hi=float(spec[1]),
                distribution="uniform",
            )
        else:
            raise ValueError(f"Invalid prior spec for '{pname}': {spec}")
    return priors


def _load_calibration_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def parse_calibration_args(args: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="calibrate",
        description="ABC-SMC calibration for Gillespie tumour simulator.",
    )
    p.add_argument(
        "--config",
        type=Path,
        required=True,
        help="JSON calibration config (priors, fixed params, CalibrationConfig).",
    )
    p.add_argument(
        "--reference",
        type=Path,
        required=True,
        help="Reference CSV (time, type, mean_N, std_N, n_replicates).",
    )
    p.add_argument(
        "--weights",
        type=Path,
        default=None,
        help="Precision weights JSON (optional).",
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("calibration_results"),
        help="Output directory for results. (default: calibration_results)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print config and exit without running simulations.",
    )
    p.add_argument(
        "--resume",
        action="store_true",
        help="Resume from last saved generation in output dir.",
    )
    p.add_argument(
        "--prior-sensitivity",
        type=int,
        default=0,
        help="Run prior sensitivity analysis with N samples (0 = skip).",
    )
    p.add_argument(
        "--particles",
        type=int,
        default=None,
        help="Override n_particles.",
    )
    p.add_argument(
        "--generations",
        type=int,
        default=None,
        help="Override n_generations.",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=None,
        help="Override n_workers.",
    )
    p.add_argument(
        "--reps",
        type=int,
        default=None,
        help="Override n_reps (stochastic replicates per particle).",
    )
    p.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Override random seed.",
    )
    p.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    return p.parse_args(args)


def calibrate_main(argv: Optional[List[str]] = None) -> None:
    args = parse_calibration_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
    )

    calib_data = _load_calibration_json(args.config)
    priors = _parse_priors(calib_data.get("priors", {}))
    fixed_params = calib_data.get("fixed_params", {})

    for pname in priors:
        if pname in fixed_params:
            logger.warning(
                "Parameter '%s' is in both priors and fixed_params. Removing from fixed_params.",
                pname,
            )
            del fixed_params[pname]

    calib_config = CalibrationConfig(
        **(calib_data.get("calibration_config", {}))
    )
    if args.particles is not None:
        calib_config = CalibrationConfig(
            **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
               "n_particles": args.particles}
        )
    if args.generations is not None:
        calib_config = CalibrationConfig(
            **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
               "n_generations": args.generations}
        )
    if args.workers is not None:
        calib_config = CalibrationConfig(
            **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
               "n_workers": args.workers}
        )
    if args.reps is not None:
        calib_config = CalibrationConfig(
            **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
               "n_reps": args.reps}
        )
    if args.seed is not None:
        calib_config = CalibrationConfig(
            **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
               "seed": args.seed}
        )
    calib_config = CalibrationConfig(
        **{**{f.name: getattr(calib_config, f.name) for f in calib_config.__class__.__dataclass_fields__.values()},
           "output_dir": str(args.output)}
    )

    base_config_dict = calib_data.get("base_config", {})
    base_config_dict.update(fixed_params)
    base_config = SimulationConfig(**base_config_dict)

    logger.info("Calibration config: %s", calib_config)
    logger.info("Priors: %s", {n: (p.lo, p.hi, p.distribution) for n, p in priors.items()})
    logger.info("Fixed params: %s", fixed_params)
    logger.info("Output dir: %s", args.output)

    if args.dry_run:
        logger.info("Dry run — exiting without simulation.")
        return

    simulator = GillespieSimulatorAdapter(base_config=base_config)
    reference = CsvReferenceAdapter(
        csv_path=str(args.reference),
        weights_path=str(args.weights) if args.weights else None,
    )
    output = CsvOutputAdapter(output_dir=str(args.output))

    rng = np.random.default_rng(calib_config.seed)
    service = AbcSmcService(
        config=calib_config,
        priors=priors,
        simulator=simulator,
        reference=reference,
        output=output,
        rng=rng,
    )

    if args.prior_sensitivity > 0:
        logger.info("Running prior sensitivity analysis (%d samples)...", args.prior_sensitivity)
        sens = service.prior_sensitivity(n_samples=args.prior_sensitivity)
        sens_path = args.output / "prior_sensitivity.csv"
        sens_path.parent.mkdir(parents=True, exist_ok=True)
        with open(sens_path, "w") as f:
            f.write("factor,n_samples,n_accepted,acceptance_rate,mean_distance,prior_rejects,sim_fails\n")
            for factor, info in sorted(sens.items()):
                f.write(
                    f"{factor},{info['n_samples']},{info['n_accepted']},"
                    f"{info['acceptance_rate']:.4f},{info['mean_distance']},"
                    f"{info['prior_rejects']},{info['sim_fails']}\n"
                )
        logger.info("Prior sensitivity results saved to %s", sens_path)

    if args.resume:
        logger.info("Resuming calibration...")
        result = service.resume()
    else:
        logger.info("Starting calibration...")
        result = service.run_abc_smc()

    logger.info("Calibration complete! %d generations.", result.n_generations())
    for s in result.all_summaries():
            logger.info(
                "  Gen %(generation)d | eps=%(epsilon)s | N=%(n_particles)d | "
                "ESS=%(ess).1f | acc=%(acceptance_rate).3f | dist=%(mean_distance).4f",
                s,
            )


if __name__ == "__main__":
    calibrate_main()
