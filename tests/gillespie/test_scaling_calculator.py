import json

from src.gillespie.application.scaling_calculator import (
    PhysicalParameters,
    dimensional_to_dimensionless,
    dimensionless_to_dimensional,
    main,
    to_gillespie,
)


def test_report_parameter_conversion_and_inverse():
    physical = PhysicalParameters(
        lambda_C=2.0, lambda_I=1.0, mu_C=0.4, mu_E=0.6,
        theta=0.2, beta=0.3, omega=0.1, K_C=100.0, K_I=50.0,
        C0=20.0, I0=10.0, E0=5.0,
    )
    dimensionless = dimensional_to_dimensionless(physical)
    assert dimensionless.gamma_C == 0.2
    assert dimensionless.sigma == 5.0
    assert dimensionless.rho == 0.5
    restored = dimensionless_to_dimensional(
        dimensionless, lambda_C=2.0, K_C=100.0, K_I=50.0
    )
    assert restored == physical


def test_gillespie_interaction_scaling_and_absolute_counts():
    result = to_gillespie(
        PhysicalParameters(
            lambda_C=1, lambda_I=1, mu_C=0.1, mu_E=0.2,
            theta=0.04, beta=0.08, omega=0.02, K_C=1, K_I=0.5,
            C0=12.4, I0=8.6, E0=2.2,
        ),
        omega_system=4000,
    )
    assert result.theta_I == 0.04
    assert result.beta == 0.08
    assert result.cell_parameters["immune"]["omega_exhaust"] == 0.02
    assert result.cell_parameters["immune"]["K"] == 1 / 8000
    assert result.cell_parameters["mutated"]["N"] == 12
    assert result.cell_parameters["immune"]["N"] == 9


def test_cli_writes_json_with_dimensionless_and_gillespie_sections(tmp_path):
    output = tmp_path / "calculator.json"
    main([
        "--output", str(output), "--lambda-c", "2", "--lambda-i", "1",
        "--mu-c", "0.4", "--mu-e", "0.6", "--theta", "0.2",
        "--beta", "0.3", "--omega", "0.1", "--k-c", "100",
        "--k-i", "50", "--c0", "20", "--i0", "10", "--e0", "5",
    ])
    payload = json.loads(output.read_text())
    assert payload["_calculator"]["dimensionless"]["gamma_C"] == 0.2
    assert payload["OMEGA"] == 4000
    assert payload["cell_parameters"]["mutated"]["N"] == 20
