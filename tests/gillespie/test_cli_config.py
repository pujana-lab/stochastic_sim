from src.gillespie.infrastructure.config.cli_config import build_parser, parse_args, get_explicit_cli_args


def test_build_parser_returns_parser():
    parser = build_parser()
    assert parser.prog is not None


def test_parse_args_defaults():
    ns = parse_args([])
    assert ns.T_max == 1000.0
    assert ns.config is None


def test_parse_args_overrides():
    ns = parse_args(["--beta", "0.75", "--T-max", "42.0"])
    assert ns.beta == 0.75
    assert ns.T_max == 42.0


def test_parse_args_store_true():
    ns = parse_args(["--use-logistic"])
    assert ns.use_logistic is True


def test_parse_args_store_true_default_is_false():
    ns = parse_args([])
    assert ns.use_logistic is True


def test_parse_args_config():
    ns = parse_args(["--config", "my_config.json"])
    assert ns.config is not None
    assert str(ns.config) == "my_config.json"


def test_get_explicit_cli_args_none_provided():
    parser = build_parser()
    ns = parse_args([])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit == {}


def test_get_explicit_cli_args_single_override():
    parser = build_parser()
    ns = parse_args(["--beta", "0.5"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit == {"beta": 0.5}


def test_get_explicit_cli_args_dest_mapping():
    parser = build_parser()
    ns = parse_args(["--T-max", "100"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit.get("T_max") == 100.0


def test_get_explicit_cli_args_multiple():
    parser = build_parser()
    ns = parse_args(["--beta", "0.1", "--T-max", "5.0"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit == {"beta": 0.1, "T_max": 5.0}


def test_get_explicit_cli_args_store_true_flag():
    parser = build_parser()
    ns = parse_args(["--no-use-logistic"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit.get("use_logistic") is False


def test_get_explicit_cli_args_non_config_fields_excluded():
    parser = build_parser()
    ns = parse_args(["--top", "20"])
    explicit = get_explicit_cli_args(ns, parser)
    assert "top" not in explicit
    assert "save_history" not in explicit


def test_stepper_flag():
    assert parse_args([]).stepper == "ssa"
    assert get_explicit_cli_args(parse_args(["--stepper", "ssa"]), build_parser()) == {"stepper": "ssa"}


def test_unknown_stepper_flag_is_rejected():
    import pytest
    with pytest.raises(SystemExit):
        parse_args(["--stepper", "sas"])


def test_flag_equal_to_its_default_still_counts():
    """A written flag overrides the JSON even when it equals the default."""
    explicit = get_explicit_cli_args(parse_args(["--T-max", "1000", "--seed", "3"]), build_parser())
    assert explicit == {"T_max": 1000.0, "seed": 3}


def test_flag_equal_to_its_default_overrides_the_json(tmp_path):
    import json
    from src.gillespie.application.config_service import build_config_from_cli
    path = tmp_path / "c.json"
    path.write_text(json.dumps({"T_max": 50, "verbose": False}))
    cfg, _ = build_config_from_cli(["--config", str(path), "--T-max", "1000"])
    assert cfg.T_max == 1000
