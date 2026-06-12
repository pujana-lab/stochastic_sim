import pytest

from src.gillespie.infrastructure.config.cli_config import build_parser, parse_args, get_explicit_cli_args


def test_build_parser_returns_parser():
    parser = build_parser()
    assert parser.prog is not None


def test_parse_args_defaults():
    ns = parse_args([])
    assert ns.N0 == 10
    assert ns.lambda0 == 0.20
    assert ns.T_max == 10.0
    assert ns.config is None


def test_parse_args_overrides():
    ns = parse_args(["--N0", "99", "--lambda0", "0.75", "--T-max", "42.0"])
    assert ns.N0 == 99
    assert ns.lambda0 == 0.75
    assert ns.T_max == 42.0


def test_parse_args_store_true():
    ns = parse_args(["--use-logistic"])
    assert ns.use_logistic is True


def test_parse_args_store_true_default_is_false():
    ns = parse_args([])
    assert ns.use_logistic is False


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
    ns = parse_args(["--N0", "50"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit == {"N0": 50}


def test_get_explicit_cli_args_dest_mapping():
    parser = build_parser()
    ns = parse_args(["--T-max", "100"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit.get("T_max") == 100.0


def test_get_explicit_cli_args_multiple():
    parser = build_parser()
    ns = parse_args(["--N0", "30", "--lambda0", "0.1", "--T-max", "5.0"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit == {"N0": 30, "lambda0": 0.1, "T_max": 5.0}


def test_get_explicit_cli_args_store_true_flag():
    parser = build_parser()
    ns = parse_args(["--use-logistic"])
    explicit = get_explicit_cli_args(ns, parser)
    assert explicit.get("use_logistic") is True


def test_get_explicit_cli_args_non_config_fields_excluded():
    parser = build_parser()
    ns = parse_args(["--top", "20"])
    explicit = get_explicit_cli_args(ns, parser)
    assert "top" not in explicit
    assert "save_history" not in explicit
