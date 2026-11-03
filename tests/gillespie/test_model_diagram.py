from scripts import model_diagram


def test_model_diagram_is_up_to_date():
    """docs/model_diagram.mmd must match the code. Fix with `make diagram`."""
    assert model_diagram.main(["--check"]) == 0


def test_every_registered_type_is_in_the_diagram():
    from src.gillespie.clone import Clone
    mmd = model_diagram.build_mermaid()
    for t in Clone._registry:
        assert f"    {t}(" in mmd
