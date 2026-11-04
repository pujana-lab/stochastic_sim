from scripts import model_diagram


def test_model_diagram_is_up_to_date():
    """docs/model_diagram.mmd must match the code. Fix with `make diagram`."""
    assert model_diagram.main(["--check"]) == 0


def test_every_registered_type_is_in_the_diagram():
    from src.gillespie.clone import Clone
    mmd = model_diagram.build_mermaid()
    for t in Clone._registry:
        assert f"    {t}(" in mmd


def test_node_formulas_follow_the_code():
    """Each node writes its rates as in clone.py: α_F multiplies r_B, β and θ_I add."""
    mmd = model_diagram.build_mermaid()
    assert "r<sub>B</sub> = λ<sub>mutated</sub> · N<sub>mutated</sub> · g · (1 + α<sub>F</sub> · N<sub>caf</sub>)" in mmd
    assert "r<sub>B</sub> = λ<sub>immune</sub> · N<sub>immune</sub> · g + β · N<sub>immune</sub> · N<sub>mutated</sub>" in mmd
    assert "r<sub>D</sub> = μ<sub>mutated</sub> · N<sub>mutated</sub> + θ<sub>I</sub> · N<sub>mutated</sub> · N<sub>immune</sub>" in mmd


def test_every_edge_symbol_appears_in_its_target_formula():
    """The arrow label must be findable in the formula of the node it acts on."""
    import re
    mmd = model_diagram.build_mermaid()
    nodes = dict(re.findall(r'^    (\w+)\("(.*)"\)$', mmd, re.M))
    for src, label, dst in re.findall(r'^    (\w+) (?:==>|--x|-\.->)\|"(.*?)"\| (\w+)$', mmd, re.M):
        owner = dst.removeprefix("conv_")
        for symbol in label.split(", "):
            assert symbol in nodes[owner], (src, symbol, owner)
