"""Generate the model diagram (Mermaid) from the real code.

The diagram is *inferred*, not drawn by hand:

* cell types        -> ``Clone._registry`` (every subclass declared with ``clone_type=``)
* default params    -> ``_DEFAULT_CELL_PARAMETERS`` in ``simulation_config.py``
* interactions      -> AST of each subclass in ``clone.py``: every
  ``pop_map.get("<type>")`` / ``pop_map["<type>"]`` read inside a rate method
  becomes an edge, and every ``self.config.<param>`` read becomes its label
* mutation target   -> ``CellTypeConfig.next_mutation``
* conversion target -> ``CellTypeConfig.converts_to`` when that field exists,
  otherwise the ``cell_parameters["<type>"]`` read in ``TumorSimulation._induce_exhaustion``

Edge semantics (by the method in which the read happens):

    crowding_pool (config)  "Y_T = N_T + N_S + ..." in T's label (crowding g(Y_T, K_T)) (or S -.-o T with --crowding-edges)
    birth_rate_effective    S ==> T    S raises T's birth (activation / support)
    death_rate_effective    S --x T    S kills T
    exhaustion_rate_effective /
    conversion_rate_effective
                            S -.-> (T -> target)   S catalyses T's conversion; the label is the
                                                   parameter that multiplies N_S

Usage::

    python -m scripts.model_diagram                 # print Mermaid to stdout
    python -m scripts.model_diagram --write         # update docs/model_diagram.mmd
                                                    # and every <!-- model-diagram --> block
    python -m scripts.model_diagram --check         # exit 1 if docs are stale (CI / pytest)
    python -m scripts.model_diagram --write --svg   # also render docs/model_diagram.svg (needs mmdc)
"""
from __future__ import annotations

import argparse
import ast
import inspect
import re
import shutil
import subprocess
import sys
from dataclasses import fields
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.gillespie import clone as clone_module  # noqa: E402
from src.gillespie import tumor_simulation as sim_module  # noqa: E402
from src.gillespie.clone import Clone  # noqa: E402
from src.gillespie.simulation_config import (  # noqa: E402
    _DEFAULT_CELL_PARAMETERS,
    CellTypeConfig,
)

MMD_PATH = ROOT / "docs" / "model_diagram.mmd"
SVG_PATH = ROOT / "docs" / "model_diagram.svg"
MARKDOWN_TARGETS = [ROOT / "docs" / "ARCHITECTURE.md", ROOT / "docs" / "PLAN_fibroblastos.md"]
BLOCK_RE = re.compile(r"(<!-- model-diagram:start -->\n)(.*?)(<!-- model-diagram:end -->)", re.S)

# Mean-field symbols of the report (C, I, E) and the
# fibroblast plan (F_N, F_T). Types without one (base) show only their name.
SYMBOLS = {
    "mutated": "C",
    "immune": "I",
    "exhausted": "E",
    "fibroblast": "F<sub>N</sub>",
    "caf": "F<sub>T</sub>",
}

# Code parameter -> symbol used in README.md (θ_I, β, ω, ν ...).
PARAM_SYMBOLS = {
    "theta_I": "θ<sub>I</sub>",
    "beta": "β",
    "omega_exhaust": "ω",
    "nu": "ν",
    "kappa": "κ",
    "C_half": "C<sub>½</sub>",
    "omega_F": "ω<sub>F</sub>",
    "alpha_F": "α<sub>F</sub>",
}


def _psym(params: list[str]) -> str:
    """Symbols of ``params``, in the order of PARAM_SYMBOLS (unknown ones last)."""
    order = list(PARAM_SYMBOLS)
    ranked = sorted(params, key=lambda p: order.index(p) if p in order else len(order))
    return ", ".join(PARAM_SYMBOLS.get(p, p) for p in ranked)

RATE_METHODS = {
    "crowding_numerator": "crowding",
    "birth_rate_effective": "birth",
    "death_rate_effective": "death",
    "exhaustion_rate_effective": "conversion",
    "conversion_rate_effective": "conversion",
}


# ── code inspection ──────────────────────────────────────────────────────────

def _pop_keys_in_order(node: ast.AST) -> list[str]:
    """Type names read from any ``pop_map`` inside ``node``, in source order."""
    found: list[tuple[int, int, str]] = []
    for sub in ast.walk(node):
        # pop_map.get("x", 0)
        if (
            isinstance(sub, ast.Call)
            and isinstance(sub.func, ast.Attribute)
            and sub.func.attr == "get"
            and _ends_with_pop_map(sub.func.value)
            and sub.args
            and isinstance(sub.args[0], ast.Constant)
        ):
            found.append((sub.lineno, sub.col_offset, sub.args[0].value))
        # pop_map["x"]
        if isinstance(sub, ast.Subscript) and _ends_with_pop_map(sub.value):
            sl = sub.slice
            if isinstance(sl, ast.Constant):
                found.append((sub.lineno, sub.col_offset, sl.value))
    return list(dict.fromkeys(key for *_, key in sorted(found)))


def _pop_keys(node: ast.AST) -> set[str]:
    """Type names read from any ``pop_map`` inside ``node``."""
    return set(_pop_keys_in_order(node))


def _ends_with_pop_map(node: ast.AST) -> bool:
    if isinstance(node, ast.Name):
        return node.id == "pop_map"
    if isinstance(node, ast.Attribute):
        return node.attr == "pop_map"
    return False


def _config_params(node: ast.AST) -> list[str]:
    """``self.config.<param>`` reads inside ``node`` (unique)."""
    out: list[str] = []
    for sub in ast.walk(node):
        if (
            isinstance(sub, ast.Attribute)
            and isinstance(sub.value, ast.Attribute)
            and sub.value.attr == "config"
            and sub.attr not in ("crowding_strategy",)
            and sub.attr not in out
        ):
            out.append(sub.attr)
    return out


def _params_by_source(fn: ast.FunctionDef) -> dict[str, list[str]]:
    """Attribute each ``self.config.<param>`` to the population(s) it multiplies.

    ``n_caf = pop_map.get("caf", 0)`` binds ``n_caf -> "caf"``; any later
    statement that uses ``n_caf`` contributes its config params to ``"caf"``.
    """
    var2key: dict[str, str] = {}
    for st in ast.walk(fn):
        if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name):
            keys = _pop_keys(st.value)
            if len(keys) == 1:
                var2key[st.targets[0].id] = next(iter(keys))
    out: dict[str, list[str]] = {}
    for st in ast.walk(fn):
        if not isinstance(st, (ast.Assign, ast.AugAssign, ast.Return)):
            continue
        used = {n.id for n in ast.walk(st) if isinstance(n, ast.Name)}
        keys = {var2key[u] for u in used if u in var2key} | _pop_keys(st)
        for k in keys:
            bucket = out.setdefault(k, [])
            bucket += [p for p in _config_params(st) if p not in bucket]
    return out


# Rate of the base ``Clone`` method that a ``super().<method>()`` call returns
# (README "Base forms"), written with the type name as in docs/NOTATION.md.
BASE_RATES = {
    "birth_rate_effective": "λ<sub>{t}</sub> · N<sub>{t}</sub> · g",
    "death_rate_effective": "μ<sub>{t}</sub> · N<sub>{t}</sub>",
    "exhaustion_rate_effective": "ω<sub>{t}</sub> · N<sub>{t}</sub>",
    "conversion_rate_effective": "ω<sub>{t}</sub> · N<sub>{t}</sub>",
}


class _RateExpr:
    """Render the return expression of a rate method in README notation.

    Local variables are inlined, ``self.N`` -> N_t, ``pop_map.get("x")`` -> N_x,
    ``self.config.p`` -> its symbol and ``super().m()`` -> the base rate.
    """

    def __init__(self, fn: ast.FunctionDef, t: str) -> None:
        self.t = t
        self.method = fn.name
        self.env = {
            st.targets[0].id: st.value
            for st in ast.walk(fn)
            if isinstance(st, ast.Assign) and len(st.targets) == 1 and isinstance(st.targets[0], ast.Name)
        }
        returns = [st for st in ast.walk(fn) if isinstance(st, ast.Return) and st.value is not None]
        # the last return in source order; earlier ones are guards (``return 0.0``)
        self.ret = max(returns, key=lambda r: r.lineno).value if returns else None

    def _resolve(self, n: ast.AST) -> ast.AST:
        seen = set()
        while isinstance(n, ast.Name) and n.id in self.env and n.id not in seen:
            seen.add(n.id)
            n = self.env[n.id]
        return n

    def _split(self, n: ast.AST, op: type) -> list[ast.AST]:
        n = self._resolve(n)
        if isinstance(n, ast.BinOp) and isinstance(n.op, op):
            return self._split(n.left, op) + self._split(n.right, op)
        return [n]

    @staticmethod
    def _is_super(n: ast.AST) -> bool:
        return (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and isinstance(n.func.value, ast.Call)
            and isinstance(n.func.value.func, ast.Name)
            and n.func.value.func.id == "super"
        )

    def _kind(self, n: ast.AST) -> int:
        """Sort key inside a product: base rate, parameters, populations, rest."""
        n = self._resolve(n)
        if self._is_super(n):
            return 0
        if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Attribute) and n.value.attr == "config":
            return 1
        if self._keys(n) and not isinstance(n, ast.BinOp):
            return 2
        return 3

    def _keys(self, n: ast.AST) -> set[str]:
        """Populations an expression depends on."""
        n = self._resolve(n)
        if self._is_super(n):
            return {self.t}
        if isinstance(n, ast.Attribute) and n.attr == "N" and isinstance(n.value, ast.Name):
            return {self.t}
        keys = set(_pop_keys(n))
        for sub in ast.iter_child_nodes(n):
            if sub is not n:
                keys |= self._keys(sub)
        return keys

    def render(self, n: ast.AST) -> str:
        n = self._resolve(n)
        if self._is_super(n):
            return BASE_RATES.get(n.func.attr, f"{n.func.attr}()").format(t=self.t)
        if isinstance(n, ast.Attribute) and n.attr == "N" and isinstance(n.value, ast.Name):
            return f"N<sub>{self.t}</sub>"
        if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Attribute) and n.value.attr == "config":
            return PARAM_SYMBOLS.get(n.attr, n.attr)
        keys = _pop_keys_in_order(n)
        if len(keys) == 1 and isinstance(n, (ast.Call, ast.Subscript)):
            return f"N<sub>{keys[0]}</sub>"
        if isinstance(n, ast.Constant):
            return f"{n.value:g}" if isinstance(n.value, float) else str(n.value)
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Mult):
            return self.product(self._split(n, ast.Mult))
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Add):
            return " + ".join(self.render(a) for a in self._split(n, ast.Add))
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Sub):
            return f"{self.render(n.left)} − {self._paren(n.right)}"
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.Div):
            return f"{self.render(n.left)} / {self._paren(n.right)}"
        return ast.unparse(n)

    def _paren(self, n: ast.AST) -> str:
        n = self._resolve(n)
        text = self.render(n)
        return f"({text})" if isinstance(n, ast.BinOp) and isinstance(n.op, (ast.Add, ast.Sub)) else text

    def product(self, factors: list[ast.AST]) -> str:
        return " · ".join(self._paren(f) for f in sorted(factors, key=self._kind))

    def formula(self) -> str | None:
        """The whole rate, e.g. ``λ_immune · N_immune · g + β · N_immune · N_mutated``."""
        if self.ret is None:
            return None
        return " + ".join(self.product(self._split(a, ast.Mult)) for a in self._split(self.ret, ast.Add))


def _returns_zero(fn: ast.FunctionDef) -> bool:
    body = [s for s in fn.body if not (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant))]
    return (
        len(body) == 1
        and isinstance(body[0], ast.Return)
        and isinstance(body[0].value, ast.Constant)
        and body[0].value.value in (0, 0.0)
    )


def inspect_clone_types() -> dict[str, dict]:
    """Per registered type: overridden rate methods, pop_map reads and config params."""
    tree = ast.parse(inspect.getsource(clone_module))
    class_nodes = {n.name: n for n in tree.body if isinstance(n, ast.ClassDef)}

    info: dict[str, dict] = {}
    for type_name, cls in Clone._registry.items():
        methods: dict[str, dict] = {}
        # walk MRO so inherited overrides (from intermediate subclasses) count too
        for klass in reversed(cls.__mro__):
            node = class_nodes.get(klass.__name__)
            if node is None or klass is Clone:
                continue
            for fn in node.body:
                if isinstance(fn, ast.FunctionDef) and fn.name in RATE_METHODS:
                    methods[fn.name] = {
                        "reads": _pop_keys(fn),
                        "reads_in_order": _pop_keys_in_order(fn),
                        "params": _config_params(fn),
                        "by_src": _params_by_source(fn),
                        "zero": _returns_zero(fn),
                        "formula": _RateExpr(fn, type_name).formula(),
                    }
        # Crowding is data (CellTypeConfig.crowding_pool), not a per-class override:
        # read it from the defaults, keeping the declared order for the "Y_X = ..." label.
        cfg = _DEFAULT_CELL_PARAMETERS.get(type_name)
        pool = tuple(getattr(cfg, "crowding_pool", {}) or {})
        if pool and "crowding_numerator" not in methods:
            methods["crowding_numerator"] = {
                "reads": set(pool), "reads_in_order": list(pool),
                "params": [], "by_src": {}, "zero": False, "formula": None,
            }
        info[type_name] = {"cls": cls, "methods": methods}
    return info


def _conversion_target(type_name: str, cfg: CellTypeConfig) -> str | None:
    if "converts_to" in {f.name for f in fields(CellTypeConfig)}:
        return getattr(cfg, "converts_to", "") or None
    # legacy: TumorSimulation._induce_exhaustion names the target type
    src = inspect.getsource(sim_module.TumorSimulation._induce_exhaustion)
    for tgt in re.findall(r"cell_parameters\[[\"'](\w+)[\"']\]", src):
        if tgt in _DEFAULT_CELL_PARAMETERS:
            return tgt
    return None


# ── Mermaid rendering ────────────────────────────────────────────────────────

def _sym(t: str) -> str:
    return SYMBOLS.get(t, t)


def _wsym(t: str) -> str:
    """Subscript of the crowding weight w_T (docs/NOTATION.md): w_C, w_I, w_T, w_N..."""
    m = re.fullmatch(r"F<sub>(\w+)</sub>", _sym(t))
    return m.group(1) if m else _sym(t)


def _legend(colors: dict[str, str], crowding_edges: bool) -> str:
    """Label of the legend node: edge kinds, rate symbols and crowding."""
    def c(kind: str, glyph: str) -> str:
        return f"<b><span style='color:{colors[kind]}'>{glyph}</span></b>"

    rows = [
        "<b>Legend</b>",
        "<b>⟶</b> transition: mutation r<sub>M</sub>, conversion r<sub>E</sub>",
        "<span style='color:#888'>●</span> transition point",
        f"{c('birth', '⟹')} source enters the target's birth r<sub>B</sub>",
        f"{c('death', '⟶✕')} source enters the target's death r<sub>D</sub>",
        f"{c('cat', '⇢')} source enters the conversion r<sub>E</sub>",
    ]
    if crowding_edges:
        rows.append(f"{c('crowding', '⇢○')} source competes for space")
    rows += [
        "arrow label: parameter of the interaction;",
        "the term is in the target's formula (in its node)",
        "r<sub>E</sub> of a type is written in the converting node",
        "g = g(Y<sub>X</sub>, K<sub>X</sub>) = 1 − Y<sub>X</sub>/K<sub>X</sub>, Y<sub>X</sub>: competing population",
    ]
    return "<br/>".join(f"<small>{r}</small>" if i else r for i, r in enumerate(rows))


def build_mermaid(crowding_edges: bool = False) -> str:
    info = inspect_clone_types()
    defaults = _DEFAULT_CELL_PARAMETERS
    lines = [
        "%% AUTO-GENERATED by scripts/model_diagram.py from src/gillespie/clone.py",
        "%% and simulation_config.py. Do not edit by hand: run `make diagram`.",
        "flowchart TB",
    ]
    edges: list[str] = []
    styles: list[str] = []
    link_styles: dict[str, list[int]] = {
        "crowding": [], "birth": [], "death": [], "cat": [], "flow": [], "anchor": [],
    }
    has_incoming: set[str] = set()
    catalysts: dict[str, list[str]] = {}

    def add(edge: str, kind: str, target: str | None = None) -> None:
        link_styles[kind].append(len(edges))
        edges.append(edge)
        if target is not None:
            has_incoming.add(target)

    for t, meta in info.items():
        cfg = defaults.get(t, CellTypeConfig())
        m = meta["methods"]
        birth_off = m.get("birth_rate_effective", {}).get("zero", False) or cfg.lambda0 == 0
        conv = m.get("exhaustion_rate_effective") or m.get("conversion_rate_effective")
        target = _conversion_target(t, cfg) if conv or cfg.omega_exhaust > 0 else None

        def rate(meth: str, symbol: str) -> str:
            formula = m[meth]["formula"] if meth in m else None
            return f"r<sub>{symbol}</sub> = {formula or BASE_RATES[meth].format(t=t)}"

        lines_lbl = []
        if not birth_off:
            lines_lbl.append(rate("birth_rate_effective", "B"))
        if cfg.mu > 0 or "death_rate_effective" in m:
            lines_lbl.append(rate("death_rate_effective", "D"))
        if target:
            meth = "conversion_rate_effective" if "conversion_rate_effective" in m else "exhaustion_rate_effective"
            lines_lbl.append(rate(meth, "E"))
        crowd = m.get("crowding_numerator")
        if crowd and not birth_off:
            weighted = {k for k, w in (getattr(cfg, "crowding_pool", {}) or {}).items() if w != 1}   # w_T = 1 not drawn
            num = " + ".join(
                (f"w<sub>{_wsym(k)}</sub> · " if k in weighted else "") + f"N<sub>{k}</sub>"
                for k in crowd["reads_in_order"]
            ) or "0"
            lines_lbl.append(f"Y<sub>{t}</sub> = {num}")
        head = f"<b>{t}</b>" + (f" ({_sym(t)})" if t in SYMBOLS else "")
        label = head + "".join(f"<br/><small>{x}</small>" for x in lines_lbl)
        lines.append(f'    {t}("{label}")')

        # mutation
        if cfg.next_mutation:
            node = f"mut_{t}"
            lines.append(f"    {node}(( ))")
            add(f"    {t} --- {node}", "flow")
            add(f'    {node} -->|"r<sub>M</sub> (ν)"| {cfg.next_mutation}', "flow", cfg.next_mutation)
            styles.append(f"    class {node} tr")

        # conversion (exhaustion / activation)
        if target:
            node = f"conv_{t}"
            lines.append(f"    {node}(( ))")
            add(f"    {t} --- {node}", "flow")
            add(f'    {node} -->|"r<sub>E</sub>"| {target}', "flow", target)
            styles.append(f"    class {node} tr")
            for src in sorted((conv or {}).get("reads", ())):
                params = _psym(conv["by_src"].get(src, []))
                # The base method multiplies by the type's own omega_exhaust (ω),
                # which is a cell parameter, not a self.config read.
                if not params and cfg.omega_exhaust > 0:
                    params = f"ω<sub>{t}</sub>"  # as written in the node formula
                label = f'|"{params}"|' if params else ""
                add(f"    {src} -.->{label} {node}", "cat")
                catalysts.setdefault(t, []).append(src)

        # interactions on rates
        for meth, role in RATE_METHODS.items():
            if meth not in m or role in ("conversion",):
                continue
            reads = m[meth]["reads"] - {t}
            for src in sorted(reads):
                params = _psym(m[meth]["by_src"].get(src, []))
                if role == "crowding":
                    if not crowding_edges:
                        continue  # already shown in the node label as "K: ..."
                    add(f'    {src} -.-o|"crowding"| {t}', "crowding", t)
                elif role == "birth":
                    add(f'    {src} ==>|"{params}"| {t}', "birth", t)
                elif role == "death":
                    add(f'    {src} --x|"{params}"| {t}', "death", t)

    # Layout only: a type nothing flows into (e.g. fibroblast) would sit in the
    # top rank, on top of the other interaction edges. Anchor it one rank below
    # the type that catalyses its conversion with an invisible link.
    for t in info:
        if t not in has_incoming and catalysts.get(t):
            add(f"    {catalysts[t][0]} ~~~ {t}", "anchor")

    colors = {"crowding": "#9aa0a6", "birth": "#2e7d32", "death": "#c62828", "cat": "#6a1b9a"}
    lines.append(f'    legend("{_legend(colors, crowding_edges)}")')

    lines += edges
    lines += [
        "    classDef tr fill:#888,stroke:#888,color:#888,width:6px,height:6px",
        "    classDef lg text-align:left,font-size:12px",
        *styles,
        "    class legend lg",
    ]
    for kind, idx in link_styles.items():
        if idx and kind in colors:
            lines.append(f"    linkStyle {','.join(map(str, idx))} stroke:{colors[kind]}")
    return "\n".join(lines) + "\n"


# ── CLI ──────────────────────────────────────────────────────────────────────

def _render_blocks(mmd: str) -> dict[Path, str]:
    out = {}
    for path in MARKDOWN_TARGETS:
        if not path.exists():
            continue
        text = path.read_text()
        if not BLOCK_RE.search(text):
            continue
        new = BLOCK_RE.sub(lambda m: m.group(1) + "```mermaid\n" + mmd + "```\n" + m.group(3), text)
        out[path] = new
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--write", action="store_true", help="update docs/model_diagram.mmd and markdown blocks")
    ap.add_argument("--check", action="store_true", help="fail if the committed diagram is stale")
    ap.add_argument("--crowding-edges", action="store_true", help="draw crowding as edges, not only in node labels")
    ap.add_argument("--svg", action="store_true", help="also render docs/model_diagram.svg with mmdc")
    args = ap.parse_args(argv)

    mmd = build_mermaid(crowding_edges=args.crowding_edges)
    blocks = _render_blocks(mmd)

    if args.check:
        stale = []
        if not MMD_PATH.exists() or MMD_PATH.read_text() != mmd:
            stale.append(MMD_PATH)
        stale += [p for p, new in blocks.items() if p.read_text() != new]
        if stale:
            print("Model diagram is stale, run `make diagram`:", *map(str, stale), sep="\n  ")
            return 1
        print("Model diagram up to date.")
        return 0

    if args.write:
        MMD_PATH.write_text(mmd)
        for p, new in blocks.items():
            p.write_text(new)
        print(f"wrote {MMD_PATH.relative_to(ROOT)}" + "".join(f", {p.relative_to(ROOT)}" for p in blocks))
        if args.svg:
            mmdc = shutil.which("mmdc")
            if mmdc is None:
                print("mmdc not found (npm i -g @mermaid-js/mermaid-cli); skipped SVG", file=sys.stderr)
            else:
                subprocess.run([mmdc, "-i", str(MMD_PATH), "-o", str(SVG_PATH), "-b", "white"], check=True)
                print(f"wrote {SVG_PATH.relative_to(ROOT)}")
        return 0

    sys.stdout.write(mmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
