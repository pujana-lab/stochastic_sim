#!/usr/bin/env python3
"""Generate clone recount and immune-escape/extinction stats from a seed folder.

Usage:
    /home/vmanso/PHD/moran/.venv/bin/python scripts/clone_recount_report.py results/multi_seed_runs/seed_0060
"""

from __future__ import annotations

import argparse
import csv
import os
import statistics
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


SPECIAL_CLONES = {"root", "mutated_root", "immune", "exhausted"}


@dataclass
class CloneRow:
    clone_id: str
    parent: str
    n: int


@dataclass
class HistoryRow:
    time: float
    clone_id: str
    n: int


def is_numeric_clone_id(clone_id: str) -> bool:
    return clone_id.isdigit()


def is_recount_clone_id(clone_id: str) -> bool:
    return clone_id == "mutated_root" or is_numeric_clone_id(clone_id)


def clone_rank_key(clone: CloneRow) -> Tuple[int, int, str]:
    # Keep numeric ordering for numeric IDs and place non-numeric tracked IDs after.
    if clone.clone_id.isdigit():
        return (-clone.n, 0, str(int(clone.clone_id)))
    return (-clone.n, 1, clone.clone_id)


def read_clones_csv(path: str) -> List[CloneRow]:
    rows: List[CloneRow] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            clone_id = str(row["clone_id"]).strip()
            parent = str(row.get("parent", "") or "").strip()
            n = int(float(row["N"]))
            rows.append(CloneRow(clone_id=clone_id, parent=parent, n=n))
    return rows


def read_history_csv(path: str) -> List[HistoryRow]:
    rows: List[HistoryRow] = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(
                HistoryRow(
                    time=float(row["time"]),
                    clone_id=str(row["clone_id"]).strip(),
                    n=int(float(row["N"])),
                )
            )
    return rows


def final_snapshot(history_rows: List[HistoryRow]) -> Tuple[float, List[HistoryRow]]:
    if not history_rows:
        return 0.0, []
    max_time = max(r.time for r in history_rows)
    final_rows = [r for r in history_rows if r.time == max_time]
    return max_time, final_rows


def immune_first_zero_time(history_rows: List[HistoryRow]) -> Optional[float]:
    times = [r.time for r in history_rows if r.clone_id == "immune" and r.n == 0]
    if not times:
        return None
    return min(times)


def pct(part: float, whole: float) -> float:
    if whole == 0:
        return 0.0
    return 100.0 * part / whole


def format_float(value: Optional[float], decimals: int = 2) -> str:
    if value is None:
        return "N/A"
    return f"{value:.{decimals}f}"


def make_markdown_table(title: str, rows: List[Tuple[str, str]]) -> str:
    lines = [f"## {title}", "| Metric | Value |", "|---|---:|"]
    for k, v in rows:
        lines.append(f"| {k} | {v} |")
    return "\n".join(lines)


def analyze(seed_dir: str, top_n: int) -> str:
    clones_path = os.path.join(seed_dir, "clones.csv")
    history_path = os.path.join(seed_dir, "history.csv")

    if not os.path.exists(clones_path):
        raise FileNotFoundError(f"Not found: {clones_path}")
    if not os.path.exists(history_path):
        raise FileNotFoundError(f"Not found: {history_path}")

    clone_rows = read_clones_csv(clones_path)
    history_rows = read_history_csv(history_path)

    numeric_clones = [r for r in clone_rows if is_recount_clone_id(r.clone_id)]
    numeric_total = len(numeric_clones)
    numeric_extinct = sum(1 for r in numeric_clones if r.n == 0)
    numeric_survivors = numeric_total - numeric_extinct

    survivor_rows = [r for r in numeric_clones if r.n > 0]
    survivor_rows_sorted = sorted(
        survivor_rows,
        key=clone_rank_key,
    )

    survivor_cells = sum(r.n for r in survivor_rows)
    survivor_mean = (survivor_cells / numeric_survivors) if numeric_survivors else 0.0
    survivor_median = statistics.median([r.n for r in survivor_rows]) if survivor_rows else 0.0

    top1_cells = survivor_rows_sorted[0].n if survivor_rows_sorted else 0
    topk_cells = sum(r.n for r in survivor_rows_sorted[:top_n])

    max_time, final_rows = final_snapshot(history_rows)
    final_by_id: Dict[str, int] = {r.clone_id: r.n for r in final_rows}

    immune_final = final_by_id.get("immune", 0)
    root_final = final_by_id.get("root", 0)
    mutated_root_final = final_by_id.get("mutated_root", 0)
    exhausted_final = final_by_id.get("exhausted", 0)

    numeric_cells_final_from_history = sum(
        n for cid, n in final_by_id.items() if is_recount_clone_id(cid)
    )
    tumor_cells_final = sum(
        n for cid, n in final_by_id.items() if cid not in {"immune", "exhausted"}
    )

    immune_zero_t = immune_first_zero_time(history_rows)

    if immune_final > 0:
        numeric_vs_immune_ratio = f"{numeric_cells_final_from_history / immune_final:.3f}"
    elif numeric_cells_final_from_history > 0:
        numeric_vs_immune_ratio = "inf"
    else:
        numeric_vs_immune_ratio = "0.000"

    summary_rows = [
        ("Seed", os.path.basename(os.path.normpath(seed_dir))),
        ("Clones created (numeric + mutated_root)", str(numeric_total)),
        ("Clones extinct (final N=0)", str(numeric_extinct)),
        ("Clones escaped", str(numeric_survivors)),
        ("% extinct", f"{pct(numeric_extinct, numeric_total):.2f}%"),
        ("% non-extinct (escape at final time)", f"{pct(numeric_survivors, numeric_total):.2f}%"),
        ("Immune time to extinction", format_float(immune_zero_t, 4) if immune_zero_t is not None else "never"),
        ("Final time", f"{max_time:.4f}"),
    ]

    burden_rows = [
        ("Immune final", str(immune_final)),
        ("Root final", str(root_final)),
        ("Mutated root final", str(mutated_root_final)),
        ("Exhausted final", str(exhausted_final)),
        ("Final recount cells (numeric + mutated_root)", str(numeric_cells_final_from_history)),
        ("Final tumor burden (excluding immune/exhausted)", str(tumor_cells_final)),
        ("Final recount/immune ratio", numeric_vs_immune_ratio),
    ]

    dominance_rows = [
        ("Recount survivors", str(numeric_survivors)),
        ("Cells in survivors", str(survivor_cells)),
        ("Mean N in survivors", f"{survivor_mean:.2f}"),
        ("Median N in survivors", f"{survivor_median:.2f}"),
        (f"Top-1 % of survivor burden", f"{pct(top1_cells, survivor_cells):.2f}%"),
        (f"Top-{top_n} % of survivor burden", f"{pct(topk_cells, survivor_cells):.2f}%"),
    ]

    top_rows = survivor_rows_sorted[:top_n]
    top_lines = [f"## Top {top_n} surviving clones", "| clone_id | Final N |", "|---:|---:|"]
    if top_rows:
        for row in top_rows:
            top_lines.append(f"| {row.clone_id} | {row.n} |")
    else:
        top_lines.append("| - | 0 |")

    blocks = [
        make_markdown_table("Summary", summary_rows),
        make_markdown_table("Burden And Immune Control", burden_rows),
        make_markdown_table("Clonal Dominance", dominance_rows),
        "\n".join(top_lines),
    ]
    return "\n\n".join(blocks)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Numeric clone recount plus immune-escape/extinction statistics.",
    )
    parser.add_argument(
        "seed_dir",
        help="Directory containing clones.csv and history.csv (e.g. results/multi_seed_runs/seed_0060)",
    )
    parser.add_argument(
        "--top",
        type=int,
        default=5,
        help="Number of top clones to show (default: 5)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    report = analyze(args.seed_dir, top_n=max(1, args.top))
    print(report)


if __name__ == "__main__":
    main()
