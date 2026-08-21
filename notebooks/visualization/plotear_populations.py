import matplotlib.pyplot as plt
import pandas as pd

## ESTE ES EL BUENO
def plot_clone_dynamics(
    file_path="results/multi_seed_runs/seed_0001",
    populations_to_include=None,  # e.g. ['root', 'mutated_clones'] or None for all
    figsize=(10, 6),
    plot_extinct: bool = False,
    log_y: bool = False,
):
    """Plots Gillespie simulation dynamics by clone_id with dynamic population filtering.

    Supports 'mutated_clones' or 'mutant_clones' keyword to include all numerical mutants.
    """
    dtypes = {
        "time": "float64",
        "type": "category",
        "clone_id": "category",
        "N": "int64",
        "rb": "float64",
        "rd": "float64",
    }

    df = pd.read_parquet(f"{file_path}/history.parquet")
    df_extinct = pd.read_parquet(f"{file_path}/clones.parquet")

    # Identify clones with no current population.
    extinct_clone_ids = set(
        df_extinct.loc[df_extinct["N"].eq(0), "clone_id"]
    )

    special_colors = {
        "root": "#e377c2",  # Pink
        "mutated_root": "#7f7f7f",  # Brown
        "immune": "#9467bd",  # Purple
        "exhausted": "#C0A822",  # Red
    }

    # Dynamic filtering logic
    if populations_to_include is not None:
        # Convert to set for faster lookup and normalise string matching
        targets = set(populations_to_include)

        # Check if user requested all numerical/non-special mutant clones
        include_all_mutants = bool(
            targets & {"mutated_clones", "mutant_clones"}
        )

        def should_keep(clone_id):
            # 1. Directly matched special population (e.g., 'root', 'exhausted')
            if clone_id in targets:
                return True
            # 2. Mutant clone (not in special_colors) when mutant wildcard flag is active
            if include_all_mutants and clone_id not in special_colors:
                return True
            return False

        # Apply filtering across unique clone IDs
        df = df[df["clone_id"].apply(should_keep)]

    fig, ax = plt.subplots(figsize=figsize)

    mutant_color = "#7f7f7f"
    extinct_color = "#d62727ff"
    mutant_legend_added = False
    extinct_legend_added = False

    for clone_name, group in df.groupby("clone_id", observed=True):
        if group.empty:
            continue

        group = group.sort_values("time")

        if plot_extinct and clone_name in extinct_clone_ids:
            color = extinct_color
            line_zorder = 3
            line_width = 1.8
            if not extinct_legend_added:
                label = "Extinct Clones"
                extinct_legend_added = True
            else:
                label = None
        elif clone_name in special_colors and clone_name not in {"mutated_root"}:
            color = special_colors[clone_name]
            label = clone_name
            line_zorder = 2
            line_width = 1.5
        else:
            color = mutant_color
            line_zorder = 2
            line_width = 1.5
            if not mutant_legend_added:
                label = "Mutant Clones"
                mutant_legend_added = True
            else:
                label = None

        ax.plot(
            group["time"],
            group["N"],
            label=label,
            color=color,
            drawstyle="steps-post",
            alpha=0.75,
            linewidth=line_width,
            zorder=line_zorder,
        )

    ax.set_xlabel("Time ($t$)", fontsize=12)
    ax.set_ylabel("Population Size ($N$)", fontsize=12)
    if log_y:
        # Note: trajectories with N=0 are not shown on a strict log scale.
        ax.set_yscale("log")
    ax.set_title(
        "Gillespie Simulation: Population Dynamics Over Time",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(title="Population Type", loc="upper left", frameon=True)
    ax.grid(True, linestyle="--", alpha=0.5)
    # ax.set_ylim(0,400)
    plt.tight_layout()
    return fig


# Example Usage:
if __name__ == "__main__":
    # Plot 'root' AND all numerical mutant clones, but exclude 'immune' and 'exhausted'
    fig = plot_clone_dynamics(
        "results/multi_seed_runs/seed_0060",
        populations_to_include= ["mutated_root","mutated_clones"],
        plot_extinct=True,
        log_y=True,
    )
    fig.savefig("results/multi_seed_runs/seed_0060/mutated_grouped_logarithmic.png", dpi=300, bbox_inches="tight")