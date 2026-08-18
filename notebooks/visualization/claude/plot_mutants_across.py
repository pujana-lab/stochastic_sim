from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


def plot_mutant_across_seeds(
    base_dir='results/multi_seed_runs',
    mutant_type_label='mutant',
    figsize=(11, 6),
    dpi=300,
):
    """Finds all seed history.csv files, extracts mutant population dynamics,

    and plots them on a single graph for cross-simulation comparison.
    """
    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64',
    }

    # Locate all history files across seed subdirectories
    history_files = sorted(Path(base_dir).glob('seed_*/history.csv'))

    if not history_files:
        raise FileNotFoundError(
            f"No history files found matching '{base_dir}/seed_*/history.csv'"
        )

    fig, ax = plt.subplots(figsize=figsize)
    aggregated_records = []

    print(f"Found {len(history_files)} simulation files. Processing...")

    for file_path in history_files:
        seed_name = file_path.parent.name  # Extracts 'seed_0001'

        # Read only required columns to improve performance
        df = pd.read_csv(
            file_path,
            dtype=dtypes,
            usecols=['time', 'type', 'N'],
        )

        # Filter strictly for mutant populations
        df_mutant = df[df['type'] == mutant_type_label]

        if df_mutant.empty:
            continue

        # Group by time step to sum all mutant clone sizes in this run
        run_agg = (
            df_mutant.groupby('time')['N']
            .sum()
            .reset_index()
            .sort_values('time')
        )
        run_agg['seed'] = seed_name

        # Plot individual simulation trajectory
        ax.plot(
            run_agg['time'],
            run_agg['N'],
            label=seed_name,
            drawstyle='steps-post',
            alpha=0.6,
            linewidth=1.2,
        )

        aggregated_records.append(run_agg)

    if not aggregated_records:
        raise ValueError(
            f"No data found for population type '{mutant_type_label}' across the simulations."
        )

    df_all_mutants = pd.concat(aggregated_records, ignore_index=True)

    # Aesthetics and formatting
    ax.set_xlabel('Time ($t$)', fontsize=12)
    ax.set_ylabel('Mutant Population Size ($N$)', fontsize=12)
    ax.set_title(
        'Mutant Population Dynamics Across Simulations',
        fontsize=14,
        fontweight='bold',
    )

    # Position legend outside if there are manageable number of seeds, or adjust styling
    if len(history_files) <= 15:
        ax.legend(
            title='Simulation',
            bbox_to_anchor=(1.02, 1),
            loc='upper left',
            fontsize=9,
        )

    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    return fig, df_all_mutants


if __name__ == '__main__':
    # Run comparison plot across all seed directories
    fig, df_mutants = plot_mutant_across_seeds(
        base_dir='results/multi_seed_runs',
        mutant_type_label='mutant',  # Change if your category label differs (e.g., 'M', 'Mutant')
    )

    # Save output
    output_path = 'results/multi_seed_runs/mutant_comparison.png'
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Comparison plot saved as '{output_path}'")