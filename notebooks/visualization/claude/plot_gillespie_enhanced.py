import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def plot_gillespie_dynamics(file_path='history.parquet', figsize=(11, 6), dpi=300):
    """
    Load and plot aggregated population dynamics from a Gillespie simulation.
    
    Parameters
    ----------
    file_path : str, default='history.parquet'
        Path to the parquet file containing simulation history.
    figsize : tuple, default=(11, 6)
        Figure size as (width, height) in inches.
    dpi : int, default=300
        Resolution for saving the figure.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    df_aggregated : pd.DataFrame
        The aggregated population data grouped by time and type.
    """
    
    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64'
    }
    
    print(f"Loading data from {file_path}...")
    df = pd.read_parquet(file_path)
    
    print("Processing and aggregating populations...")
    df_aggregated = df.groupby(['time', 'type'], observed=False)['N'].sum().reset_index()
    
    print("Plotting aggregated population dynamics...")
    fig, ax = plt.subplots(figsize=figsize)
    
    for type_name, group in df_aggregated.groupby('type', observed=False):
        group = group.sort_values('time')
        
        ax.plot(
            group['time'],
            group['N'],
            label=type_name,
            drawstyle='steps-post',
            alpha=0.85,
            linewidth=1.5
        )
    
    ax.set_xlabel('Time ($t$)', fontsize=12)
    ax.set_ylabel('Total Population Size ($N$)', fontsize=12)
    ax.set_title('Gillespie Simulation: Aggregated Population Dynamics', fontsize=14, fontweight='bold')
    
    ax.legend(title='Population Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    return fig, df_aggregated


def plot_mutant_trajectories(file_path='history.parquet', mutant_types=None, figsize=(11, 6), dpi=300):
    """
    Plot only mutant population trajectories from a Gillespie simulation.
    
    Parameters
    ----------
    file_path : str, default='history.parquet'
        Path to the parquet file containing simulation history.
    mutant_types : list, optional
        List of mutant type names to plot. If None, plots all types except 'WT' (wild-type).
    figsize : tuple, default=(11, 6)
        Figure size as (width, height) in inches.
    dpi : int, default=300
        Resolution for saving the figure.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The matplotlib figure object.
    df_mutants : pd.DataFrame
        The aggregated mutant population data.
    """
    
    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64'
    }
    
    print(f"Loading data from {file_path}...")
    df = pd.read_parquet(file_path)
    
    print("Processing and aggregating populations...")
    df_aggregated = df.groupby(['time', 'type'], observed=False)['N'].sum().reset_index()
    
    # Filter for mutants only
    if mutant_types is None:
        # Exclude wild-type (WT) if present, assume WT is the first type or named 'WT'
        all_types = df_aggregated['type'].unique()
        mutant_types = [t for t in all_types if t != 'WT']
    
    df_mutants = df_aggregated[df_aggregated['type'].isin(mutant_types)].copy()
    
    if df_mutants.empty:
        print("Warning: No mutant populations found!")
        return None, None
    
    print(f"Plotting mutant trajectories for types: {mutant_types}")
    fig, ax = plt.subplots(figsize=figsize)
    
    for type_name, group in df_mutants.groupby('type', observed=False):
        group = group.sort_values('time')
        
        ax.plot(
            group['time'],
            group['N'],
            label=type_name,
            drawstyle='steps-post',
            alpha=0.85,
            linewidth=1.5
        )
    
    ax.set_xlabel('Time ($t$)', fontsize=12)
    ax.set_ylabel('Total Mutant Population Size ($N$)', fontsize=12)
    ax.set_title('Gillespie Simulation: Mutant Population Trajectories', fontsize=14, fontweight='bold')
    
    ax.legend(title='Mutant Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    return fig, df_mutants


# Usage examples
if __name__ == "__main__":
    # Example 1: Plot all populations
    fig1, df1 = plot_gillespie_dynamics('results/multi_seed_runs/seed_0001/history.parquet')
    fig1.savefig('results/multi_seed_runs/seed_0001/gillespie_all_populations.png', dpi=300, bbox_inches='tight')
    print("All populations plot saved\n")
    
    # Example 2: Plot only mutants
    fig2, df2 = plot_mutant_trajectories('results/multi_seed_runs/seed_0001/history.parquet')
    fig2.savefig('results/multi_seed_runs/seed_0001/gillespie_mutant_trajectories.png', dpi=300, bbox_inches='tight')
    print("Mutant trajectories plot saved")
