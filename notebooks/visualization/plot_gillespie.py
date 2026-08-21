import pandas as pd
import matplotlib.pyplot as plt



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
    
    Examples
    --------
    >>> fig, df = plot_gillespie_dynamics('history.parquet')
    >>> fig.savefig('my_plot.png', dpi=300, bbox_inches='tight')
    
    >>> fig, df = plot_gillespie_dynamics('history.parquet')
    >>> fig.savefig('my_plot.pdf')
    >>> plt.show()
    """
    
    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64'
    }
    
    print("Loading data...")
    df = pd.read_parquet(file_path)
    
    print("Processing and aggregating populations...")
    # Agrupamos por tiempo y tipo para sumar las N de todos los clones en cada instante
    df_aggregated = df.groupby(['time', 'type'], observed=False)['N'].sum().reset_index()
    
    print("Plotting aggregated population dynamics...")
    fig, ax = plt.subplots(figsize=figsize)
    
    # Graficamos los datos consolidados por tipo
    for type_name, group in df_aggregated.groupby('type', observed=False):
        # Aseguramos el orden cronológico de los eventos
        group = group.sort_values('time')
        
        ax.plot(
            group['time'],
            group['N'],
            label=type_name,
            drawstyle='steps-post',  # Mantiene el comportamiento de saltos estocásticos
            alpha=0.85,
            linewidth=1.5
        )
    
    # Formato estético del gráfico
    ax.set_xlabel('Time ($t$)', fontsize=12)
    ax.set_ylabel('Total Population Size ($N$)', fontsize=12)
    ax.set_title('Gillespie Simulation: Aggregated Population Dynamics', fontsize=14, fontweight='bold')
    
    ax.legend(title='Population Type', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    return fig, df_aggregated


# Usage examples
if __name__ == "__main__":
    # Example 1: Save with default settings
    fig, df = plot_gillespie_dynamics('results/multi_seed_runs/seed_0020/history.parquet')
    fig.savefig('results/multi_seed_runs/seed_0020/gillespie_aggregated_populations.png', dpi=300, bbox_inches='tight')
    print("Plot saved as 'gillespie_aggregated_populations.png'")
    
    # Example 2: Save with different format
    # fig, df = plot_gillespie_dynamics('history.parquet')
    # fig.savefig('gillespie_aggregated_populations.pdf')
    # print("Plot saved as PDF")
    
    # Example 3: Display and save
    # fig, df = plot_gillespie_dynamics('history.parquet')
    # plt.show()
    # fig.savefig('output.png')
    
    # Example 4: Access the aggregated data
    # fig, df = plot_gillespie_dynamics('history.parquet')
    # print(df.head())
    # print(df.info())