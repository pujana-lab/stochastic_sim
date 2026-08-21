import glob
import pandas as pd
import matplotlib.pyplot as plt
## ESTE ES EL BUENO
def plot_mutant_dynamics_multi_seed(
    file_pattern='results/multi_seed_runs/seed_*/history.parquet',
    mutant_label='mutated',
    figsize=(11, 6),
    dpi=300
):
    """
    Carga y grafica la población mutante total a lo largo del tiempo
    para múltiples corridas/semillas de una simulación de Gillespie.

    Parameters
    ----------
    file_pattern : str, default='results/multi_seed_runs/seed_*/history.parquet'
        Patrón glob para encontrar los archivos history.parquet de cada corrida.
    mutant_label : str, default='mutated'
        Valor en la columna 'type' que identifica a la población mutante.
    figsize : tuple, default=(11, 6)
        Dimensiones de la figura (ancho, alto) en pulgadas.
    dpi : int, default=300
        Resolución de la figura guardada.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Objeto de la figura creada.
    df_consolidated : pd.DataFrame
        DataFrame consolidado con el tiempo, población mutante total y semilla.
    """
    files = sorted(glob.glob(file_pattern))
    if not files:
        print(f"No se encontraron archivos que coincidan con el patrón: '{file_pattern}'")
        return None, None

    print(f"Se encontraron {len(files)} archivos. Procesando...")

    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64'
    }

    fig, ax = plt.subplots(figsize=figsize)
    all_aggregated = []

    for file_path in files:
        # Extraer el identificador de la semilla (ej. 'seed_0001') desde la ruta
        path_parts = file_path.replace('\\', '/').split('/')
        seed_label = path_parts[-2] if len(path_parts) > 1 else file_path

        df = pd.read_parquet(file_path)

        # 1. Filtrar SOLO la población mutante
        df_mutant = df[df['type'] == mutant_label]

        if df_mutant.empty:
            continue

        # 2. Agrupar por tiempo para consolidar la N total de clones mutantes
        df_agg = df_mutant.groupby('time', observed=False)['N'].sum().reset_index()
        df_agg = df_agg.sort_values('time')
        df_agg['seed'] = seed_label

        all_aggregated.append(df_agg)

        # 3. Graficar la línea correspondiente a esta corrida
        ax.plot(
            df_agg['time'],
            df_agg['N'],
            label=seed_label,
            drawstyle='steps-post',  # Comportamiento estocástico por saltos
            alpha=0.7,
            linewidth=1.2
        )

    # Formato estético del gráfico
    ax.set_xlabel('Time ($t$)', fontsize=12)
    ax.set_ylabel('Total Mutant Population Size ($N_{mutant}$)', fontsize=12)
    ax.set_title('Gillespie Simulation: Total Mutant Population Trajectories', fontsize=14, fontweight='bold')
    
    # Mostrar leyenda únicamente si el número de semillas no saturará el gráfico
    if 1 <= len(files) <= 15:
        ax.legend(title='Run / Seed', bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9)
    
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    df_consolidated = pd.concat(all_aggregated, ignore_index=True) if all_aggregated else pd.DataFrame()
    return fig, df_consolidated


# Ejemplos de uso
if __name__ == "__main__":

    # Ejemplo 1: Procesar todas las semillas del patrón indicado
    fig, df_mutants = plot_mutant_dynamics_multi_seed(
        file_pattern='results/multi_seed_runs/seed_*/history.parquet'
    )
    
    if fig is not None:
        fig.savefig('mutant_trajectories_multi_seed.png', dpi=300, bbox_inches='tight')
        print("Gráfico guardado como 'mutant_trajectories_multi_seed.png'")
        plt.show()