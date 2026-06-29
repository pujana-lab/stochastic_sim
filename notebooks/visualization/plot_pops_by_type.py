import pandas as pd
import matplotlib.pyplot as plt

file_path = 'history.csv'

dtypes = {
    'time': 'float64',
    'type': 'category',
    'clone_id': 'category',
    'N': 'int64',
    'rb': 'float64',
    'rd': 'float64'
}

print("Loading data...")
df = pd.read_csv(file_path, dtype=dtypes)

print("Processing and aggregating populations...")
# 1. Agrupamos por tiempo y tipo para sumar las N de todos los clones en cada instante
df_aggregated = df.groupby(['time', 'type'], observed=False)['N'].sum().reset_index()

print("Plotting aggregated population dynamics...")
fig, ax = plt.subplots(figsize=(11, 6))

# 2. Ahora graficamos los datos ya consolidados por tipo
for type_name, group in df_aggregated.groupby('type', observed=False):
    # Aseguramos el orden cronológico de los eventos
    group = group.sort_values('time')
    
    # Opcional: Si sigue habiendo demasiados puntos y el guardado es lento, 
    # puedes descomentar la siguiente línea para optimizar:
    # group = group.iloc[::5]

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

output_filename = 'gillespie_aggregated_populations.png'
fig.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as '{output_filename}'")