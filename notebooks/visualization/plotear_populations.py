import pandas as pd
import matplotlib.pyplot as plt

file_path = "history.csv"
#TODO: same as other: COnvertir en funcion para que me saque solo la figura y luego llamarla desde fuera y guardarlo como quiera 

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

print("Plotting population dynamics...")
fig, ax = plt.subplots(figsize=(10, 6))

# Definimos colores específicos para las poblaciones clave
special_colors = {
    'root': '#e377c2',          # Rosa
    'mutated_root': '#8c564b',  # Marrón
    'immune': '#9467bd',        # Morado
    'exhausted': '#d62728'      # Rojo
}
mutant_color = '#7f7f7f'        # Gris para el resto de mutantes numéricas

mutant_legend_added = False

for type_name, group in df.groupby('clone_id'):
    group = group.sort_values('time')
    
    # Determinar color y etiqueta
    if type_name in special_colors:
        color = special_colors[type_name]
        label = type_name
    else:
        color = mutant_color
        # Solo añadimos la etiqueta una vez para no duplicarla en la leyenda
        if not mutant_legend_added:
            label = 'Mutant Clones (Numerical)'
            mutant_legend_added = True
        else:
            label = None # Al ser None, matplotlib lo ignora en la leyenda

    ax.plot(
        group['time'], 
        group['N'], 
        label=label, 
        color=color,
        drawstyle='steps-post',  
        alpha=0.75, 
        linewidth=1.5
    )

ax.set_xlabel('Time ($t$)', fontsize=12)
ax.set_ylabel('Population Size ($N$)', fontsize=12)
ax.set_title('Gillespie Simulation: Population Dynamics Over Time', fontsize=14, fontweight='bold')

# La leyenda ahora será compacta y limpia
ax.legend(title='Population Type', loc='upper left', fontsize=10, frameon=True)
ax.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
output_filename = 'gillespie_populations_grouped.png'
fig.savefig(output_filename, dpi=300, bbox_inches='tight')
print(f"Plot successfully saved as '{output_filename}'")