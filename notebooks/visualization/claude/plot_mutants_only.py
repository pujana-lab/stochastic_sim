"""
Simple script to plot only mutant trajectories from a single Gillespie simulation.

Usage:
    python plot_mutants_only.py
    
    Then modify the file_path in the main section as needed.
"""

import pandas as pd
import matplotlib.pyplot as plt
from plot_gillespie_enhanced import plot_mutant_trajectories


if __name__ == "__main__":
    # Configure these paths according to your needs
    file_path = 'results/multi_seed_runs/seed_0001/history.csv'
    output_path = 'results/multi_seed_runs/seed_0001/mutants_only.png'
    
    print(f"Processing: {file_path}")
    print("-" * 60)
    
    # Plot mutant trajectories
    fig, df_mutants = plot_mutant_trajectories(
        file_path=file_path,
        mutant_types=None,  # None = all non-WT types
        figsize=(11, 6),
        dpi=300
    )
    
    if fig:
        print(f"\nSaving to: {output_path}")
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print("✓ Plot saved successfully!")
        
        print(f"\nMutant population summary:")
        print(df_mutants.groupby('type')['N'].describe())
    else:
        print("✗ Failed to create plot")
