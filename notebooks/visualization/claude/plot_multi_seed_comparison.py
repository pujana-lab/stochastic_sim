import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


NON_MUTANT_TYPES = {'WT', 'root', 'mutated_root', 'immune', 'exhausted', 'wildtype'}


def _is_mutant_row(row):
    """Keep only mutant rows; ignore the non-mutant control populations."""
    type_name = str(row.get('type', '')).strip()
    clone_id = str(row.get('clone_id', '')).strip()
    return type_name not in NON_MUTANT_TYPES and clone_id not in NON_MUTANT_TYPES


def load_mutant_trajectories_from_seeds(base_dir='results/multi_seed_runs', 
                                        mutant_types=None,
                                        pattern='seed_*'):
    """
    Load mutant trajectories from multiple seed runs.
    
    Parameters
    ----------
    base_dir : str, default='results/multi_seed_runs'
        Base directory containing seed subdirectories.
    mutant_types : list, optional
        List of mutant type names to extract. If None, extracts all non-WT types.
    pattern : str, default='seed_*'
        Pattern to match seed directories.
    
    Returns
    -------
    trajectories : dict
        Dictionary with structure {mutant_type: {seed: data_array}}
    seed_dirs : list
        List of found seed directories (sorted).
    """
    
    base_path = Path(base_dir)
    seed_dirs = sorted(base_path.glob(pattern))
    
    print(f"Found {len(seed_dirs)} seed directories")
    
    dtypes = {
        'time': 'float64',
        'type': 'category',
        'clone_id': 'category',
        'N': 'int64',
        'rb': 'float64',
        'rd': 'float64'
    }
    
    trajectories = {}
    all_times = None
    
    for seed_dir in seed_dirs:
        history_file = seed_dir / 'history.csv'
        seed_name = seed_dir.name
        
        if not history_file.exists():
            print(f"  ⚠ Skipping {seed_name}: history.csv not found")
            continue
        
        try:
            print(f"  Loading {seed_name}...")
            df = pd.read_csv(history_file, dtype=dtypes)

            # Keep only mutant rows. This excludes WT/control populations and the
            # main root/immune/exhausted families that are not mutants.
            mutant_df = df[df.apply(_is_mutant_row, axis=1)].copy()
            if mutant_df.empty:
                print(f"  ⚠ No mutant data in {seed_name}")
                continue

            # Aggregate by time and type
            df_agg = mutant_df.groupby(['time', 'type'], observed=False)['N'].sum().reset_index()

            # Determine mutant types if not specified
            if mutant_types is None and not trajectories:
                all_types = df_agg['type'].unique()
                mutant_types = sorted([t for t in all_types if str(t) not in NON_MUTANT_TYPES])
                print(f"  Detected mutant types: {mutant_types}\n")

            # Store data for each mutant type
            for mut_type in (mutant_types or []):
                if mut_type not in trajectories:
                    trajectories[mut_type] = {}

                type_data = df_agg[df_agg['type'] == mut_type].sort_values('time')

                if not type_data.empty:
                    trajectories[mut_type][seed_name] = type_data[['time', 'N']].values

                    # Keep track of all unique time points for later interpolation
                    if all_times is None:
                        all_times = type_data['time'].values
                    else:
                        all_times = np.unique(np.concatenate([all_times, type_data['time'].values]))

        except Exception as e:
            print(f"  Error loading {seed_name}: {e}")
            continue
    
    return trajectories, seed_dirs, sorted(mutant_types) if mutant_types else []


def plot_all_trajectories(base_dir='results/multi_seed_runs',
                          mutant_types=None,
                          figsize=(12, 7),
                          alpha=0.4,
                          dpi=300):
    """
    Plot all mutant trajectories from all seed runs on same figure.
    Each seed's trajectory is a semi-transparent line.
    
    Parameters
    ----------
    base_dir : str, default='results/multi_seed_runs'
        Base directory containing seed subdirectories.
    mutant_types : list, optional
        List of mutant types to plot. If None, auto-detect.
    figsize : tuple, default=(12, 7)
        Figure size as (width, height) in inches.
    alpha : float, default=0.4
        Alpha transparency for individual trajectories.
    dpi : int, default=300
        Resolution for saving.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    trajectories : dict
        The trajectory data.
    """
    
    trajectories, seed_dirs, detected_types = load_mutant_trajectories_from_seeds(
        base_dir=base_dir,
        mutant_types=mutant_types
    )
    
    if not trajectories:
        print("No data to plot!")
        return None, None
    
    mutant_types = detected_types if mutant_types is None else mutant_types
    
    fig, axes = plt.subplots(
        nrows=len(mutant_types),
        ncols=1,
        figsize=(figsize[0], figsize[1] * len(mutant_types) / 2),
        sharex=True
    )
    
    # Handle single subplot case
    if len(mutant_types) == 1:
        axes = [axes]
    
    print(f"\nPlotting trajectories for {len(mutant_types)} mutant type(s)...")
    
    for ax, mut_type in zip(axes, mutant_types):
        if mut_type not in trajectories:
            print(f"  Skipping {mut_type}: no data")
            continue
        
        type_trajectories = trajectories[mut_type]
        print(f"  Plotting {mut_type} ({len(type_trajectories)} trajectories)...")
        
        # Plot each seed's trajectory
        for seed_name, data in type_trajectories.items():
            times = data[:, 0]
            populations = data[:, 1]
            
            ax.plot(
                times,
                populations,
                label=seed_name,
                drawstyle='steps-post',
                alpha=alpha,
                linewidth=1.0,
                color='steelblue'
            )
        
        ax.set_ylabel(f'{mut_type}\n($N$)', fontsize=11)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_yscale('log')  # Optional: use log scale to see small populations
    
    axes[-1].set_xlabel('Time ($t$)', fontsize=12)
    fig.suptitle(
        f'Mutant Population Trajectories Across {len(seed_dirs)} Simulations',
        fontsize=14,
        fontweight='bold',
        y=0.995
    )
    
    plt.tight_layout()
    
    return fig, trajectories


def plot_mean_trajectories_with_ci(base_dir='results/multi_seed_runs',
                                    mutant_types=None,
                                    figsize=(12, 7),
                                    ci=95,
                                    dpi=300):
    """
    Plot mean trajectories with confidence intervals (shaded regions).
    
    Parameters
    ----------
    base_dir : str, default='results/multi_seed_runs'
        Base directory containing seed subdirectories.
    mutant_types : list, optional
        List of mutant types to plot. If None, auto-detect.
    figsize : tuple, default=(12, 7)
        Figure size as (width, height) in inches.
    ci : float, default=95
        Confidence interval percentage (e.g., 95 for 95% CI).
    dpi : int, default=300
        Resolution for saving.
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object.
    stats : dict
        Dictionary with mean and CI data for each type.
    """
    
    trajectories, seed_dirs, detected_types = load_mutant_trajectories_from_seeds(
        base_dir=base_dir,
        mutant_types=mutant_types
    )
    
    if not trajectories:
        print("No data to plot!")
        return None, None
    
    mutant_types = detected_types if mutant_types is None else mutant_types
    
    # Create a common time grid by finding all unique time points
    all_times_set = set()
    for mut_type in mutant_types:
        for seed_traj in trajectories[mut_type].values():
            all_times_set.update(seed_traj[:, 0])
    
    common_times = np.array(sorted(all_times_set))
    
    stats = {}
    
    fig, axes = plt.subplots(
        nrows=len(mutant_types),
        ncols=1,
        figsize=(figsize[0], figsize[1] * len(mutant_types) / 2),
        sharex=True
    )
    
    if len(mutant_types) == 1:
        axes = [axes]
    
    print(f"\nComputing statistics for {len(mutant_types)} mutant type(s)...")
    
    for ax, mut_type in zip(axes, mutant_types):
        if mut_type not in trajectories:
            continue
        
        type_trajectories = trajectories[mut_type]
        
        # Interpolate all trajectories to common time grid
        populations_at_times = []
        
        for seed_name, data in type_trajectories.items():
            times = data[:, 0]
            pops = data[:, 1]
            
            # Interpolate using step-post (forward-fill)
            interp_pops = np.interp(common_times, times, pops, left=np.nan, right=pops[-1])
            populations_at_times.append(interp_pops)
        
        populations_array = np.array(populations_at_times)
        
        # Compute mean and CI
        mean_pop = np.nanmean(populations_array, axis=0)
        
        # Compute percentiles for CI
        alpha_ci = (100 - ci) / 2
        lower_ci = np.nanpercentile(populations_array, alpha_ci, axis=0)
        upper_ci = np.nanpercentile(populations_array, 100 - alpha_ci, axis=0)
        
        stats[mut_type] = {
            'times': common_times,
            'mean': mean_pop,
            'lower_ci': lower_ci,
            'upper_ci': upper_ci,
            'n_simulations': len(type_trajectories)
        }
        
        print(f"  {mut_type}: mean computed from {len(type_trajectories)} simulations")
        
        # Plot mean with confidence interval
        ax.plot(
            common_times,
            mean_pop,
            label=f'{mut_type} (mean)',
            drawstyle='steps-post',
            linewidth=2.0,
            color='darkred'
        )
        
        ax.fill_between(
            common_times,
            lower_ci,
            upper_ci,
            step='post',
            alpha=0.3,
            color='darkred',
            label=f'{ci}% CI'
        )
        
        ax.set_ylabel(f'{mut_type}\n($N$)', fontsize=11)
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_yscale('log')
        ax.legend(loc='best', fontsize=9)
    
    axes[-1].set_xlabel('Time ($t$)', fontsize=12)
    fig.suptitle(
        f'Mutant Trajectories: Mean ± {ci}% CI ({len(seed_dirs)} simulations)',
        fontsize=14,
        fontweight='bold',
        y=0.995
    )
    
    plt.tight_layout()
    
    return fig, stats


# Usage examples
if __name__ == "__main__":
    print("="*70)
    print("MULTI-SEED MUTANT TRAJECTORY ANALYSIS")
    print("="*70)
    
    # Example 1: Plot all individual trajectories
    print("\n[1/2] Plotting all individual trajectories...")
    fig1, traj = plot_all_trajectories(
        base_dir='results/multi_seed_runs',
        alpha=0.3
    )
    if fig1:
        fig1.savefig('results/all_mutant_trajectories.png', dpi=300, bbox_inches='tight')
        print("✓ Saved: results/all_mutant_trajectories.png\n")
    
    # Example 2: Plot mean with confidence intervals
    print("[2/2] Plotting mean trajectories with confidence intervals...")
    fig2, stats = plot_mean_trajectories_with_ci(
        base_dir='results/multi_seed_runs',
        ci=95
    )
    if fig2:
        fig2.savefig('results/mutant_trajectories_mean_ci.png', dpi=300, bbox_inches='tight')
        print("✓ Saved: results/mutant_trajectories_mean_ci.png")
