#!/usr/bin/env python3
"""
Interactive plotting wizard for Gillespie simulation analysis.
Guides you through choosing and generating plots for mutant trajectories.

Usage:
    python plot_wizard.py
"""

import sys
from pathlib import Path
from plot_gillespie_enhanced import plot_gillespie_dynamics, plot_mutant_trajectories
from plot_multi_seed_comparison import (
    plot_all_trajectories,
    plot_mean_trajectories_with_ci,
    load_mutant_trajectories_from_seeds
)


def menu_main():
    """Main menu."""
    print("\n" + "="*70)
    print("GILLESPIE SIMULATION PLOTTING WIZARD")
    print("="*70)
    print("\nWhat would you like to do?\n")
    print("1. Plot a SINGLE simulation")
    print("2. Compare MULTIPLE simulations (different seeds)")
    print("3. Exit")
    
    choice = input("\nEnter your choice (1-3): ").strip()
    return choice


def menu_single():
    """Menu for single simulation plotting."""
    print("\n" + "-"*70)
    print("SINGLE SIMULATION PLOTTING")
    print("-"*70 + "\n")
    
    file_path = input("Enter path to history.parquet (default: results/multi_seed_runs/seed_0001/history.parquet): ").strip()
    if not file_path:
        file_path = "results/multi_seed_runs/seed_0001/history.parquet"
    
    if not Path(file_path).exists():
        print(f"✗ File not found: {file_path}")
        return
    
    print("\nWhat would you like to plot?\n")
    print("1. All populations (WT + mutants)")
    print("2. Only mutants")
    print("3. Back to main menu")
    
    choice = input("\nEnter your choice (1-3): ").strip()
    
    if choice == "1":
        print("\nGenerating plot of ALL populations...")
        fig, df = plot_gillespie_dynamics(file_path, figsize=(11, 6))
        
        output = input("Enter output filename (default: all_populations.png): ").strip()
        if not output:
            output = "all_populations.png"
        
        fig.savefig(output, dpi=300, bbox_inches='tight')
        print(f"✓ Plot saved: {output}\n")
        
    elif choice == "2":
        print("\nGenerating plot of MUTANT ONLY populations...")
        fig, df = plot_mutant_trajectories(file_path, figsize=(11, 6))
        
        if fig:
            output = input("Enter output filename (default: mutants_only.png): ").strip()
            if not output:
                output = "mutants_only.png"
            
            fig.savefig(output, dpi=300, bbox_inches='tight')
            print(f"✓ Plot saved: {output}\n")
        else:
            print("✗ No mutant populations found.\n")
    
    elif choice == "3":
        return
    
    else:
        print("✗ Invalid choice")


def menu_multiple():
    """Menu for multiple simulations comparison."""
    print("\n" + "-"*70)
    print("MULTIPLE SIMULATIONS COMPARISON")
    print("-"*70 + "\n")
    
    base_dir = input("Enter base directory (default: results/multi_seed_runs): ").strip()
    if not base_dir:
        base_dir = "results/multi_seed_runs"
    
    if not Path(base_dir).exists():
        print(f"✗ Directory not found: {base_dir}")
        return
    
    seed_dirs = list(Path(base_dir).glob("seed_*"))
    print(f"\nFound {len(seed_dirs)} seed directories")
    
    if not seed_dirs:
        print("✗ No seed directories found matching pattern 'seed_*'")
        return
    
    print("\nWhat would you like to do?\n")
    print("1. Overlay all trajectories (see variability)")
    print("2. Plot mean ± confidence interval (publication-ready)")
    print("3. Both (generate 2 plots)")
    print("4. Back to main menu")
    
    choice = input("\nEnter your choice (1-4): ").strip()
    
    if choice == "1":
        print("\nGenerating overlay plot of all trajectories...")
        alpha = input("Enter transparency level (0.1-1.0, default 0.4): ").strip()
        try:
            alpha = float(alpha) if alpha else 0.4
        except ValueError:
            alpha = 0.4
        
        fig, _ = plot_all_trajectories(base_dir=base_dir, alpha=alpha)
        
        if fig:
            output = input("Enter output filename (default: trajectories_overlay.png): ").strip()
            if not output:
                output = "trajectories_overlay.png"
            
            fig.savefig(output, dpi=300, bbox_inches='tight')
            print(f"✓ Plot saved: {output}\n")
    
    elif choice == "2":
        print("\nGenerating mean trajectory plot with confidence intervals...")
        ci = input("Enter confidence interval % (default 95): ").strip()
        try:
            ci = int(ci) if ci else 95
        except ValueError:
            ci = 95
        
        fig, stats = plot_mean_trajectories_with_ci(base_dir=base_dir, ci=ci)
        
        if fig:
            output = input("Enter output filename (default: trajectories_mean_ci.png): ").strip()
            if not output:
                output = "trajectories_mean_ci.png"
            
            fig.savefig(output, dpi=300, bbox_inches='tight')
            print(f"✓ Plot saved: {output}\n")
            
            print("Statistical summary:")
            for mut_type, data in stats.items():
                print(f"  {mut_type}: {data['n_simulations']} simulations")
    
    elif choice == "3":
        print("\nGenerating both plots...")
        
        # Plot 1: Overlay
        print("  [1/2] Overlay plot...")
        fig1, _ = plot_all_trajectories(base_dir=base_dir, alpha=0.3)
        if fig1:
            fig1.savefig("trajectories_overlay.png", dpi=300, bbox_inches='tight')
            print("       ✓ Saved: trajectories_overlay.png")
        
        # Plot 2: Mean + CI
        print("  [2/2] Mean with CI plot...")
        fig2, stats = plot_mean_trajectories_with_ci(base_dir=base_dir, ci=95)
        if fig2:
            fig2.savefig("trajectories_mean_ci.png", dpi=300, bbox_inches='tight')
            print("       ✓ Saved: trajectories_mean_ci.png\n")
    
    elif choice == "4":
        return
    
    else:
        print("✗ Invalid choice")


def main():
    """Main program loop."""
    while True:
        choice = menu_main()
        
        if choice == "1":
            menu_single()
        
        elif choice == "2":
            menu_multiple()
        
        elif choice == "3":
            print("\nGoodbye!\n")
            sys.exit(0)
        
        else:
            print("✗ Invalid choice. Please try again.")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\n✗ Interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
