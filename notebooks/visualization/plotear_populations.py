import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp


def system_ode(t, x, p):
    """Sistema de Ecuaciones Diferenciales Adimensional."""
    C, T, E = x
    dC = C * (1.0 - C) - p["gamma"] * C - p["sigma"] * C * T
    dT = p["rho"] * T * (1.0 - T - E) + (p["beta"] - p["omega"]) * C * T
    dE = p["omega"] * C * T - p["mu"] * E
    return [dC, dT, dE]


def plot_clone_dynamics(
    file_path="results/multi_seed_runs/seed_0001",
    populations_to_include=None,
    figsize=(10, 6),
    plot_extinct: bool = False,
    log_y: bool = False,
    plot_ode: bool = True,
    system_size: float = 1000.0,
):
    """Plots Gillespie simulation dynamics with optional ODE overlay."""
    df = pd.read_parquet(f"{file_path}/history.parquet")
    df_extinct = pd.read_parquet(f"{file_path}/clones.parquet")

    extinct_clone_ids = set(df_extinct.loc[df_extinct["N"].eq(0), "clone_id"])

    special_colors = {
        "root": "#e377c2",  # Pink (Tumor / C)
        "mutated_root": "#7f7f7f",  # Brown
        "immune": "#9467bd",  # Purple (Efectoras / T)
        "exhausted": "#C0A822",  # Gold/Yellow (Exhaustas / E)
    }

    if populations_to_include is not None:
        targets = set(populations_to_include)
        include_all_mutants = bool(
            targets & {"mutated_clones", "mutant_clones"}
        )

        def should_keep(clone_id):
            if clone_id in targets:
                return True
            if include_all_mutants and clone_id not in special_colors:
                return True
            return False

        df = df[df["clone_id"].apply(should_keep)]

    fig, ax = plt.subplots(figsize=figsize)

    mutant_color = "#7f7f7f"
    extinct_color = "#d62727ff"
    mutant_legend_added = False
    extinct_legend_added = False

    # 1. Graficar las trayectorias estocásticas de Gillespie
    for clone_name, group in df.groupby("clone_id", observed=True):
        if group.empty:
            continue

        group = group.sort_values("time")

        if plot_extinct and clone_name in extinct_clone_ids:
            color = extinct_color
            line_zorder = 3
            line_width = 1.8
            label = "Extinct Clones" if not extinct_legend_added else None
            extinct_legend_added = True
        elif (
            clone_name in special_colors and clone_name not in {"mutated_root"}
        ):
            color = special_colors[clone_name]
            label = f"Gillespie: {clone_name}"
            line_zorder = 2
            line_width = 1.5
        else:
            color = mutant_color
            line_zorder = 2
            line_width = 1.5
            label = "Mutant Clones" if not mutant_legend_added else None
            mutant_legend_added = True

        ax.plot(
            group["time"],
            group["N"],
            label=label,
            color=color,
            drawstyle="steps-post",
            alpha=0.6,
            linewidth=line_width,
            zorder=line_zorder,
        )

    # 2. Integrar y Graficar la Solución Determinista (ODE)
    if plot_ode:
        t_max = df["time"].max() if not df.empty else 200.0
        t_eval = np.linspace(0, t_max, 1000)

        # Parámetros exactos indicados
        p = {
            "rho": 0.1,
            "sigma": 0.8,
            "mu": 0.4,
            "beta": 0.45,
            "gamma": 0.3,
            "omega": 0.9,
        }
        x0 = [0.1, 1.2, 0.0]  # [C0, T0, E0] adimensionales

        # Resolver ODE
        sol = solve_ivp(
            system_ode,
            (0, t_max),
            x0,
            args=(p,),
            t_eval=t_eval,
            method="RK45",
        )

        # Mapeo de soluciones escaladas por el tamaño de sistema N_scale
        ode_curves = {
            "root": sol.y[0] * system_size,  # C -> Tumor
            "immune": sol.y[1] * system_size,  # T -> Efectoras
            "exhausted": sol.y[2] * system_size,  # E -> Exhaustas
        }

        # Dibujar líneas discontinuas negras/oscuras para comparar
        for pop_key, y_values in ode_curves.items():
            if (
                populations_to_include is None
                or pop_key in populations_to_include
            ):
                color = special_colors.get(pop_key, "black")
                ax.plot(
                    sol.t,
                    y_values,
                    linestyle="--",
                    color=color,
                    linewidth=2.2,
                    alpha=0.9,
                    label=f"ODE: {pop_key}",
                    zorder=5,
                )

    ax.set_xlabel("Time ($t$)", fontsize=12)
    ax.set_ylabel("Population Size ($N$)", fontsize=12)
    if log_y:
        ax.set_yscale("log")
    ax.set_title(
        "Gillespie vs ODE Deterministic Dynamics",
        fontsize=14,
        fontweight="bold",
    )
    ax.legend(title="Population Type", loc="upper right", frameon=True)
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    return fig


if __name__ == "__main__":
    fig = plot_clone_dynamics(
        "results/gamma_omega_sweep/roj_verde_2",
        populations_to_include=None,
        plot_extinct=False,
        log_y=False,
        plot_ode=True,
        system_size=1000.0,
    )
    fig.savefig(
        "results/gamma_omega_sweep/roj_verde_2/gillespie_vs_ode.svg",
        dpi=300,
        bbox_inches="tight",
        format="svg",
    )