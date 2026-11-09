"""
Omega-convergence figure: Gillespie ensemble mean vs deterministic ODE solution
for the tumor-immune model from the poster.

ODE system (non-dimensionalized):
  dC/dt = C(1-C) - gamma_C * C - sigma * C * I
  dI/dt = rho * I * (1 - I - E) + (beta - omega) * C * I
  dE/dt = omega * C * I - mu_E * E

Stochastic CTMC (system size Omega):
  N_C = Omega * C, N_I = Omega * I, N_E = Omega * E
  Propensities (-> converge to ODE as Omega -> inf):
  1.  C birth (logistic):  a1 = N_C * (1 - N_C/Omega)     -> dN_C +1
  2.  C natural death:     a2 = gamma_C * N_C              -> dN_C -1
  3.  C killed by I:       a3 = sigma * N_C * N_I / Omega  -> dN_C -1
  4.  I intrinsic birth:   a4 = rho * N_I                  -> dN_I +1
  5.  I crowding death:    a5 = rho * N_I*(N_I+N_E)/Omega  -> dN_I -1
  6.  I activation by C:   a6 = beta * N_C * N_I / Omega   -> dN_I +1
  7.  I exhaustion by C:   a7 = omega * N_C * N_I / Omega  -> dN_I -1, dN_E +1
  8.  E clearance:         a8 = mu_E * N_E                 -> dN_E -1
"""

from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D

# ─── Parameters ──────────────────────────────────────────────────────────────
PARAMS = dict(
    gamma_C = 0.08,
    sigma   = 0.6,
    rho     = 0.9,
    beta    = 0.45,
    omega   = 0.25,
    mu_E    = 0.15,
)

T_MAX   = 60.0
N_TRAJ  = 60           # trajectories per Omega
OMEGAS  = [25, 150, 800]   # system sizes to compare
C0, I0, E0 = 0.06, 0.35, 0.0   # initial fractions

SEED = 42

# ─── ODE system ──────────────────────────────────────────────────────────────
def ode_rhs(t, y, p):
    C, I, E = y
    C, I, E = max(C, 0), max(I, 0), max(E, 0)
    dC = C * (1 - C) - p['gamma_C'] * C - p['sigma'] * C * I
    dI = p['rho'] * I * (1 - I - E) + (p['beta'] - p['omega']) * C * I
    dE = p['omega'] * C * I - p['mu_E'] * E
    return [dC, dI, dE]

def run_ode(p, t_eval):
    sol = solve_ivp(
        lambda t, y: ode_rhs(t, y, p),
        [0, T_MAX], [C0, I0, E0],
        t_eval=t_eval, method='RK45',
        rtol=1e-9, atol=1e-11, dense_output=False,
    )
    return sol.t, sol.y  # y[0]=C, y[1]=I, y[2]=E

# ─── Gillespie ───────────────────────────────────────────────────────────────
def run_gillespie(Omega: int, p: dict, t_eval: np.ndarray, rng: np.random.Generator):
    """Run one SSA trajectory, return interpolated fractions at t_eval."""
    gamma_C = p['gamma_C']
    sigma   = p['sigma']
    rho     = p['rho']
    beta    = p['beta']
    omega   = p['omega']
    mu_E    = p['mu_E']

    NC = int(round(C0 * Omega))
    NI = int(round(I0 * Omega))
    NE = int(round(E0 * Omega))

    t = 0.0
    t_list = [0.0]
    C_list = [NC / Omega]
    I_list = [NI / Omega]
    E_list = [NE / Omega]

    while t < T_MAX:
        # Propensities
        a1 = max(NC * (1 - NC / Omega), 0.0)   # C birth
        a2 = gamma_C * NC                        # C death
        a3 = sigma * NC * NI / Omega             # C killed
        a4 = rho * NI                            # I birth
        a5 = rho * NI * (NI + NE) / Omega       # I crowding
        a6 = beta  * NC * NI / Omega             # I activation
        a7 = omega * NC * NI / Omega             # I exhaustion
        a8 = mu_E  * NE                          # E clearance

        a0 = a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8
        if a0 <= 0:
            break

        tau = -np.log(rng.random()) / a0
        t += tau
        if t > T_MAX:
            break

        r = rng.random() * a0
        cum = 0.0
        for ai, delta in [
            (a1, ( 1,  0,  0)),
            (a2, (-1,  0,  0)),
            (a3, (-1,  0,  0)),
            (a4, ( 0,  1,  0)),
            (a5, ( 0, -1,  0)),
            (a6, ( 0,  1,  0)),
            (a7, ( 0, -1,  1)),
            (a8, ( 0,  0, -1)),
        ]:
            cum += ai
            if r < cum:
                NC += delta[0]
                NI += delta[1]
                NE += delta[2]
                break

        NC = max(NC, 0)
        NI = max(NI, 0)
        NE = max(NE, 0)

        t_list.append(t)
        C_list.append(NC / Omega)
        I_list.append(NI / Omega)
        E_list.append(NE / Omega)

    # Ensure we have the final state at T_MAX
    t_list.append(T_MAX)
    C_list.append(C_list[-1])
    I_list.append(I_list[-1])
    E_list.append(E_list[-1])

    t_arr = np.array(t_list)
    C_interp = np.interp(t_eval, t_arr, np.array(C_list))
    I_interp = np.interp(t_eval, t_arr, np.array(I_list))
    E_interp = np.interp(t_eval, t_arr, np.array(E_list))
    return C_interp, I_interp, E_interp

# ─── Run everything ───────────────────────────────────────────────────────────
t_eval = np.linspace(0, T_MAX, 400)
t_ode, y_ode = run_ode(PARAMS, t_eval)
C_ode, I_ode, E_ode = y_ode

print("ODE done. Running Gillespie ensembles...")

ensembles = {}
for Omega in OMEGAS:
    rng = np.random.default_rng(SEED + Omega)
    trajs = []
    for _ in range(N_TRAJ):
        C_t, I_t, E_t = run_gillespie(Omega, PARAMS, t_eval, rng)
        trajs.append((C_t, I_t, E_t))
    ensembles[Omega] = trajs
    print(f"  Omega={Omega}: {N_TRAJ} trajectories done")

# ─── Figure ──────────────────────────────────────────────────────────────────
COL_TUMOR  = "#D04C3E"
COL_IMMUNE = "#2A6B9C"

BG        = "#FAFAF8"
GRID_CLR  = "#E4E4E0"
TICK_CLR  = "#555550"
LABEL_CLR = "#222220"

ALPHA_TRAJ = 0.10
ALPHA_BAND = 0.18
LW_ODE     = 2.0
LW_MEAN    = 2.0
LW_TRAJ    = 0.5

fig_w = 13.0
fig_h = 4.4
fig, axes = plt.subplots(1, 3, figsize=(fig_w, fig_h),
                         facecolor=BG, constrained_layout=True)

omega_labels = {
    25:  r"$\Omega = 25$",
    150: r"$\Omega = 150$",
    800: r"$\Omega = 800$",
}

for ax, Omega in zip(axes, OMEGAS):
    trajs = ensembles[Omega]
    C_stack = np.array([tr[0] for tr in trajs])
    I_stack = np.array([tr[1] for tr in trajs])

    C_mean = C_stack.mean(axis=0)
    I_mean = I_stack.mean(axis=0)
    C_std  = C_stack.std(axis=0)
    I_std  = I_stack.std(axis=0)

    for C_t, I_t, _ in trajs:
        ax.plot(t_eval, C_t, color=COL_TUMOR,  lw=LW_TRAJ, alpha=ALPHA_TRAJ, rasterized=True)
        ax.plot(t_eval, I_t, color=COL_IMMUNE, lw=LW_TRAJ, alpha=ALPHA_TRAJ, rasterized=True)

    ax.fill_between(t_eval, C_mean - C_std, C_mean + C_std,
                    color=COL_TUMOR, alpha=ALPHA_BAND)
    ax.fill_between(t_eval, I_mean - I_std, I_mean + I_std,
                    color=COL_IMMUNE, alpha=ALPHA_BAND)

    ax.plot(t_eval, C_mean, color=COL_TUMOR,  lw=LW_MEAN, zorder=5)
    ax.plot(t_eval, I_mean, color=COL_IMMUNE, lw=LW_MEAN, zorder=5)

    ax.plot(t_ode, C_ode, color=COL_TUMOR,  lw=LW_ODE, ls='--', dashes=(5, 2.5), zorder=6)
    ax.plot(t_ode, I_ode, color=COL_IMMUNE, lw=LW_ODE, ls='--', dashes=(5, 2.5), zorder=6)

    ax.set_facecolor(BG)
    ax.set_xlim(0, T_MAX)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("Tiempo (u.a.)", fontsize=10, color=LABEL_CLR, labelpad=4)
    ax.tick_params(colors=TICK_CLR, labelsize=9)
    for spine in ax.spines.values():
        spine.set_edgecolor(GRID_CLR)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator(2))
    ax.grid(True, color=GRID_CLR, lw=0.7, zorder=0)
    ax.grid(True, which='minor', color=GRID_CLR, lw=0.3, zorder=0)
    ax.set_title(omega_labels[Omega], fontsize=13, fontweight='bold',
                 color=LABEL_CLR, pad=6)

axes[0].set_ylabel("Fracción de población", fontsize=10, color=LABEL_CLR, labelpad=6)
for ax in axes[1:]:
    ax.set_yticklabels([])

legend_handles = [
    Line2D([0], [0], color=COL_TUMOR,  lw=LW_MEAN, label="Células tumorales $C$"),
    Line2D([0], [0], color=COL_IMMUNE, lw=LW_MEAN, label="Células inmunes $I$"),
    Line2D([0], [0], color='#666662',  lw=LW_ODE, ls='--', dashes=(5, 2.5),
           label="ODE determinista"),
    Line2D([0], [0], color='#888884',  lw=3, alpha=0.55,
           label=r"Media $\pm 1\sigma$ (ensemble)"),
]
fig.legend(handles=legend_handles, loc='lower center', ncol=4,
           fontsize=9, frameon=True, framealpha=0.9,
           facecolor=BG, edgecolor=GRID_CLR,
           bbox_to_anchor=(0.5, -0.04))

fig.suptitle(
    r"Convergencia del algoritmo de Gillespie a la solución determinista con el tamaño del sistema $\Omega$",
    fontsize=11.5, color=LABEL_CLR, y=1.02,
)

out_pdf = "/mnt/user-data/outputs/omega_convergence.pdf"
out_png = "/mnt/user-data/outputs/omega_convergence.png"
fig.savefig(out_pdf, dpi=200, bbox_inches='tight', facecolor=BG)
fig.savefig(out_png, dpi=200, bbox_inches='tight', facecolor=BG)
print(f"Saved: {out_pdf}")
print(f"Saved: {out_png}")
