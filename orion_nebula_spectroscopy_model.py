#!/usr/bin/env python3
"""Simple, reproducible blackbody/Strömgren and optical-spectrum model for M42.

This is an educational one-zone model, not a substitute for a photoionization
code (e.g. Cloudy). The stellar ionizing photon rate and hydrogen Strömgren
radius are computed from first principles; ionic profiles and non-Balmer line
ratios are representative Orion-like inputs for visualization.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# ---------------------------- Physical parameters ---------------------------
T_STAR = 39_000.0                 # K, theta^1 Ori C
R_STAR = 10.0 * 6.957e10          # cm
N_H = 2_000.0                      # cm^-3 (uniform pure-H approximation)
T_E = 10_000.0                     # K, ionized gas
ALPHA_B = 2.6e-13                  # cm^3 s^-1, case B at 10^4 K
PC_CM = 3.085677581491367e18
H_PLANCK = 6.62607015e-27          # erg s
K_B = 1.380649e-16                 # erg K^-1
C_LIGHT = 2.99792458e10            # cm s^-1
EV_ERG = 1.602176634e-12

# Planck photon integral: Q = 4*pi*R^2 * 2*pi/c^2*(kT/h)^3 * int x^2/(exp(x)-1) dx
x0 = 13.6 * EV_ERG / (K_B * T_STAR)
x = np.linspace(x0, 100.0, 500_000)
planck_integrand = x*x / np.expm1(x)
photon_surface_flux = (2.0*np.pi/C_LIGHT**2) * (K_B*T_STAR/H_PLANCK)**3 * np.trapz(planck_integrand, x)
Q_H = 4.0*np.pi*R_STAR**2 * photon_surface_flux
R_STROM_CM = (3.0*Q_H/(4.0*np.pi*ALPHA_B*N_H**2))**(1.0/3.0)
R_STROM_PC = R_STROM_CM / PC_CM

# --------------------------- Radial ionization ------------------------------
r_pc = np.linspace(0.0, 1.08*R_STROM_PC, 900)
r_norm = r_pc / R_STROM_PC
# Finite, narrow front (0.6% of RS wide) avoids a numerically discontinuous edge.
h_ion = 1.0 / (1.0 + np.exp(np.clip((r_norm-1.0)/0.006, -100, 100)))
# Representative stratification: high-ionization O++ interior, O+ toward edge.
o_iii = 1.0 / (1.0 + np.exp(np.clip((r_norm-0.66)/0.085, -100, 100)))
o_ii = 1.0 - o_iii

# ---------------------- Representative optical spectrum ---------------------
# Intrinsic line fluxes relative to H-beta=1.  Balmer decrement is case-B;
# metal/He ratios are representative Orion-like values, not solved abundances.
line_data = [
    (4861.33, 1.00, r"H$\beta$"),
    (4958.91, 1.33, r"[O III] 4959"),
    (5006.84, 4.00, r"[O III] 5007"),
    (5875.62, 0.12, r"He I 5876"),
    (6548.05, 0.10, r"[N II] 6548"),
    (6562.80, 2.86, r"H$\alpha$"),
    (6583.45, 0.30, r"[N II] 6584"),
    (6716.44, 0.12, r"[S II] 6716"),
    (6730.82, 0.11, r"[S II] 6731"),
]
wave = np.linspace(4500.0, 7000.0, 30_000)
# Smooth, weak nebular continuum (arbitrary shape, scaled to 0.02 at H-beta).
continuum = 0.018 * (wave/4861.33)**-1.4 + 0.002
sigma_A = 1.15
spectrum = continuum.copy()
for lam, flux, _ in line_data:
    spectrum += flux*np.exp(-0.5*((wave-lam)/sigma_A)**2)/(sigma_A*np.sqrt(2*np.pi))
# Scale line-integrated fluxes to the requested convention H-beta=1.0.
# Spectrum is shown as flux density; annotations use intrinsic integrated ratios.

# --------------------------------- Figure ------------------------------------
plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10,
                     "figure.dpi": 150, "savefig.dpi": 220})
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12.0, 9.0), constrained_layout=True)

ax1.plot(r_pc, h_ion, color="#242f9b", lw=2.6, label=r"$n(\mathrm{H}^+)/n_H$")
ax1.plot(r_pc, o_iii, color="#d1495b", lw=2.2, label=r"$\mathrm{O}^{++}/\mathrm{O}$ (representative)")
ax1.plot(r_pc, o_ii, color="#edae49", lw=2.2, label=r"$\mathrm{O}^{+}/\mathrm{O}$ (representative)")
ax1.axvline(R_STROM_PC, color="0.25", ls="--", lw=1.3)
ax1.text(R_STROM_PC-0.012*R_STROM_PC, 0.48, rf"$R_S={R_STROM_PC:.3f}$ pc",
         rotation=90, va="center", ha="right", color="0.2")
ax1.set(xlim=(0, 1.08*R_STROM_PC), ylim=(-0.03, 1.06),
        xlabel="Radius from $\\theta^1$ Ori C (pc)", ylabel="Fraction",
        title="1. Strömgren ionization structure")
ax1.grid(alpha=0.22)
ax1.legend(loc="center left", frameon=True, fontsize=9)
ax1.text(0.985, 0.04, rf"$T_*=39{chr(44)}000$ K; $n_H={N_H:.0f}$ cm$^{{-3}}$; $Q(H^0)={Q_H:.2e}$ s$^{{-1}}$",
         transform=ax1.transAxes, ha="right", va="bottom", fontsize=9,
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.9, edgecolor="0.8"))

ax2.plot(wave, spectrum, color="#243b53", lw=0.85, label="Synthetic spectrum")
ax2.plot(wave, continuum, color="#7b8fa1", lw=1.0, ls="--", label="Nebular continuum")
# Labels are staggered vertically so the close H-alpha/[N II] and [S II] blends remain legible.
label_y = {4861.33: 0.72, 4958.91: 1.12, 5006.84: 2.08, 5875.62: 0.63,
           6548.05: 0.83, 6562.80: 1.43, 6583.45: 0.93, 6716.44: 0.66, 6730.82: 0.50}
for lam, flux, label in line_data:
    ax2.axvline(lam, color="#d1495b", lw=0.65, alpha=0.42)
    ax2.annotate(label, xy=(lam, label_y[lam]), xytext=(0, 4), textcoords="offset points",
                 rotation=90, ha="center", va="bottom", fontsize=8, color="#7d2633")
ax2.set(xlim=(4500, 7000), ylim=(0, 2.65), xlabel="Rest wavelength (Å)",
        ylabel=r"Relative flux density (H$\beta$ integrated flux = 1)",
        title="2. Synthetic Orion optical spectrum (intrinsic, unreddened)")
ax2.grid(alpha=0.20)
ax2.xaxis.set_minor_locator(AutoMinorLocator(5))
ax2.legend(loc="upper right", frameon=True, fontsize=9)
ax2.text(0.012, 0.97,
         r"Case-B H$\alpha$/H$\beta$ = 2.86; [O III] 5007/4959 = 3.01; [N II] 6584/6548 = 3.0; [S II] 6716/6731 = 1.09",
         transform=ax2.transAxes, ha="left", va="top", fontsize=8.2,
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white", alpha=0.9, edgecolor="0.8"))
fig.suptitle("Orion Nebula (M42): ionization balance and representative line spectrum", fontsize=15, weight="bold")
# Resolve output path robustly:
#   1) optional CLI argument:  python3 orion_nebula_spectroscopy_model.py [output.png]
#   2) otherwise save next to this script (works locally and in /mnt/data sandboxes)
import sys
from pathlib import Path

if len(sys.argv) > 1:
    output = Path(sys.argv[1]).expanduser()
else:
    try:
        output = Path(__file__).resolve().parent / "orion_nebula_spectroscopy_model.png"
    except NameError:  # e.g. run in an interactive/embedded context without __file__
        output = Path.cwd() / "orion_nebula_spectroscopy_model.png"

output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, bbox_inches="tight")
plt.close(fig)

print(f"Q(H0) = {Q_H:.6e} photons s^-1")
print(f"Strömgren radius = {R_STROM_PC:.6f} pc ({R_STROM_CM:.6e} cm)")
print(f"H-alpha/H-beta intrinsic = {2.86:.2f}; [S II] 6716/6731 = {0.12/0.11:.3f}")
print(f"Saved figure: {output}")
