# Orion Nebula (M42) Photoionization & Synthetic Spectrum Explorer

A lightweight, reproducible educational model illustrating Strömgren sphere physics, radial ionization stratification, and representative optical emission-line spectra for the Orion Nebula (M42) ionized by the Trapezium star $\theta^1$ Ori C.

---

## 🌌 Astrophysical Background
H II regions like the Orion Nebula (M42) are prototypical astrophysical plasmas where ultraviolet radiation ($h\nu \ge 13.6\text{ eV}$) from hot, massive stars ionizes surrounding interstellar gas. 

This project explores the fundamental connection between:
- The stellar ionizing photon emission rate ($Q_{\text{H}}$) from an O-type star ($\theta^1$ Ori C, $T_\star \approx 39,000\text{ K}$).
- The balance between photoionization and volume recombination in a pure-hydrogen gas (Strömgren radius $R_S$).
- Radial ionization stratification (transition from high-ionization $\text{O}^{++}$ inner regions to outer $\text{O}^+$ frontiers).
- The resulting optical emission spectrum featuring prominent Balmer lines and forbidden transitions ([O III], [N II], [S II]).

---

## 🔬 Model Description
The accompanying Python script implements an **illustrative forward model**:
1. **Ionizing Photon Rate ($Q_{\text{H}}$):** Computes $Q_{\text{H}} \approx 7.5 \times 10^{48}\text{ s}^{-1}$ by numerical integration of the Planck blackbody flux above the hydrogen ionization threshold (13.6 eV).
2. **Strömgren Sphere:** Computes the equilibrium radius $R_S \approx 0.44\text{ pc}$ for uniform density $n_{\text{H}} = 2{,}000\text{ cm}^{-3}$ and Case-B recombination ($\alpha_B = 2.6 \times 10^{-13}\text{ cm}^3\text{ s}^{-1}$ at $T_e = 10,000\text{ K}$).
3. **Ionization Front:** Uses smooth parameterizations to visualize relative ionization fractions $x(\text{H}^+)$, $x(\text{O}^{++})$, and $x(\text{O}^+)$ across the Strömgren boundary.
4. **Synthetic Optical Spectrum:** Synthesizes representative emission lines with Gaussian profiles ($\sigma = 1.15\text{ \AA}$) and an illustrative continuum.

> **Note on Model Scope:** This is a pedagogical forward model designed for visualization and physical illustration. It is not a full photoionization code (such as Cloudy), nor does it solve empirical temperatures, densities, or chemical abundances from observed data.

---

## 🚀 How to Run

### Requirements
- Python 3.8+
- NumPy
- Matplotlib

Install dependencies:
```bash
pip install numpy matplotlib
