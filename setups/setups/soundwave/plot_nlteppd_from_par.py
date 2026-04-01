#!/usr/bin/env python3
"""
NLTEPPD-style velocity / density / energy figures: load **physics from a FARGO .par**
(same file you passed to fargoboth) and binary dumps from ``OutputDir`` in that par.

Usage (from fargoboth repo root)::

  python setups/soundwave/plot_nlteppd_from_par.py --par setups/soundwave/soundwave.par
  python setups/soundwave/plot_nlteppd_from_par.py -p setups/soundwave/soundwave.par \\
      --project-root . --fig-dir nlteppd
  # or positional / env:
  python setups/soundwave/plot_nlteppd_from_par.py setups/soundwave/soundwave.par
  export SOUNDWAVE_PAR=setups/soundwave/soundwave.par
  python setups/soundwave/plot_nlteppd_from_par.py --project-root .

Requires: ``DT``, ``Ntot`` in the par for time base (optional; dt defaults to 0.01,
ntot defaults to number of ``gasdens*.dat`` found).
"""
import argparse
import os
import sys

_THIS = os.path.dirname(os.path.abspath(__file__))
if _THIS not in sys.path:
    sys.path.insert(0, _THIS)

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from analysis import (
    analysis_kwargs_from_par,
    eigen,
    parse_soundwave_par,
    resolve_par_path,
    _count_dust_species,
    _linear_gas_T_amplitude,
    _load_fargo_binary_field,
    _real_dtype_from_variables_par,
)


def _try_latex_rc():
    try:
        plt.rcParams.update(
            {
                "text.usetex": True,
                "font.family": "serif",
                "font.size": 16,
                "axes.labelsize": 18,
                "axes.titlesize": 20,
                "lines.linewidth": 2.5,
                "xtick.labelsize": 14,
                "ytick.labelsize": 14,
                "legend.fontsize": 14,
                "figure.figsize": (6, 5),
                "figure.dpi": 120,
                "axes.grid": False,
            }
        )
        fig, ax = plt.subplots()
        ax.set_title(r"$\mathrm{test}$")
        fig.canvas.draw()
        plt.close(fig)
    except Exception:
        plt.rcParams.update(
            {
                "font.size": 14,
                "figure.figsize": (6, 5),
                "figure.dpi": 120,
            }
        )
        print(
            "(LaTeX matplotlib backend unavailable — using default mathtext; "
            "install a TeX distro or set text.usetex False.)",
            file=sys.stderr,
        )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "-p",
        "--par",
        dest="par",
        metavar="FILE",
        help="Soundwave .par (same as for fargoboth). Else positional, else SOUNDWAVE_PAR.",
    )
    ap.add_argument(
        "par_positional",
        nargs="?",
        metavar="PAR",
        help="Optional par path if --par is not used",
    )
    ap.add_argument(
        "--project-root",
        default=".",
        help="Directory OutputDir in the par is relative to (default: cwd)",
    )
    ap.add_argument(
        "--fig-dir",
        default="nlteppd",
        help="Directory for velocity_eigen.pdf, density_eigen.pdf, energy_eigen.pdf",
    )
    args = ap.parse_args()
    par = resolve_par_path(args.par, args.par_positional)
    if not par:
        ap.error(
            "Give the par file: --par PATH, PATH as positional argument, or set SOUNDWAVE_PAR"
        )

    par_abs = os.path.abspath(par)
    project_root = os.path.abspath(args.project_root)
    kw = analysis_kwargs_from_par(par_abs)
    raw = parse_soundwave_par(par_abs)
    path = os.path.join(project_root, raw["output_dir_rel"]) + os.sep

    if not os.path.isdir(path):
        print("Output directory does not exist: %s" % path, file=sys.stderr)
        sys.exit(1)

    _try_latex_rc()
    os.makedirs(args.fig_dir, exist_ok=True)

    # --- from .par / analysis_kwargs (replaces hand-set Ts, eps, lambda_n, …) ---
    gamma = kw["gamma"]
    cpdg = kw["cpdg"]
    Ts = kw["ts"]
    Th = kw["th"]
    eps = kw["eps"]
    lambda_n = kw["lambda_complex"]
    L = kw["L"]
    rhog = kw["rhog"]
    Tg0 = kw["Tg"]
    amp = kw["amplitude"]
    tdust = kw["tdust"]
    densgr = kw["densgr"]
    densgi = kw["densgi"]
    dt = kw["dt"]
    delta_rho = densgr + 1j * densgi

    kmode = 2.0 * np.pi / L
    R_MU = 1.0
    cg = gamma * R_MU / (gamma - 1.0)
    cd = cpdg * cg
    cs = np.sqrt(gamma * Tg0 * R_MU)

    rgas = delta_rho / rhog
    vgas = 1j * lambda_n / kmode * delta_rho / rhog
    vdust = vgas / (1.0 + lambda_n * Ts)
    rdust = eps * delta_rho / rhog / (1.0 + lambda_n * Ts)
    # δT_g / T_g and δT_d / T_g-style factors (same as condinit.c)
    Tgas = _linear_gas_T_amplitude(
        lambda_n, gamma, cs, kmode, rhog, delta_rho, eps, Ts
    )
    Tdust = Tgas / (1.0 + lambda_n * Th)

    domain_z = np.loadtxt(path + "domain_z.dat")[3:-3]
    nx = 0
    x = domain_z
    x_face = float(x[nx])
    x_center = 0.5 * (float(x[nx]) + float(x[nx + 1])) if len(x) > nx + 1 else x_face

    n_frames_disk = 0
    for i in range(10000):
        if os.path.isfile(path + "gasdens%d.dat" % i):
            n_frames_disk = i + 1
        else:
            break
    ntot = kw.get("ntot", n_frames_disk)
    ntot = min(ntot, n_frames_disk)
    if ntot <= 0:
        print("No gasdens*.dat in %s" % path, file=sys.stderr)
        sys.exit(1)

    n_cells = len(x)
    real_dtype_hint = _real_dtype_from_variables_par(path.rstrip(os.sep))
    d0 = path + "gasdens0.dat"
    if real_dtype_hint is None and os.path.isfile(d0):
        nb = os.path.getsize(d0)
        if n_cells * 4 == nb:
            real_dtype_hint = np.float32
        elif n_cells * 8 == nb:
            real_dtype_hint = np.float64

    n_dust = _count_dust_species(path.rstrip(os.sep))

    time = np.arange(ntot) * dt
    time2 = np.linspace(time.min(), time.max(), ntot)

    colors = ["cornflowerblue", "orange", "gray", "green", "red"]

    # -------- velocity --------
    F3Dvgas = []
    for i in range(ntot):
        vz = _load_fargo_binary_field(
            path + "gasvz%d.dat" % i, n_cells, real_dtype_hint
        )
        F3Dvgas.append(vz[min(nx, len(vz) - 1)])
    plt.figure()
    plt.scatter(time[::10], np.array(F3Dvgas[::10]) / amp, color=colors[1], label="Gas")
    if n_dust > 0:
        for m in range(n_dust):
            F3Dvdust = []
            for i in range(ntot):
                vz = _load_fargo_binary_field(
                    path + "dust%dvz%d.dat" % (m + 1, i), n_cells, real_dtype_hint
                )
                F3Dvdust.append(vz[min(nx, len(vz) - 1)])
            plt.scatter(
                time[::10],
                np.array(F3Dvdust[::10]) / amp,
                label="Dust %d" % (m + 1),
                color=colors[0],
            )
    plt.plot(
        time2,
        eigen(vgas.real, vgas.imag, kmode, x_face, lambda_n, time2),
        color="black",
    )
    if n_dust > 0:
        plt.plot(
            time2,
            eigen(vdust.real, vdust.imag, kmode, x_face, lambda_n, time2),
            color="gray",
        )
    plt.legend(frameon=False)
    plt.ylabel("Velocity")
    plt.xlabel("Time")
    plt.xlim(0, 5)
    plt.tight_layout()
    plt.savefig(os.path.join(args.fig_dir, "velocity_eigen.pdf"))
    plt.close()

    # -------- density --------
    F3Drgas = []
    for i in range(ntot):
        rho = _load_fargo_binary_field(
            path + "gasdens%d.dat" % i, n_cells, real_dtype_hint
        )
        F3Drgas.append(rho[min(nx, len(rho) - 1)] - rhog)
    plt.figure()
    plt.scatter(time[::10], np.array(F3Drgas[::10]) / amp, label="Gas", color=colors[1])
    plt.plot(
        time2,
        eigen(rgas.real, rgas.imag, kmode, x_center, lambda_n, time2),
        color="black",
    )
    if n_dust > 0:
        m = 0
        F3Drdust = []
        for i in range(ntot):
            rho = _load_fargo_binary_field(
                path + "dust%ddens%d.dat" % (m + 1, i), n_cells, real_dtype_hint
            )
            F3Drdust.append(rho[min(nx, len(rho) - 1)] - eps * rhog)
        plt.scatter(
            time[::10],
            np.array(F3Drdust[::10]) / amp,
            label="Dust %d" % (m + 1),
            color=colors[0],
        )
        plt.plot(
            time2,
            eigen(rdust.real, rdust.imag, kmode, x_center, lambda_n, time2),
            label="Dust",
            color="gray",
        )
    plt.legend(frameon=False)
    plt.ylabel("Density")
    plt.xlabel("Time")
    plt.xlim(0, 5)
    plt.tight_layout()
    plt.savefig(os.path.join(args.fig_dir, "density_eigen.pdf"))
    plt.close()

    # -------- energy (same thermo as condinit / your notebook) --------
    e0g = cg / gamma * rhog * Tg0
    F3Degas = []
    for i in range(ntot):
        e = _load_fargo_binary_field(
            path + "gasenergy%d.dat" % i, n_cells, real_dtype_hint
        )
        F3Degas.append(e[min(nx, len(e) - 1)] - e0g)
    plt.figure()
    plt.scatter(time, np.array(F3Degas), color="orange", label="Gas")
    plt.ylabel("Energy")

    T_eigen = amp * eigen(Tgas.real, Tgas.imag, kmode, x_center, lambda_n, time2)
    rho_eigen = amp * eigen(rgas.real, rgas.imag, kmode, x_center, lambda_n, time2)
    e_eigen = (cg / gamma) * (rhog * T_eigen + Tg0 * rho_eigen)

    plt.plot(time2, e_eigen, color="black")

    if n_dust > 0:
        m = 0
        e0d = tdust * cd * eps * rhog
        F3Dedust = []
        for i in range(ntot):
            e = _load_fargo_binary_field(
                path + "dust%denergy%d.dat" % (m + 1, i), n_cells, real_dtype_hint
            )
            F3Dedust.append(e[min(nx, len(e) - 1)] - e0d)
        plt.scatter(time, np.array(F3Dedust), color=colors[0], label="Dust")
        T_eigen_dust = amp * eigen(
            Tdust.real, Tdust.imag, kmode, x_center, lambda_n, time2
        )
        rho_eigen_dust = amp * eigen(
            rdust.real, rdust.imag, kmode, x_center, lambda_n, time2
        )
        e_eigen_dust = cd * (eps * rhog * T_eigen_dust + Tg0 * rho_eigen_dust)
        plt.plot(time2, e_eigen_dust, color="gray")
    plt.xlabel("Time")
    plt.legend(frameon=False)
    plt.xlim(0, 5)
    plt.tight_layout()
    plt.savefig(os.path.join(args.fig_dir, "energy_eigen.pdf"))
    plt.close()

    print("Wrote PDFs under %s" % os.path.abspath(args.fig_dir))
    print("  par: %s" % par_abs)
    print("  data: %s" % path)
    print("  lambda = %s" % lambda_n)


if __name__ == "__main__":
    main()
