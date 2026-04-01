"""
Soundwave post-processing: read a FARGO ``.par``, load ``OutputDir`` dumps, optional PDF.

Notebook-style loops over ``gasdens{i}``, ``gasvz{i}``, ``gasenergy{i}``, and when present
``dust*dens``, ``dust*vz``, ``dust*energy``. Physics (Γ, CPDG, Ts, Th, ε, λ, amplitude, L,
TDUST, …) comes from the same ``.par`` via ``analysis_kwargs_from_par``.
"""
import os
import sys
import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_PLT = True
except ImportError:
    HAS_PLT = False


# ---------------------------------------------------------------------------
# .par → physics (replace hand-set Ts, Th, eps, lambda_n, L, amp in a notebook)
# ---------------------------------------------------------------------------

def parse_soundwave_par(path):
    """
    Read soundwave-related keys from a FARGO parameter file.
    Required: GAMMA, CPDG, TS_SOUNDWAVE, TH_SOUNDWAVE, EPS_SOUNDWAVE.
    Optional: REAL_EIGEN, IMAG_EIGEN, AMPLITUDE, RHOG, TGAS, DENSGR, DENSGI, ZMIN, ZMAX, OUTPUTDIR, …
    """
    key_to_field = {
        "GAMMA": ("gamma", float),
        "CPDG": ("cpdg", float),
        "TS_SOUNDWAVE": ("ts", float),
        "TH_SOUNDWAVE": ("th", float),
        "EPS_SOUNDWAVE": ("eps", float),
        "OUTPUTDIR": ("output_dir_rel", str),
        "SOUNDWAVE_USE_PAR_EIGEN": ("soundwave_use_par_eigen", int),
        "REAL_EIGEN": ("real_eigen", float),
        "IMAG_EIGEN": ("imag_eigen", float),
        "AMPLITUDE": ("amplitude", float),
        "RHOG": ("rhog", float),
        "TGAS": ("tg", float),
        "DENSGR": ("densgr", float),
        "DENSGI": ("densgi", float),
        "ZMIN": ("zmin", float),
        "ZMAX": ("zmax", float),
        "TDUST": ("tdust", float),
        "DT": ("dt", float),
        "NTOT": ("ntot", int),
    }
    raw = {}
    with open(path) as f:
        for line in f:
            line = line.split("#")[0].strip()
            if not line:
                continue
            parts = line.replace("\t", " ").split()
            if len(parts) < 2:
                continue
            name = parts[0].upper()
            if name not in key_to_field:
                continue
            field, typ = key_to_field[name]
            if typ is float:
                raw[field] = float(parts[1])
            elif typ is int:
                raw[field] = int(float(parts[1]))
            else:
                val = " ".join(parts[1:]).strip()
                if val.endswith("/"):
                    val = val[:-1]
                raw[field] = val
    need = ("gamma", "cpdg", "ts", "th", "eps")
    missing = [k for k in need if k not in raw]
    if missing:
        raise ValueError(
            "parse_soundwave_par(%s): missing %s (need GAMMA CPDG TS_SOUNDWAVE TH_SOUNDWAVE EPS_SOUNDWAVE)"
            % (path, missing)
        )
    if "output_dir_rel" not in raw:
        raw["output_dir_rel"] = "outputs/soundwave"
    return raw


def lambda_from_dispersion_relation(gamma, cpdg, ts, th, eps):
    """λ from the Python dispersion model (same as condinit when not using par eigen)."""
    from eigenvalue import get_lambda_for_condinit
    return get_lambda_for_condinit(
        gamma=gamma, ts=ts, th=th, cdg=cpdg, eps1=eps
    )


def analysis_kwargs_from_par(par_path):
    """
    Everything ``run_analysis`` needs, parsed from the same ``.par`` you passed to fargoboth.
    """
    p = parse_soundwave_par(par_path)
    use_eigen = int(p.get("soundwave_use_par_eigen", 1))
    if use_eigen != 0 and "real_eigen" in p and "imag_eigen" in p:
        lam = complex(p["real_eigen"], p["imag_eigen"])
    else:
        lam = lambda_from_dispersion_relation(
            p["gamma"], p["cpdg"], p["ts"], p["th"], p["eps"]
        )
    zmin = float(p.get("zmin", 0.0))
    zmax = float(p.get("zmax", 1.0))
    return {
        "gamma": p["gamma"],
        "cpdg": p["cpdg"],
        "ts": p["ts"],
        "th": p["th"],
        "eps": p["eps"],
        "lambda_complex": lam,
        "amplitude": float(p.get("amplitude", 1e-4)),
        "rhog": float(p.get("rhog", 1.0)),
        "Tg": float(p.get("tg", 1.0)),
        "densgr": float(p.get("densgr", 1.0)),
        "densgi": float(p.get("densgi", 0.0)),
        "L": zmax - zmin,
        "tdust": float(p.get("tdust", 1.0)),
        "dt": float(p.get("dt", 0.01)),
        **({"ntot": int(p["ntot"])} if "ntot" in p else {}),
    }


def run_analysis_from_par(par_path, project_root=".", save_figures=True, figure_dir=None):
    """
    Convenience: ``OutputDir`` in the par is taken relative to ``project_root``, then
    ``run_analysis`` is called with kwargs from the same file.
    """
    par_path = os.path.abspath(par_path)
    kw = analysis_kwargs_from_par(par_path)
    p = parse_soundwave_par(par_path)
    output_dir = os.path.join(os.path.abspath(project_root), p["output_dir_rel"])
    keys = (
        "gamma", "cpdg", "ts", "th", "eps", "lambda_complex",
        "amplitude", "rhog", "Tg", "densgr", "densgi", "L", "tdust",
        "dt",
    )
    call_kw = {k: kw[k] for k in keys}
    if "ntot" in kw:
        call_kw["ntot"] = kw["ntot"]
    return run_analysis(
        output_dir,
        save_figures=save_figures,
        figure_dir=figure_dir,
        **call_kw,
    )


# ---------------------------------------------------------------------------
# NLTEPPD-style eigen (real part of exp(i k x - i omega t), lambda = -i omega)
# ---------------------------------------------------------------------------

def eigen(fr, fi, kmode, x, lambda_n, time):
    lr = lambda_n.real
    li = lambda_n.imag
    return np.exp(lr * time) * (
        fr * np.cos(kmode * x + li * time) - fi * np.sin(kmode * x + li * time)
    )


def _time_scatter_indices(nsteps, max_points=50):
    nsteps = int(nsteps)
    if nsteps <= 0:
        return np.array([], dtype=int)
    if nsteps <= max_points:
        return np.arange(nsteps, dtype=int)
    return np.unique(np.round(np.linspace(0, nsteps - 1, max_points)).astype(int))


def _real_dtype_from_variables_par(output_dir):
    path = os.path.join(output_dir, "variables.par")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.split("#")[0].strip()
                if not line:
                    continue
                parts = line.replace("\t", " ").split()
                if not parts or parts[0].upper() != "REALTYPE":
                    continue
                if len(parts) < 2:
                    continue
                val = parts[1].lower()
                if "float32" in val or val == "float":
                    return np.float32
                return np.float64
    except OSError:
        pass
    return None


def _load_fargo_binary_field(path, n_cells_hint, real_dtype_hint):
    nbytes = os.path.getsize(path)

    def _read(dt, count):
        return np.fromfile(path, dtype=dt, count=count).astype(np.float64, copy=False)

    if n_cells_hint is not None:
        if n_cells_hint * 4 == nbytes:
            return _read(np.float32, n_cells_hint)
        if n_cells_hint * 8 == nbytes:
            return _read(np.float64, n_cells_hint)
    if real_dtype_hint is not None:
        w = int(np.dtype(real_dtype_hint).itemsize)
        if w > 0 and nbytes % w == 0:
            return _read(real_dtype_hint, nbytes // w)
    if nbytes % 8 == 0:
        return _read(np.float64, nbytes // 8)
    if nbytes % 4 == 0:
        return _read(np.float32, nbytes // 4)
    return np.fromfile(path, dtype=np.float64).astype(np.float64, copy=False)


def _linear_gas_T_amplitude(lambda_n, gamma, cs, kmode, rhog, delta_rhog, eps, ts):
    sum_term = eps / (1.0 + lambda_n * ts)
    return -(1.0 - (gamma * (-lambda_n**2) / (cs**2 * kmode**2) * (1.0 + sum_term))) * (
        delta_rhog / rhog
    )


def _count_dust_species(output_dir):
    """Number of dust fluids with ``dust1vz0.dat``, ``dust2vz0.dat``, …"""
    n = 0
    while os.path.isfile(os.path.join(output_dir, "dust%dvz0.dat" % (n + 1))):
        n += 1
    return n


# ---------------------------------------------------------------------------
# Main: load series like a notebook, compare to analytic, optional scatter PDF
# ---------------------------------------------------------------------------

def run_analysis(
    output_dir,
    gamma=1.4,
    cpdg=1.0,
    ts=0.01,
    th=0.1,
    eps=0.5,
    lambda_complex=None,
    amplitude=1e-4,
    rhog=1.0,
    Tg=1.0,
    densgr=1.0,
    densgi=0.0,
    L=1.0,
    tdust=1.0,
    ntot=None,
    dt=None,
    save_figures=True,
    figure_dir=None,
):
    """
    ``output_dir``: folder with ``domain_z.dat``, ``gasdens0.dat``, … (your ``OutputDir``).

    Pass Γ, CPDG, Ts, Th, ε, λ, etc. from ``analysis_kwargs_from_par(my.par)`` or by hand.

    If ``dust1vz0.dat`` exists, also loads ``dust*dens``, ``dust*vz``, ``dust*energy`` (NLTEPPD-style),
    compares to the same linear solution as ``condinit.c``, and plots gas + dust on one PDF
    (scatter + analytic curves). ``max_dev_*`` then include gas and all dust species.
    """
    if lambda_complex is None:
        lambda_complex = lambda_from_dispersion_relation(gamma, cpdg, ts, th, eps)

    kmode = 2 * np.pi / L
    cs = np.sqrt(gamma * Tg)
    cp_gas = gamma / (gamma - 1.0)
    cp_dust = cpdg * cp_gas
    delta_rhog = densgr + 1j * densgi
    vgas = 1j * lambda_complex / kmode * delta_rhog / rhog
    Tgas_factor = _linear_gas_T_amplitude(
        lambda_complex, gamma, cs, kmode, rhog, delta_rhog, eps, ts
    )
    delta_Tg = Tgas_factor * Tg
    # One dust species (same Ts, Th, ε as in par — matches condinit when SOUNDWAVE_USE_PAR_TS)
    vdust_c = vgas / (1.0 + lambda_complex * ts)
    rdust_c = eps * delta_rhog / (1.0 + lambda_complex * ts)
    tdust_c = delta_Tg / (1.0 + lambda_complex * th)

    # --- domain: x = np.loadtxt("domain_z.dat")[3:-3]  (notebook) ---
    domain_path = os.path.join(output_dir, "domain_z.dat")
    if not os.path.exists(domain_path):
        return {"error": "domain_z.dat not found in %s" % output_dir}
    x_raw = np.loadtxt(domain_path)
    if x_raw.ndim > 1:
        x_raw = x_raw[:, 3] if x_raw.shape[1] > 3 else x_raw[:, 0]
    x_int = x_raw[3:-3] if len(x_raw) > 6 else x_raw

    nx = 0
    if len(x_int) < 2:
        x_face = float(x_int[0])
        x_center = x_face
    else:
        if nx + 1 >= len(x_int):
            nx = max(0, len(x_int) - 2)
        x_face = float(x_int[nx])
        x_center = 0.5 * (float(x_int[nx]) + float(x_int[nx + 1]))

    n_frames = 0
    for i in range(2000):
        if os.path.exists(os.path.join(output_dir, "gasdens%d.dat" % i)):
            n_frames = i + 1
        else:
            break
    if n_frames == 0:
        return {"error": "No gasdens*.dat in %s" % output_dir}

    n_cells = len(x_int)
    real_dtype_hint = _real_dtype_from_variables_par(output_dir)
    dens0 = os.path.join(output_dir, "gasdens0.dat")
    if real_dtype_hint is None and os.path.isfile(dens0):
        nb = os.path.getsize(dens0)
        if n_cells * 4 == nb:
            real_dtype_hint = np.float32
        elif n_cells * 8 == nb:
            real_dtype_hint = np.float64

    if dt is None:
        dt = 0.01
    if ntot is not None:
        dt = (L / 64.0) * 0.3 / 10.0
    time = np.arange(n_frames) * dt

    j = min(nx, max(0, n_cells - 1))
    n_dust = _count_dust_species(output_dir)

    F3D_rho, F3D_v, F3D_e = [], [], []
    F3D_rd = [[] for _ in range(n_dust)]
    F3D_vd = [[] for _ in range(n_dust)]
    F3D_ed = [[] for _ in range(n_dust)]

    for i in range(n_frames):
        rho = _load_fargo_binary_field(
            os.path.join(output_dir, "gasdens%d.dat" % i), n_cells, real_dtype_hint
        )
        vz = _load_fargo_binary_field(
            os.path.join(output_dir, "gasvz%d.dat" % i), n_cells, real_dtype_hint
        )
        jj = min(j, len(rho) - 1, len(vz) - 1)
        F3D_rho.append(rho[jj] - rhog)
        F3D_v.append(vz[jj])
        epath = os.path.join(output_dir, "gasenergy%d.dat" % i)
        if os.path.exists(epath):
            e = _load_fargo_binary_field(epath, n_cells, real_dtype_hint)
            e0g = rhog * Tg * (cp_gas / gamma)
            F3D_e.append(e[min(jj, len(e) - 1)] - e0g)
        else:
            F3D_e.append(np.nan)

        for m in range(n_dust):
            did = m + 1
            rd = _load_fargo_binary_field(
                os.path.join(output_dir, "dust%ddens%d.dat" % (did, i)),
                n_cells,
                real_dtype_hint,
            )
            vd = _load_fargo_binary_field(
                os.path.join(output_dir, "dust%dvz%d.dat" % (did, i)),
                n_cells,
                real_dtype_hint,
            )
            jj2 = min(jj, len(rd) - 1, len(vd) - 1)
            F3D_rd[m].append(rd[jj2] - eps * rhog)
            F3D_vd[m].append(vd[jj2])
            edp = os.path.join(output_dir, "dust%denergy%d.dat" % (did, i))
            if os.path.exists(edp):
                en = _load_fargo_binary_field(edp, n_cells, real_dtype_hint)
                e0d = tdust * cp_dust * eps * rhog
                F3D_ed[m].append(en[min(jj2, len(en) - 1)] - e0d)
            else:
                F3D_ed[m].append(np.nan)

    F3D_rho = np.asarray(F3D_rho)
    F3D_v = np.asarray(F3D_v)
    F3D_e = np.asarray(F3D_e)
    F3D_rd = [np.asarray(a) for a in F3D_rd]
    F3D_vd = [np.asarray(a) for a in F3D_vd]
    F3D_ed = [np.asarray(a) for a in F3D_ed]

    # --- Analytic (same structure as condinit / NLTEPPD notebook) ---
    rho_anal = amplitude * eigen(densgr, densgi, kmode, x_center, lambda_complex, time)
    v_anal = amplitude * eigen(vgas.real, vgas.imag, kmode, x_face, lambda_complex, time)
    T_eigen = amplitude * eigen(delta_Tg.real, delta_Tg.imag, kmode, x_center, lambda_complex, time)
    e_anal = (cp_gas / gamma) * (rhog * T_eigen + Tg * rho_anal)

    rho_dust_anal = amplitude * eigen(
        rdust_c.real, rdust_c.imag, kmode, x_center, lambda_complex, time
    )
    v_dust_anal = amplitude * eigen(
        vdust_c.real, vdust_c.imag, kmode, x_face, lambda_complex, time
    )
    T_dust_eigen = amplitude * eigen(
        tdust_c.real, tdust_c.imag, kmode, x_center, lambda_complex, time
    )
    e_dust_anal = cp_dust * (
        eps * rhog * T_dust_eigen + Tg * rho_dust_anal
    )

    dev_rho = np.abs(F3D_rho - rho_anal)
    dev_v = np.abs(F3D_v - v_anal)
    dev_e = np.abs(F3D_e - e_anal) if np.any(np.isfinite(F3D_e)) else np.array([0.0])

    max_rho = [float(np.nanmax(dev_rho))]
    max_v = [float(np.nanmax(dev_v))]
    max_e = [float(np.nanmax(dev_e))] if np.any(np.isfinite(F3D_e)) else [0.0]
    for m in range(n_dust):
        max_rho.append(float(np.nanmax(np.abs(F3D_rd[m] - rho_dust_anal))))
        max_v.append(float(np.nanmax(np.abs(F3D_vd[m] - v_dust_anal))))
        if np.any(np.isfinite(F3D_ed[m])):
            max_e.append(float(np.nanmax(np.abs(F3D_ed[m] - e_dust_anal))))
    result = {
        "max_dev_rho": max(max_rho),
        "max_dev_v": max(max_v),
        "max_dev_e": max(max_e) if max_e else 0.0,
        "amplitude": amplitude,
        "lambda_used": lambda_complex,
        "n_dust_species": n_dust,
        "figure": None,
    }

    # --- PDF ---
    fig_dir = figure_dir or output_dir
    if save_figures and HAS_PLT:
        os.makedirs(fig_dir, exist_ok=True)
        idx = _time_scatter_indices(n_frames)
        t_sc = time[idx]
        amp = amplitude if amplitude else 1.0

        fig, axes = plt.subplots(3, 1, figsize=(7, 8), sharex=True)

        def scatter_ok(ax, y, label, color, ec):
            yv = np.asarray(y) / amp
            m = np.isfinite(t_sc) & np.isfinite(yv[idx])
            if np.any(m):
                ax.scatter(
                    t_sc[m],
                    yv[idx][m],
                    s=30,
                    c=color,
                    edgecolors=ec,
                    linewidths=0.5,
                    label=label,
                    zorder=3,
                )

        # Same panel order as typical NLTEPPD script: velocity, density, energy
        ax = axes[0]
        scatter_ok(ax, F3D_v, "Gas", "tab:orange", "darkred")
        for m in range(n_dust):
            scatter_ok(ax, F3D_vd[m], "Dust %d" % (m + 1), "cornflowerblue", "navy")
        if n_dust > 0:
            ax.plot(time, v_anal / amp, "k-", lw=1.5, label="Gas (linear)", zorder=1)
            ax.plot(time, v_dust_anal / amp, color="gray", lw=1.5, label="Dust (linear)", zorder=1)
        ax.set_ylabel(r"$v_z / A$")
        ax.set_title("Velocity")
        ax.legend(loc="best", fontsize=8, frameon=False)

        ax = axes[1]
        scatter_ok(ax, F3D_rho, "Gas", "tab:orange", "darkred")
        for m in range(n_dust):
            scatter_ok(ax, F3D_rd[m], "Dust %d" % (m + 1), "cornflowerblue", "navy")
        if n_dust > 0:
            ax.plot(time, rho_anal / amp, "k-", lw=1.5, label="Gas (linear)", zorder=1)
            ax.plot(time, rho_dust_anal / amp, color="gray", lw=1.5, label="Dust (linear)", zorder=1)
        ax.set_ylabel(r"$\delta\rho / A$")
        ax.set_title("Density perturbation")
        ax.legend(loc="best", fontsize=8, frameon=False)

        ax = axes[2]
        scatter_ok(ax, F3D_e, "Gas", "tab:orange", "darkred")
        for m in range(n_dust):
            scatter_ok(ax, F3D_ed[m], "Dust %d" % (m + 1), "cornflowerblue", "navy")
        if n_dust > 0:
            ax.plot(time, e_anal / amp, "k-", lw=1.5, label="Gas (linear)", zorder=1)
            ax.plot(time, e_dust_anal / amp, color="gray", lw=1.5, label="Dust (linear)", zorder=1)
        ax.set_ylabel(r"$\delta e / A$")
        ax.set_title("Energy perturbation")
        ax.set_xlabel("Time")
        ax.legend(loc="best", fontsize=8, frameon=False)

        fig.suptitle(
            "$\\gamma$=%.3g  CPDG=%.3g  $T_s$=%.3g  $T_h$=%.3g  $\\varepsilon$=%.3g\n"
            "$\\lambda$ = %.5g %+.5g i  ($n_{\\mathrm{dust}}$=%d)"
            % (
                gamma,
                cpdg,
                ts,
                th,
                eps,
                lambda_complex.real,
                lambda_complex.imag,
                n_dust,
            ),
            fontsize=10,
        )
        fig.tight_layout()
        fig_path = os.path.join(fig_dir, "soundwave_comparison.pdf")
        plt.savefig(fig_path, bbox_inches="tight")
        plt.close()
        result["figure"] = fig_path

    return result


def resolve_par_path(par_flag=None, par_positional=None, env_var="SOUNDWAVE_PAR"):
    """
    Pick the .par file path from, in order: ``par_flag`` (``--par``), positional path,
    then environment ``SOUNDWAVE_PAR``. Returns ``None`` if none set.
    """
    for candidate in (par_flag, par_positional):
        if candidate is not None and str(candidate).strip():
            return str(candidate).strip()
    env = os.environ.get(env_var, "").strip()
    return env or None


def main_cli(argv=None):
    """CLI: ``python analysis.py --par FILE.par [--project-root DIR]``."""
    import argparse

    argv = sys.argv[1:] if argv is None else argv
    ap = argparse.ArgumentParser(
        description="Load OutputDir from the par and write soundwave_comparison.pdf (+ metrics)."
    )
    ap.add_argument(
        "-p",
        "--par",
        metavar="FILE",
        help="Soundwave .par (same as for fargoboth). If omitted, use positional or SOUNDWAVE_PAR.",
    )
    ap.add_argument(
        "par_positional",
        nargs="?",
        metavar="PAR",
        help="Optional par path when -p/--par is not used",
    )
    ap.add_argument(
        "--project-root",
        default=".",
        help="Directory OutputDir in the par is relative to (default: .)",
    )
    ap.add_argument(
        "--no-figure",
        action="store_true",
        help="Do not write soundwave_comparison.pdf",
    )
    args = ap.parse_args(argv)
    par = resolve_par_path(args.par, args.par_positional)
    if not par:
        ap.error(
            "Give the par file: --par PATH, or PATH as positional argument, or set SOUNDWAVE_PAR"
        )
    res = run_analysis_from_par(
        par,
        project_root=args.project_root,
        save_figures=not args.no_figure,
    )
    if "error" in res:
        print(res["error"], file=sys.stderr)
        return 1
    print(
        "max_dev_rho=%.6e  max_dev_v=%.6e  max_dev_e=%.6e"
        % (res["max_dev_rho"], res["max_dev_v"], res["max_dev_e"])
    )
    if res.get("figure"):
        print("figure: %s" % res["figure"])
    return 0


if __name__ == "__main__":
    sys.exit(main_cli())
