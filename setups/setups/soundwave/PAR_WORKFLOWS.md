# Soundwave: parameter file workflows

## Eigenvalue in the `.par` file (no `soundwave_init.dat`)

The setup reads the linear eigenvalue from the parameter file when **`SOUNDWAVE_USE_PAR_EIGEN`** is **1**:

| Variable | Meaning |
|----------|---------|
| **`REAL_EIGEN`** | Real part of λ |
| **`IMAG_EIGEN`** | Imaginary part of λ |
| **`SOUNDWAVE_USE_PAR_EIGEN`** | `1` = use REAL/IMAG from par; `0` = legacy: read `soundwave_init.dat` if present |

With **`SOUNDWAVE_USE_PAR_EIGEN 1`**, you can switch between runs by **only changing `.par` files** (and running `./fargoboth yourcase.par`). **No recompile** between cases. Rebuild **`make SETUP=soundwave`** only after changing **`condinit.c`** or other source.

Other soundwave controls include **`SOUNDWAVE_USE_PAR_TS`**, **`TS_SOUNDWAVE`**, **`TH_SOUNDWAVE`**, **`EPS_SOUNDWAVE`**, **`TDUST`**, **`Amplitude`**, background fields (**`RHOG`**, **`TGAS`**, **`DENSGR`**, **`DENSGI`**), and **`ZMIN`** / **`ZMAX`** (domain). See the committed **`soundwave.par`** template.

---

## Several hand-maintained `.par` files

1. Copy **`soundwave.par`** to e.g. **`setups/soundwave/case_a.par`** / **`case_b.par`**.
2. Set **`OutputDir`**, physics, and eigen lines as needed.
3. Run:

```bash
./fargoboth setups/soundwave/case_a.par
```

---

## Post-run plots (Python)

**Option A — one combined PDF** (gas + dust, same logic as integrated `analysis.run_analysis`):

```bash
cd /path/to/fargoboth
python setups/soundwave/analysis.py --par setups/soundwave/soundwave.par --project-root .
# short: -p FILE.par ; or pass FILE.par as the only positional argument;
# or: export SOUNDWAVE_PAR=setups/soundwave/soundwave.par  then  python setups/soundwave/analysis.py --project-root .
```

Or from Python:

```python
import sys
sys.path.insert(0, "setups/soundwave")
from analysis import run_analysis_from_par

run_analysis_from_par("setups/soundwave/soundwave.par", project_root=".")
```

**Option B — three NLTEPPD-style PDFs** (velocity / density / energy, LaTeX-style rcParams if TeX is installed):

```bash
python setups/soundwave/plot_nlteppd_from_par.py --par setups/soundwave/soundwave.par --project-root .
# or: python setups/soundwave/plot_nlteppd_from_par.py setups/soundwave/soundwave.par
```

Physics (**Γ, CPDG, Ts, Th, ε, λ, Amplitude, RHOG, TGAS, DENSGR, DENSGI, Zmin/Zmax→L, TDUST, OutputDir**) and time step (**DT**, optional **Ntot**) are read from the **same `.par`** file.

**Minimal in a notebook** — only pull numbers from the par:

```python
import sys
sys.path.insert(0, "setups/soundwave")
from analysis import analysis_kwargs_from_par, parse_soundwave_par

par = "setups/soundwave/soundwave.par"
kw = analysis_kwargs_from_par(par)
raw = parse_soundwave_par(par)
# kw["lambda_complex"], kw["gamma"], kw["ts"], …
# data directory (relative to repo root):
out = raw["output_dir_rel"]
```

(`project_root` is the repo root; **`OutputDir`** in the par is relative to that.)

---

## Legacy: `soundwave_init.dat`

Set **`SOUNDWAVE_USE_PAR_EIGEN`** **`0`** in the `.par` file. Then condinit reads **`soundwave_init.dat`** (five numbers: Re λ, Im λ, Ts, Th, ε) if present. Use **`write_soundwave_init.py`** or edit by hand.

---

## `condinit.c`

Edit when you need different **default** physics or init logic. Normal case switching uses **only** different `.par` files + **`REAL_EIGEN` / `IMAG_EIGEN`** (when par-eigen mode is on).
