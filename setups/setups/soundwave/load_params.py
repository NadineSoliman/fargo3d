"""
Load soundwave parameters from a ``.par`` file (for scripts / notebooks).

Example::

  from load_params import get_soundwave_params
  p = get_soundwave_params()  # default: setups/soundwave/soundwave.par
  gamma, cpdg, ts, th, eps = p["gamma"], p["cpdg"], p["ts"], p["th"], p["eps"]
"""
import os


def get_soundwave_params(par_path=None):
    """Parse ``soundwave.par`` (or ``par_path``) for Γ, CPDG, Ts, Th, ε."""
    if par_path is None:
        par_path = os.path.join(os.path.dirname(__file__), "soundwave.par")
    from analysis import parse_soundwave_par

    raw = parse_soundwave_par(par_path)
    return {
        "gamma": raw["gamma"],
        "cpdg": raw["cpdg"],
        "ts": raw["ts"],
        "th": raw["th"],
        "eps": raw["eps"],
    }
