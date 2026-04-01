"""
Eigenvalue computation for the linear sound wave (NLDTPPD).
Matches the 6x6 matrix from Mathematica NLDTPPD.nb:
  M = { {0, 0, -I ws, 0, 0, 0},
        {0, 0, 0, -I eps1 ws, 0, 0},
        {-I ws/gamma, 0, -eps1/ts, eps1/ts, -I ws/gamma, 0},
        {0, 0, 1/ts, -1/ts, 0, 0},
        {0, 0, -I (gamma-1) ws, 0, -gamma cdg eps1/th, gamma cdg eps1/th},
        {0, 0, 0, 0, 1/th, -1/th} };
With R=1, Tg=1, gamma, cs = sqrt(gamma*R*Tg), ws = cs*2*pi/L, L=1, ts, th, cdg, eps1.
Returns the oscillatory + decaying eigenvalue (Re < 0, Im != 0), choosing Im > 0 to match Mathematica Eigenvalues[].
"""
import numpy as np


def build_M(gamma=1.4, R=1.0, Tg=1.0, L=1.0, ts=0.01, th=0.1, cdg=1.0, eps1=0.5):
    cs = np.sqrt(gamma * R * Tg)
    ws = cs * 2 * np.pi / L
    j = 1j
    M = np.array([
        [0, 0, -j * ws, 0, 0, 0],
        [0, 0, 0, -j * eps1 * ws, 0, 0],
        [-j * ws / gamma, 0, -eps1 / ts, eps1 / ts, -j * ws / gamma, 0],
        [0, 0, 1 / ts, -1 / ts, 0, 0],
        [0, 0, -j * (gamma - 1) * ws, 0, -gamma * cdg * eps1 / th, gamma * cdg * eps1 / th],
        [0, 0, 0, 0, 1 / th, -1 / th],
    ], dtype=complex)
    return M


def compute_eigenvalues(gamma=1.4, R=1.0, Tg=1.0, L=1.0, ts=0.01, th=0.1, cdg=1.0, eps1=0.5):
    M = build_M(gamma=gamma, R=R, Tg=Tg, L=L, ts=ts, th=th, cdg=cdg, eps1=eps1)
    w, v = np.linalg.eig(M)
    return w


def select_oscillatory_decay(w, prefer_mathematica=True):
    """
    Select an eigenvalue that is oscillatory (Im != 0) and decaying (Re < 0).
    prefer_mathematica: if True (default), prefer Im > 0 to match Eigenvalues in
    Mathematica NLDTPPD.nb (complex conjugate pair gives ±Im; NB picks +Im).
    If False, prefer Im < 0.
    """
    candidates = []
    for lam in w:
        re, im = lam.real, lam.imag
        if re >= 0:
            continue  # skip growing or marginal
        if np.abs(im) < 1e-10:
            continue  # skip purely decaying
        candidates.append(lam)
    if not candidates:
        return None
    if prefer_mathematica:
        with_pos_im = [c for c in candidates if c.imag > 0]
        if with_pos_im:
            # Largest positive Im first (principal oscillatory mode)
            with_pos_im.sort(key=lambda z: (-z.imag, z.real))
            return with_pos_im[0]
    else:
        with_neg_im = [c for c in candidates if c.imag < 0]
        if with_neg_im:
            with_neg_im.sort(key=lambda z: (z.imag, z.real))
            return with_neg_im[0]
    return candidates[0]


def get_lambda_for_condinit(gamma=1.4, R=1.0, Tg=1.0, L=1.0, ts=0.01, th=0.1, cdg=1.0, eps1=0.5):
    w = compute_eigenvalues(gamma=gamma, R=R, Tg=Tg, L=L, ts=ts, th=th, cdg=cdg, eps1=eps1)
    lam = select_oscillatory_decay(w, prefer_mathematica=True)
    if lam is None:
        raise ValueError("No oscillatory+decaying eigenvalue found. Eigenvalues: %s" % w)
    return lam


def matrix_to_human_readable(M, col_width=24):
    """6x6 complex matrix: one row per line, each entry as (real + imag*I)."""
    n, m = M.shape
    lines = []
    lines.append("        " + "".join(str(j).center(col_width) for j in range(m)))
    for i in range(n):
        parts = ["[%d]    " % i]
        for j in range(m):
            z = M[i, j]
            if abs(z.imag) < 1e-14:
                s = "% .10g" % z.real
            elif abs(z.real) < 1e-14:
                s = "% .10g*I" % z.imag
            else:
                s = "(% .8g%+.8g*I)" % (z.real, z.imag)
            parts.append(s.ljust(col_width))
        lines.append("".join(parts))
    return "\n".join(lines)


def residual_eigenvalue(M, lam):
    """
    If lam is an exact eigenvalue, (M - lam*I) is singular.
    Return smallest singular value of (M - lam*I) (should be ~0 if lam is in spectrum).
    """
    n = M.shape[0]
    A = M - lam * np.eye(n, dtype=complex)
    s = np.linalg.svd(A, compute_uv=False)
    return float(s[-1])


def print_matrix_report(
    gamma=1.4,
    R=1.0,
    Tg=1.0,
    L=1.0,
    ts=0.01,
    th=0.1,
    cdg=1.0,
    eps1=0.5,
    paper_lambda=None,
):
    """
    Print cs, ws, full M, all eigenvalues, and (if given) SVD residual for a candidate lambda.
    paper_lambda: e.g. complex(-0.17149106273322487, 5.73870605564313330) for NLDTPPD.nb with cdg=1.
    """
    cs = np.sqrt(gamma * R * Tg)
    ws = cs * 2 * np.pi / L
    M = build_M(gamma=gamma, R=R, Tg=Tg, L=L, ts=ts, th=th, cdg=cdg, eps1=eps1)
    w = compute_eigenvalues(gamma=gamma, R=R, Tg=Tg, L=L, ts=ts, th=th, cdg=cdg, eps1=eps1)

    print("=" * 72)
    print("Soundwave 6x6 matrix M (same structure as Mathematica NLDTPPD.nb)")
    print("=" * 72)
    print("Parameters:  gamma=%.12g  R=%.12g  Tg=%.12g  L=%.12g" % (gamma, R, Tg, L))
    print("             ts=%.12g  th=%.12g  cdg=%.12g  eps1=%.12g" % (ts, th, cdg, eps1))
    print("Derived:     cs = sqrt(gamma*R*Tg) = %.12g" % cs)
    print("             ws = cs*2*pi/L       = %.12g" % ws)
    print()
    print("M (rows 0..5, cols 0..5); each entry shown as (real + imag*I) or real-only:")
    print(matrix_to_human_readable(M))
    print()
    print("All eigenvalues of M (numpy.linalg.eig):")
    for k, z in enumerate(sorted(w, key=lambda z: (z.real, z.imag))):
        print("  [%d]  %.18g %+.18g*I" % (k, z.real, z.imag))
    print()
    lam_sel = select_oscillatory_decay(w, prefer_mathematica=True)
    if lam_sel is not None:
        r = residual_eigenvalue(M, lam_sel)
        print("Selected oscillatory+decaying (Im>0, Re<0):  %.18g %+.18g*I" % (lam_sel.real, lam_sel.imag))
        print("  min singular value of (M - lambda*I): %.3e  (should be ~0 if exact)" % r)
    if paper_lambda is not None:
        r_p = residual_eigenvalue(M, paper_lambda)
        print()
        print("Paper / NB candidate lambda:")
        print("  %.18g %+.18g*I" % (paper_lambda.real, paper_lambda.imag))
        print("  min singular value of (M - lambda*I): %.3e" % r_p)
        if r_p > 1e-6:
            print("  NOTE: large residual usually means this lambda is NOT an eigenvalue for THIS M.")
            print("        The published value is for the NB parameters (often cdg=1, not CPDG=0.01).")
    print("=" * 72)


if __name__ == "__main__":
    paper = -0.17149106273322487 + 5.73870605564313330j

    print("\n>>> Case A: Mathematica / paper parameters (cdg = 1.0)\n")
    print_matrix_report(gamma=1.4, ts=0.01, th=0.1, cdg=1.0, eps1=0.5, paper_lambda=paper)

    print("\n>>> Case B: Typical par file CPDG = 0.01 (different M -> different eigenvalues)\n")
    print_matrix_report(gamma=1.4, ts=0.01, th=0.1, cdg=0.01, eps1=0.5, paper_lambda=paper)
