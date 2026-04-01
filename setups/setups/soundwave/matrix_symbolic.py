"""
Symbolic 6x6 soundwave matrix M in SymPy (NLDTPPD.nb structure).
Parameters are symbols — nothing is substituted numerically.

  M = [ [0, 0, -I*ws, 0, 0, 0],
        [0, 0, 0, -I*eps1*ws, 0, 0],
        [-I*ws/gamma, 0, -eps1/ts, eps1/ts, -I*ws/gamma, 0],
        [0, 0, 1/ts, -1/ts, 0, 0],
        [0, 0, -I*(gamma-1)*ws, 0, -gamma*cdg*eps1/th, gamma*cdg*eps1/th],
        [0, 0, 0, 0, 1/th, -1/th] ]

Run:  python matrix_symbolic.py
Requires: sympy
"""
try:
    import sympy as sp
    from sympy import I, pprint
    from sympy.printing.latex import latex
except ImportError:
    raise SystemExit("Install sympy: pip install sympy")


def build_M_symbolic():
    """Return (M, symbols_dict) with all parameters as SymPy symbols."""
    gamma, ws, ts, th, cdg, eps1 = sp.symbols(
        "gamma w_s t_s t_h c_dg epsilon_1", real=True
    )
    # Matrix rows/cols 0..5
    M = sp.Matrix(
        [
            [0, 0, -I * ws, 0, 0, 0],
            [0, 0, 0, -I * eps1 * ws, 0, 0],
            [
                -I * ws / gamma,
                0,
                -eps1 / ts,
                eps1 / ts,
                -I * ws / gamma,
                0,
            ],
            [0, 0, 1 / ts, -1 / ts, 0, 0],
            [
                0,
                0,
                -I * (gamma - 1) * ws,
                0,
                -gamma * cdg * eps1 / th,
                gamma * cdg * eps1 / th,
            ],
            [0, 0, 0, 0, 1 / th, -1 / th],
        ]
    )
    syms = {
        "gamma": gamma,
        "w_s": ws,
        "t_s": ts,
        "t_h": th,
        "c_dg": cdg,
        "epsilon_1": eps1,
    }
    return M, syms


def main():
    M, syms = build_M_symbolic()
    print("Symbols (real):")
    for name, s in syms.items():
        print("  %s" % s)
    print()
    print("M (6x6), symbolic:")
    print()
    pprint(M, use_unicode=True)
    print()
    print("LaTeX (paste into notebook / paper):")
    print(latex(M))
    print()
    # Optional: eigenvalues are generally huge as RootOf; show charpoly degree note
    lam = sp.symbols("lambda", complex=True)
    I6 = sp.eye(6)
    char = (M - lam * I6).det()
    print("Characteristic polynomial det(M - lambda*I) is degree %s in lambda."
          % sp.degree(char, lam))
    print("(Expand with char.expand() if needed; roots are symbolic.)")


if __name__ == "__main__":
    main()
