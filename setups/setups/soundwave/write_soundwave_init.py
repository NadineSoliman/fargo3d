#!/usr/bin/env python
"""
Write soundwave_init.dat and optionally update .par for a given parameter set.
Usage:
  python write_soundwave_init.py [--gamma 1.4] [--cpdg 0.01] [--ts 0.1] [--th 0.1] [--eps 2.24] [--outdir .]
Reads/writes from the setup directory or --outdir (for test runs).
"""
import argparse
import os

from eigenvalue import get_lambda_for_condinit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", type=float, default=1.4)
    ap.add_argument("--cpdg", type=float, default=0.01)
    ap.add_argument("--ts", type=float, default=0.1)
    ap.add_argument("--th", type=float, default=0.1)
    ap.add_argument("--eps", type=float, default=2.24)
    ap.add_argument("--R", type=float, default=1.0)
    ap.add_argument("--Tg", type=float, default=1.0)
    ap.add_argument("--L", type=float, default=1.0)
    ap.add_argument("--outdir", default=None, help="Directory to write soundwave_init.dat (default: script dir)")
    ap.add_argument("--write_par", action="store_true", help="Also write/update Gamma and CPDG in soundwave.par")
    args = ap.parse_args()

    outdir = args.outdir or os.path.dirname(os.path.abspath(__file__))
    os.makedirs(outdir, exist_ok=True)

    lam = get_lambda_for_condinit(
        gamma=args.gamma, R=args.R, Tg=args.Tg, L=args.L,
        ts=args.ts, th=args.th, cdg=args.cpdg, eps1=args.eps,
    )
    init_path = os.path.join(outdir, "soundwave_init.dat")
    with open(init_path, "w") as f:
        # lambda_re lambda_im Ts Th eps (one dust)
        f.write("%.18g %.18g %.18g %.18g %.18g\n" % (
            lam.real, lam.imag, args.ts, args.th, args.eps,
        ))
    print("Wrote %s (lambda = %g %+g*I)" % (init_path, lam.real, lam.imag))

    if args.write_par:
        par_path = os.path.join(outdir, "soundwave.par")
        if os.path.exists(par_path):
            with open(par_path) as f:
                lines = f.readlines()
            new_lines = []
            for line in lines:
                if line.strip().startswith("Gamma") or line.strip().startswith("Gamma\t"):
                    new_lines.append("Gamma\t\t%g\n" % args.gamma)
                elif "CPDG" in line.split()[0] if line.split() else False:
                    new_lines.append("CPDG\t\t%g\n" % args.cpdg)
                else:
                    new_lines.append(line)
            with open(par_path, "w") as f:
                f.writelines(new_lines)
            print("Updated Gamma and CPDG in %s" % par_path)
        else:
            print("No existing soundwave.par in %s, skipping --write_par" % outdir)


if __name__ == "__main__":
    main()
