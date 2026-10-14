"""FKS2 Table 7 upper ends L* (E_pi(x) <= B(x) for 1 <= log x <= L*) under corrected Buthe inputs.
eps_pi,num from FKS2 Table 4 (log-linear interpolation, fks2_table7.T4).  FKS2 use Buthe 2016 Thm 1 (printed
E2) for 10^19 <= x < e^2073 and FKS1 Thm 1.1 beyond, so the correction is applied only for L < 2073.
Correction = additive increase of eps_psi at L (evaluated at the crossing): fixed parameters (the
printed-E2 optimum) or re-optimised, for E2 with log(3c) and for Bhattacharjee's corrected E2.
Uniform +2.6% / +8% perturbations are applied only to L < 2073 as well."""
import math, sys
sys.path.insert(0, '.')
sys.stdout.reconfigure(encoding="utf-8")
from fks2_table7 import eps_pi, ROWS, crossing
from e2_variants import delta, optimise
H = 3.0e12
CUT = 2073.0
def deltas(L):
    Lr = round(L, 1)
    r0 = optimise(Lr, 'printed', H, c0=33.3, a0=max(0.002, min(0.12, 6.0 / Lr)))
    d0, c, T, a, _ = r0
    out = {}
    for v in ('log3c', 'thesis'):
        out['fixed_' + v] = delta(Lr, c, T, a, v)[0] - d0
        out['reopt_' + v] = optimise(Lr, v, H, c0=c + 1.5, T0=T, a0=a)[0] - d0
    return out, d0
if __name__ == "__main__":
    cols = ['fixed_log3c', 'reopt_log3c', 'fixed_thesis', 'reopt_thesis']
    print(f"{'B(x)':24s} {'stated':>8s} {'repro':>9s} | " + ' '.join(f'{c:>13s}' for c in cols) + f" | {'+2.6%':>9s} {'+8%':>9s}")
    for name, B, stated in ROWS:
        Ls = crossing(B, eps_pi)
        if Ls < CUT:
            ds, d0 = deltas(Ls)
        else:
            ds, d0 = {k: 0.0 for k in cols}, float('nan')
        res = [crossing(B, lambda l, d=ds[k]: eps_pi(l) + (d if l < CUT else 0.0)) for k in cols]
        L26 = crossing(B, lambda l: eps_pi(l) * (1.026 if l < CUT else 1.0))
        L8 = crossing(B, lambda l: eps_pi(l) * (1.08 if l < CUT else 1.0))
        print(f"{name:24s} {stated:8.2f} {Ls:9.3f} | " + ' '.join(f'{x:13.3f}' for x in res) + f" | {L26:9.3f} {L8:9.3f}"
              + (f"   [rel. increase at L*: " + ', '.join(f'{k}={ds[k]/eps_pi(Ls):+.2%}' for k in cols) + "]" if Ls < CUT else "   [L* beyond 2073: FKS1 range, no correction]"), flush=True)
