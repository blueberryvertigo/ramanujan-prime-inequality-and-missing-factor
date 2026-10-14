"""
Theorem B' (classical-shape zero-free region): if zeta(beta + i gamma) != 0 whenever |gamma| > H0 = 3e12 and
beta > 1 - 1/(R log|gamma|), then FKS Theorem 1.1 (fks_eps.py) holds with that R.  Combined with
[CH, Cor 1.3] it gives psi-bounds on [L_A - 1, L_C]; bisection for the largest R for which (C) holds
on [L_A - 1, L_C], so that with Theorems A and C, N_R = 38,358,837,683.

Rows: eps(x0) at log x0 = L for L on a grid of step STEP, each valid for all x >= x0 (FKS Thm 1.1).
FKS state Theorem 1.1 for log x0 > 1000; the ingredients (their Props 3.4, 3.6, 3.8, 3.11, 3.14)
require only x > e^50, 3 log x < T < sqrt(x)/3 and T > t0, all checked in fks_eps.eps_parts, and the
monotonicity in x of Prop. 3.14 holds for x > exp(R e^2).  We use rows from log x0 = 700.
argv: L_C (top of the open range covered by Theorem C)
"""
import sys, time
from mpmath import mpf
import fks_eps as F
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
LA = 753.985
L_C = float(sys.argv[1])
if len(sys.argv) > 2: F.DENSITY['mode'] = sys.argv[2]   # KLN | KLN618 | BEL618
ROW_LO, STEP = 700, 3

SIGMA2 = [0.95, 0.96, 0.97, 0.98, 0.985, 0.990, 0.993, 0.996]

def best(L, R, Lmax):
    """least eps over sigma_2: Table 7 constants at the 0.99x points, Table 8 (interval) elsewhere;
    t0 is chosen so that the row is valid on [L, Lmax]."""
    out = []
    for s2 in SIGMA2:
        try:
            out.append(F.eps(L, s2, N=4000, K=40000, R=R, Lmax=Lmax))
        except AssertionError:
            pass
    return min(out) if out else None

def rows_fn(R):
    rows = []
    grid_L = list(range(ROW_LO, 1200, STEP)) + list(range(1200, int(L_C) + 2, 25))
    for i, L in enumerate(grid_L):
        Lnext = grid_L[i + 1] if i + 1 < len(grid_L) else L + 25
        e = best(L, R, Lnext)
        if e is not None:
            rows.append((L, Lnext, I(repr(float(e)))))
    def f(l):
        cand = [v for L, Ln, v in rows if L <= l.a and l.b <= Ln]
        return cand[-1] if cand else I(1)
    return f, float(ROW_LO)

def closes(R):
    f, fmin = rows_fn(R)
    B = CH.FnBounds([(CH.ch_cor13, 0.0), (f, fmin)])
    e, fl = CH.run(B, LA - 1, L_C + 0.5, 0.5)
    return fl is None, fl

if __name__ == "__main__":
    t0 = time.time()
    lo, hi = 1.0, 3.0             # lo closes, hi fails
    assert closes(lo)[0]
    for _ in range(12):
        mid = (lo + hi) / 2
        ok, fl = closes(mid)
        print(f"   [B'] R = {mid:.5f}: {'closes' if ok else 'fails at ' + str(fl)}   ({time.time()-t0:.0f} s)", flush=True)
        if ok: lo = mid
        else: hi = mid
    print(f"[B'] zero-free region beta <= 1 - 1/(R log|gamma|) for |gamma| > 3e12 with R <= {lo:.5f} suffices; "
          f"fails at {hi:.5f}   ({time.time() - t0:.0f} s)")
