"""
Theorem D (hybrid): if RH holds up to height T, then for log x in [L_A - 1, L_top]
  |psi(x) - x| <= x * min( [CH, Prop 9.1](T),  FKS Thm 1.1 with verification height H0 = T  (fks_eps),
                           [CH, Cor 1.3], FKS Table 5 / Cor 1.4 ),
and (C) is checked on [L_A - 1, L_top]; Theorem C covers L >= L_top (argv[3]).  Bisection in T.
argv: R  tag  L_top
"""
import math, sys, time
from mpmath import mpf
import fks_eps as F
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
R = float(sys.argv[1]); tag = sys.argv[2]; L_TOP = float(sys.argv[3])
if len(sys.argv) > 4: F.DENSITY['mode'] = sys.argv[4]   # KLN | KLN618 | BEL618
LA = 753.985
ROW_LO, STEP = 1700, 4

def fks_rows_fn(T):
    rows = []
    for L in range(ROW_LO, int(L_TOP) + 2, STEP):
        try:
            e, s2 = F.best_eps_table7(L, N=4000, K=40000, R=R, H0v=T, Lmax=L + STEP)
        except ValueError:
            continue
        rows.append((L, I(repr(float(e)))))
    last = rows[-1][0] if rows else ROW_LO
    def f(l):
        if l.b > last + STEP: return I(1)
        cand = [v for L, v in rows if L <= l.a]
        return cand[-1] if cand else I(1)
    return f, float(ROW_LO)

def closes(T):
    g, gmin = CH.ch_prop91(T)
    f, fmin = fks_rows_fn(T)
    B = CH.FnBounds([(CH.ch_cor13, 0.0), (g, gmin), (f, fmin)])
    e, fl = CH.run(B, LA - 1, L_TOP + 0.5, 0.5)
    return fl is None, fl

if __name__ == "__main__":
    t0 = time.time()
    lo_T, hi_T = 1.5e14, 6e14
    assert closes(hi_T)[0]
    for _ in range(11):
        mid = (lo_T * hi_T) ** 0.5
        ok, fl = closes(mid)
        print(f"   [D-hybrid, R={R}] T = {mid:.6g}: {'closes' if ok else 'fails at ' + str(fl)}   ({time.time()-t0:.0f} s)", flush=True)
        if ok: hi_T = mid
        else: lo_T = mid
    print(f"[D-hybrid, R={R}, {tag}] RH to height {hi_T:.6g} suffices; fails at {lo_T:.6g}   ({time.time() - t0:.0f} s)")
