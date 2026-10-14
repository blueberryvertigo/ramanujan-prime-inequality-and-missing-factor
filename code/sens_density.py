"""
SENSITIVITY STUDY (hypothetical; not a theorem): how far would Theorem C's threshold move if the explicit
zero-density bound N(sigma, T) of FKS/KLN (second block of Table 7, a1 = 0.618) were smaller by a factor K?
Everything else as in fks_fine_R.py (MTY R = 5.558691, H0 = 3e12).  Both uses of the density bound are scaled:
fks_eps.Nt (eps4) and fks_eps.B0 (eps3).  Rows every STEP in log x0, each valid on [L, L + STEP] (Lmax).
argv: K  LO  [STEP]
"""
import sys, time
from mpmath import mpf
import fks_eps as F
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
K = float(sys.argv[1]); LO = int(sys.argv[2]); STEP = int(sys.argv[3]) if len(sys.argv) > 3 else 2
HI = 3040
F.DENSITY['mode'] = 'KLN618'
_Nt, _B0 = F.Nt, F.B0
F.Nt = lambda s, t: _Nt(s, t) / K
F.B0 = lambda s, U, V: _B0(s, U, V) / K
t0 = time.time()
rows = []
for L in range(LO, HI + 1, STEP):
    try:
        e, s2 = F.best_eps_table7(L, N=10000, K=100000, R=5.558691, Lmax=L + STEP)
    except (ValueError, AssertionError):
        continue
    rows.append((L, I(repr(float(e)))))
print(f"  K = {K}: {len(rows)} rows in {time.time() - t0:.0f} s", flush=True)
def f(l):
    if l.b > HI + STEP: return I(1)
    cand = [v for L, v in rows if L <= l.a and l.b <= L + STEP]
    return cand[-1] if cand else I(1)
B = CH.FnBounds([(f, float(LO))], kmax=10050)
last_fail = None
for a, b in V.grid(LO, HI + 0.5, 0.5, 'B'):
    if not CH.ok_adaptive(B, mpf(a), mpf(b)): last_fail = b
print(f"[sensitivity K={K}] FKS rows alone verify (C) for L in [{float(last_fail) if last_fail else LO}, {HI + 0.5}] "
      f"(grid 0.5; Theorem C threshold would move from 3035.00 to about this)   ({time.time() - t0:.0f} s)", flush=True)
