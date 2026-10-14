"""Non-certified cross-check of verify_ramanujan.py: recompute the main term
L^6 [ j(L-1)/L - j(L)^2 ] directly from mpmath's li (no asymptotic expansion), and the
error budget with plain high-precision floats, and compare with the interval enclosures."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mpm
from mpmath import mp, iv
import verify_ramanujan as V

mp.dps = 120
def j_direct(l):
    x = mpm.exp(l)
    J = mpm.li(x) - mpm.li(2) + 2/mpm.log(2)          # J(x) = x/L + int_2^x dt/log^2 t
    return J/x

print(" L      direct L^6[j(L-1)/L-j(L)^2]   interval main+corr enclosure            inside?")
for L in [43, 50, 80, 103, 150, 300, 642, 1000, 3158]:
    Lm = mpm.mpf(L)
    direct = Lm**6*(j_direct(Lm-1)/Lm - j_direct(Lm)**2)
    LI = iv.mpf(L)
    mg, main, corr, err = V.margin(LI, V.d_U if L <= 700 else V.d_Q)
    enc = main + corr
    inside = enc.a <= direct <= enc.b
    print(f"{L:5d}  {mpm.nstr(direct, 18):>24}   [{mpm.nstr(mpm.mpf(enc.a), 18)}, {mpm.nstr(mpm.mpf(enc.b), 18)}]  {inside}")

# Simple-minded budget check: (R) needs roughly 3*eps < L^-4 (1 + 14/L) ; compare with the
# rigorous endpoint 642.83 found by the interval run.
eps = mpm.mpf('1.98895e-12')
f = lambda L: (1 + 14/L)/(3*L**4) - eps
print("heuristic crossover for eps(600,700):", mpm.nstr(mpm.findroot(f, 640), 8))
