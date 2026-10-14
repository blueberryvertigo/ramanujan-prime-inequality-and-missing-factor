"""Largest theta such that a zero-free half-plane Re s > theta, with refereed inputs only, gives (C) on [L_lo, 3036]
(L_lo = 672.5 < L_A' = 672.849, so Proposition A' (thesis-corrected) and Theorem C close the rest).
eps_theta(L) = 200 L^2/T + 2 x^{-1/2} log^2(H/2pi)/(4pi) + 2 x^{theta-1} log^2(2.45T/2pi)/(4pi), log T = L/2 - 2.
Decreasing in L for L >= 110 when theta <= 0.99 (checked numerically below).  argv: L_lo"""
import sys, time
from mpmath import iv, mpf
import ch_ramanujan as CH, verify_ramanujan as V
from qrh_ramanujan import LOGH2, PI
I = V.I
L_LO = float(sys.argv[1]) if len(sys.argv) > 1 else 672.5
def make(theta):
    th = I(theta)
    def eps(L):
        L = I(L); lT = L / 2 - 2
        return (200 * L**2 * iv.exp(-lT) + 2 * iv.exp(-L / 2) * LOGH2 / (4 * PI)
                + 2 * iv.exp((th - 1) * L) * (lT + iv.log(I('2.45') / (2 * PI)))**2 / (4 * PI))
    lmin = max(110.0, 2 / (1 - theta) + 6)          # e^{(theta-1)L}(L/2-2.94)^2 decreases for L > 2/(1-theta) + 5.9
    vals = [eps(L).b for L in range(int(lmin), 3100, 7)]
    assert all(a >= b for a, b in zip(vals, vals[1:])), "not monotone"
    return (lambda l: I(eps(I(l.a)).b)), lmin
def closes(theta):
    f, lmin = make(theta)
    B = CH.FnBounds([(f, lmin)], kmax=3100)
    e, fl = CH.run(B, L_LO, 3036.0, 0.5)
    return fl is None, fl
t0 = time.time()
lo, hi = 0.875, 0.99
assert closes(lo)[0]
for _ in range(14):
    mid = (lo + hi) / 2
    ok, fl = closes(mid)
    print(f"   theta = {mid:.6f}: {'closes' if ok else 'fails at ' + str(fl)}   ({time.time()-t0:.0f} s)", flush=True)
    if ok: lo = mid
    else: hi = mid
print(f"[B-ref] Re s > theta zero-free with theta <= {lo:.6f} suffices (refereed inputs); fails at {hi:.6f}")
