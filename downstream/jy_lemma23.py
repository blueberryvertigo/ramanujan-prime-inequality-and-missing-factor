"""Johnston-Yang (arXiv:2204.01980v2; JMAA 2023), Lemma 2.3: |psi(x)-x| <= 1.570e-12 x for x >= e^2000,
from Buthe 2016 Thm 1 with H = 3,000,175,332,800, alpha = 0, c = 33.6, eps = 1.12e-11 (printed E2, no log 3c).
Section 3.1 uses it on e^2000 < x <= e^2488 to get Theorem 1.1:
    |psi(x)-x| <= f(log x) x,   f(L) = 9.39 L^1.515 exp(-0.8274 sqrt L),
and Section 3.2 starts at x0 = e^2488.  Compute the lemma value as printed, with log(3c) at their
parameters, re-optimised with log(3c) (alpha = 0 and alpha free), and the largest L with f(L) >= value."""
import sys, math
sys.path.insert(0, '.')
from zeros_model import delta, optimise
from mpmath import mp, mpf, exp, sqrt, findroot
mp.dps = 30
H = 3000175332800.0
f = lambda L: mpf('9.39') * mpf(L)**mpf('1.515') * exp(-mpf('0.8274') * sqrt(mpf(L)))
cross = lambda e: findroot(lambda L: f(L) - mpf(e), 2485)
print("f(2000) =", mp.nstr(f(2000), 6), "  f(2488) =", mp.nstr(f(2488), 8))
T = 33.6 / 1.12e-11
p0 = delta(2000, 33.6, T, 0.0, False)[0]; p1 = delta(2000, 33.6, T, 0.0, True)[0]
print(f"Lemma 2.3 parameters (c=33.6, T=c/eps={T:.6e}, alpha=0): no log3c {p0:.6e};  with log3c {p1:.6e} (+{p1/p0-1:.2%})")
best0 = min((delta(2000, 32 + 0.01 * i, H, 0.0, False)[0], 32 + 0.01 * i) for i in range(400))
best1 = min((delta(2000, 33 + 0.01 * i, H, 0.0, True)[0], 33 + 0.01 * i) for i in range(400))
print(f"alpha = 0, T = H: optimum no log3c {best0[0]:.6e} (c={best0[1]:.2f});  optimum with log3c {best1[0]:.6e} (c={best1[1]:.2f})")
r1 = optimise(2000, True, H, c0=35.0, T0=H, a0=0.003)
print(f"alpha free, with log3c: {r1[0]:.6e} (c={r1[1]:.3f}, alpha={r1[3]:.4g})")
for name, e in [("printed 1.570e-12", 1.570e-12), ("with log3c, c=33.6, alpha=0", p1),
                ("re-opt, alpha=0", best1[0]), ("re-opt, alpha free", r1[0]),
                ("1.570e-12 x 1.026", 1.570e-12 * 1.026), ("1.570e-12 x 1.08", 1.570e-12 * 1.08)]:
    print(f"{name:30s} eps={e:.5e}: f(L) >= eps exactly for L <= {mp.nstr(cross(e), 7)}  (gap below 2488: {mp.nstr(2488 - cross(e), 4)})")
# what A would be needed at L = 2488 with B, C unchanged
for name, e in [("printed", 1.570e-12), ("re-opt alpha free", r1[0])]:
    print(f"A needed at L=2488 for {name}: {mp.nstr(mpf(e) / (f(2488) / mpf('9.39')), 6)}")
