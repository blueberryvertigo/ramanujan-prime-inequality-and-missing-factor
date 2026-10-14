"""Lichtman-Pomerance (arXiv:1705.02442, JNT 2018), Prop. 2.1: |theta(x)-x|/x < 1.2e-8 (x>e^45), 1.2e-9 (x>e^50),
2.9e-10 (x>e^55), from Buthe 2016 Table 1 (H = 3.061e10): |psi-x|/x <= 1.118e-8, 1.165e-9, 2.885e-10, plus
psi - theta <= 1.02 x^1/2 + 3 x^1/3.  Recompute those Table 1 rows with the zero-sum model (as printed, and
with log(3c) at Buthe's parameters, and re-optimised with log(3c), T <= 3.061e10), and the theta margins."""
import sys, math
sys.path.insert(0, '.')
from zeros_model import delta, optimise
H = 3.061e10
rows = [(45, 25, 3.5e9, 0.11, 1.11742e-8, 1.2e-8), (50, 30, H, 0.11, 1.16465e-9, 1.2e-9), (55, 30, H, 0.10, 2.88434e-10, 2.9e-10)]
for b, c, T, a, printed, target in rows:
    m0 = delta(b, c, T, a, False)[0]; m1 = delta(b, c, T, a, True)[0]
    r0 = optimise(b, False, H, c0=c, T0=T, a0=a); r1 = optimise(b, True, H, c0=c + 1, T0=T, a0=a)
    pt = 1.02 * math.exp(-b / 2) + 3 * math.exp(-2 * b / 3)
    fixed = printed * m1 / m0; reopt = r1[0]   # least corrected value (model reproduces the printed rows to <0.03%)
    print(f"e^{b}: printed {printed:.5e}; model at Buthe's params {m0:.5e} ({m0/printed-1:+.2%}); with log3c same params x{m1/m0:.5f};"
          f" model optimum {r0[0]:.5e} -> with log3c {r1[0]:.5e} (x{r1[0]/r0[0]:.5f}, c {r0[1]:.2f}->{r1[1]:.2f}, T {r0[2]:.3e}->{r1[2]:.3e})")
    for name, e in [("printed", printed), ("fixed-param corr.", fixed), ("re-opt corr.", reopt), ("+2.6%", 1.026 * printed), ("+8%", 1.08 * printed)]:
        th = e + pt
        print(f"    {name:18s} eps_psi={e:.5e}  eps_theta<= {th:.5e}  vs stated {target:.1e}: {'OK' if th < target else 'FAILS'} (room {1 - th/target:+.2%})")
