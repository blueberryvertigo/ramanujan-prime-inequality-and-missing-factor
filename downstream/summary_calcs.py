"""Margin checks for the downstream report (uses e2_variants.py; values printed are the least
re-optimised bounds unless stated).  Sections:
 A. Johnston-Yang Lemma 2.3 / Theorem 1.1 crossing below e^2488.
 B. Lichtman-Pomerance Prop. 2.1 (theta bounds at e^45, e^50, e^55; RH height 3.061e10).
 C. BKLNW Table 15 maxima near b = 2300 and Cor. 11.1 (0.024334), on BKLNW's 25-step grid.
 D. Bordignon-Starichkova constant 4.47e-9 (= M1 (1 + 3/(2*1000)) for x > e^1000).
 E. FKS1 Lemma 5.2 margins (Corollary 1.4 in the Buthe range).
 F. JRT Lemma 8 (1.994e-8 for x >= e^40) and the BKLNW 10^19 row (Cor. 2.1's 1.93378e-8)."""
import sys, math
sys.path.insert(0, '.')
sys.stdout.reconfigure(encoding="utf-8")
from e2_variants import optimise, delta
from mpmath import mp, mpf, exp, sqrt, findroot
mp.dps = 30
H_PT, H_B1, H_BK = 3000175332800.0, 3.061e10, 2445999556030.0
V = ('printed', 'log3c', 'thesis')
print("A. Johnston-Yang: f(L) = 9.39 L^1.515 exp(-0.8274 sqrt L); Lemma 2.3 value 1.570e-12 (x >= e^2000)")
f = lambda L: mpf('9.39') * mpf(L)**mpf('1.515') * exp(-mpf('0.8274') * sqrt(mpf(L)))
cross = lambda e: findroot(lambda L: f(L) - mpf(e), 2485)
print(f"   f(2488) = {mp.nstr(f(2488), 7)};  printed 1.570e-12 covers L <= {mp.nstr(cross(1.570e-12), 7)}")
for v in V[1:]:
    r = optimise(2000, v, H_PT, c0=35, a0=0.003)
    rf = delta(2000, 33.6, 33.6 / 1.12e-11, 0.0, v)[0]
    print(f"   {v:7s}: at JY's params (c=33.6, alpha=0) {rf:.5e} -> covers L <= {mp.nstr(cross(rf), 7)};"
          f"  re-optimised {r[0]:.5e} -> covers L <= {mp.nstr(cross(r[0]), 7)}")
print("B. Lichtman-Pomerance Prop 2.1: eps_theta <= eps_psi + 1.02 e^{-b/2} + 3 e^{-2b/3}")
for b, c, T, a, target in [(45, 25, 3.5e9, 0.11, 1.2e-8), (50, 30, H_B1, 0.11, 1.2e-9), (55, 30, H_B1, 0.10, 2.9e-10)]:
    pt = 1.02 * math.exp(-b / 2) + 3 * math.exp(-2 * b / 3)
    for v in V:
        r = optimise(b, v, H_B1, c0=c, T0=T, a0=a)
        th = r[0] + pt
        print(f"   e^{b} {v:7s}: least eps_psi {r[0]:.5e}, eps_theta <= {th:.5e} vs {target:.1e}: {'OK' if th < target else 'FAILS'} ({1 - th/target:+.2%})")
print("C. BKLNW Table 15: B_k(2300,2325) = eps(2300)*2325^k (+ negligible), printed eps(2300) = 4.4627e-9/2325")
e_pr = 4.4627e-9 / 2325
B2300 = {1: 4.4627e-9, 2: 1.0376e-5, 3: 2.4124e-2, 4: 56.088, 5: 1.3040e5}
Mk = {1: 4.4627e-9, 2: 1.0376e-5, 3: 2.4334e-2, 4: 57.184, 5: 1.3441e5}
for v in V[1:]:
    r = optimise(2300, v, H_BK, c0=35, a0=0.0025)
    s = r[0] / e_pr
    print(f"   {v:7s}: eps(2300) {r[0]:.5e} (x{s:.5f}); " + "; ".join(
        f"k={k}: {B2300[k]*s:.5g} {'>' if B2300[k]*s > Mk[k] else '<='} M_k={Mk[k]:.5g}" for k in range(1, 6)))
print("D. Bordignon-Starichkova: M1 (1 + 1.5/1000) vs stated 4.47e-9")
for v in V[1:]:
    r = optimise(2300, v, H_BK, c0=35, a0=0.0025)
    m1 = 4.4627e-9 * r[0] / e_pr
    print(f"   {v:7s}: M1 -> {m1:.5e};  M1 (1+1.5e-3) = {m1 * 1.0015:.4e} vs 4.47e-9")
print("E. FKS1 Lemma 5.2: eps(e^1000) < 2 (2100)^1.5 exp(-0.8476836 sqrt 2100)")
g = lambda L: 2 * L**1.5 * math.exp(-0.8476836 * math.sqrt(L))
for v in V[1:]:
    r = optimise(1000, v, H_BK, c0=35, a0=0.005)
    print(f"   {v:7s}: eps(e^1000) {r[0]:.5e} vs g(2100) = {g(2100):.5e}: room {1 - r[0]/g(2100):+.1%};"
          f"  eps(e^40)=1.93378e-8 vs g(1000) = {g(1000):.4e}")
print("F. 10^19 row (BKLNW Cor 2.1 / Table 8 b=40 / JRT Lemma 8)")
for v in V:
    r = optimise(19 * math.log(10), v, H_BK, c0=30, T0=3.5e9, a0=0.12)
    print(f"   {v:7s}: least eps(10^19) {r[0]:.5e} (c={r[1]:.2f}, T={r[2]:.3e}); vs 1.93378e-8: {'<=' if r[0] <= 1.93378e-8 else '>'}")
