"""BKLNW Table 8 rows b = 19 log 10 (= the value carried by the b=40 row), 45, 50, 100, 200.
Buthe 2016 Thm 1, E3 bounded with ell <= 1 and Lemma 3, i.e.
    E3 <= 2/sqrt(x0) * log(T/2pi)^2/(4pi) + 2/x0,   T = c/eps free, 5000 <= T <= H.
Minimise over (c, T, alpha) with and without log(3c) (floats; mu_c, nu_c by mpmath quadrature, cached).
Descriptive: how much of each row the E2 term carries, and the corrected/uncorrected ratio."""
import math
from functools import lru_cache
from mpmath import mp, mpf, besseli, quad, sinh as msinh, sqrt as msqrt
mp.dps = 25
@lru_cache(maxsize=None)
def munu(c, a):
    c = mpf(c); a = mpf(a)
    eta = lambda t: c * besseli(0, c * msqrt(1 - t * t)) / (2 * msinh(c))
    m = mpf(1) / 2 - quad(eta, [0, a])
    n = besseli(1, c) / (2 * msinh(c)) - a / 2 + quad(lambda t: (a - t) * eta(t), [0, a])
    return float(m), float(n)
H = 2445999556030.0
def delta(b, c, T, a, log3c):
    eps = c / T; lx0 = b - a * eps; x0 = math.exp(lx0)
    m, n = munu(c, a); B0 = eps * x0 * math.exp(-eps) * n / (2 * m)
    if B0 <= 1: return None
    E1 = math.exp(2 * eps) * (eps + lx0) * (2 * eps * n / math.log(B0) + 2.01 * eps / math.sqrt(x0)
          + math.log(math.log(2) + 2 * lx0) / (2 * x0)) + math.expm1(a * eps)
    E2 = 0.16 * (1 + 1 / x0) / math.sinh(c) * math.exp(0.71 * math.sqrt(c * eps)) * (math.log(3 * c) if log3c else 1) * math.log(T)
    E3 = 2 / math.sqrt(x0) * math.log(T / (2 * math.pi))**2 / (4 * math.pi) + 2 / x0
    return math.exp(a * eps) * (E1 + E2 + E3), (E1, E2, E3)
def opt(b, log3c, cs, Ts, As):
    best = None
    for c in cs:
        for T in Ts:
            for a in As:
                r = delta(b, c, T, a, log3c)
                if r and (best is None or r[0] < best[0]): best = (r[0], c, T, a, r[1])
    return best
def refine(b, log3c):
    cs = [10 + 0.5 * i for i in range(51)]
    Ts = sorted(set([min(H, 10 ** (6 + 0.05 * k)) for k in range(140)]))
    As = [0.0, 0.0003, 0.001, 0.003, 0.01, 0.02, 0.03, 0.06, 0.1]
    r = opt(b, log3c, cs, Ts, As)
    c0, T0, a0 = r[1], r[2], r[3]
    cs = [c0 - 0.5 + 0.05 * i for i in range(21)]
    Ts = sorted(set([min(H, T0 * 10 ** (-0.05 + 0.005 * k)) for k in range(21)]))
    As = sorted(set([max(0.0, a0 * f) for f in (0.5, 0.7, 0.85, 1, 1.2, 1.4, 2)] + [0.0]))
    return opt(b, log3c, cs, Ts, As)
if __name__ == "__main__":
    for b, printed in [(19 * math.log(10), 1.93378e-8), (45, 1.09073e-8), (50, 1.11990e-9), (100, 2.45299e-12), (200, 2.18154e-12)]:
        r0 = refine(b, False); r1 = refine(b, True)
        print(f"b={b:.4f} printed {printed:.5e}: no-log3c {r0[0]:.5e} (c={r0[1]:.2f}, T={r0[2]:.3e}, a={r0[3]:.4g}; "
              f"E1={r0[4][0]:.3e} E2={r0[4][1]:.3e} E3={r0[4][2]:.3e});  with-log3c {r1[0]:.5e} (c={r1[1]:.2f}, T={r1[2]:.3e}, a={r1[3]:.4g}) "
              f"ratio corrected/printed {r1[0]/printed:.5f}", flush=True)
