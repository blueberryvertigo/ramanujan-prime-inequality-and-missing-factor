"""Re-optimise Buthe's (c, alpha) with log(3c) restored, evaluating mu_c, nu_c by 40-digit quadrature
(crosscheck_buthe.delta_quad; the interval enclosures in buthe_bound.py are loose for larger alpha and
would overstate the minimum).  Descriptive: how far above the printed values the corrected optimum lies."""
from mpmath import mp, mpf, besseli, quad, sinh, log, exp, sqrt, pi
mp.dps = 30
def eta(c, t): return c * besseli(0, c * sqrt(1 - t * t)) / (2 * sinh(c))
def mu(c, a): return mpf(1) / 2 - quad(lambda t: eta(c, t), [0, a])
def nu(c, a):   # |nu_c(a)| = |nu_c(0)| - int_0^a mu = |nu_c(0)| - a/2 + int_0^a (a - t) eta(t) dt
    return besseli(1, c) / (2 * sinh(c)) - mpf(a) / 2 + quad(lambda t: (a - t) * eta(c, t), [0, a])
def delta_quad(b, c, T, a, log3c=True):
    eps = mpf(c) / T; x0 = exp(b); m = mu(c, a); n = nu(c, a)
    B0 = eps * x0 * exp(-eps) * n / (2 * m)
    E1 = (exp(2 * eps) * (eps + b) * (2 * eps * n / log(B0) + mpf('2.01') * eps / sqrt(x0)
          + log(log(2) + 2 * b) / (2 * x0)) + (exp(a * eps) - 1))
    E2 = (mpf('0.16') * (1 + 1 / x0) / sinh(c) * exp(mpf('0.71') * sqrt(c * eps))
          * (log(3 * c) if log3c else 1) * log(mpf(c) / eps))
    E3 = 2 / sqrt(x0) * log(T / (2 * pi))**2 / (4 * pi) + 2 / x0
    return exp(a * eps) * (E1 + E2 + E3)
ROWS = [("Buthe Table 2, e^1000", 1000, 2.445e12, 1.94751e-12), ("Buthe Table 2, e^4000", 4000, 2.445e12, 1.90866e-12),
        ("BKLNW Table 8, b=200", 200, 2.445999556030e12, 2.18154e-12), ("BKLNW Table 8, b=600", 600, 2.445999556030e12, 1.98894e-12),
        ("Buthe Table 1, e^500", 500, 3.061e10, 1.47067e-10)]
for name, b, T, p in ROWS:
    best = None
    c0 = 30 if T < 1e11 else 33
    for c in [c0 + 0.25 * i for i in range(25)]:
        for j in range(31):
            a = 10 ** (-3.5 + 0.1 * j)             # 3.2e-4 ... 0.32
            d = float(delta_quad(b, c, T, a, log3c=True))
            if best is None or d < best[0]: best = (d, c, a)
    print(f"{name}: printed {p:.5e}; least corrected value {best[0]:.5e} at c={best[1]}, alpha={best[2]:.4g} "
          f"(+{100 * (best[0] / p - 1):.2f}%)", flush=True)
