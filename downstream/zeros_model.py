"""Buthe 2016 Thm 1 with E3 = 2/sqrt(x0) * sum_{0<gamma<c/eps} ell_{c,eps}(gamma)/gamma + 2/x0, the zero sum
approximated by the first 10 zeros exactly plus the Riemann-von Mangoldt density (1/2pi) log(t/2pi) dt
beyond t = 50 (S(t) neglected; error well below 1% of the sum).  NOT rigorous: used only to reproduce
the small-b rows (b < 60) of BKLNW Table 8 / FKS2 Table 2 and to measure the size of the log(3c) effect."""
import math
from functools import lru_cache
from mpmath import mp, mpf, besseli, quad, sinh as msinh, sqrt as msqrt, log as mlog, pi as mpi
mp.dps = 20
G10 = [14.134725142, 21.022039639, 25.010857580, 30.424876126, 32.935061588,
       37.586178159, 40.918719012, 43.327073281, 48.005150881, 49.773832478]
@lru_cache(maxsize=None)
def munu(c, a):
    c = mpf(c); a = mpf(a)
    eta = lambda t: c * besseli(0, c * msqrt(1 - t * t)) / (2 * msinh(c))
    m = mpf(1) / 2 - quad(eta, [0, a])
    n = besseli(1, c) / (2 * msinh(c)) - a / 2 + quad(lambda t: (a - t) * eta(t), [0, a])
    return float(m), float(n)
def ell(c, eps, g):
    u = g * eps
    if u >= c: return 0.0
    r = math.sqrt(c * c - u * u)
    # c/sinh(c) * sinh(r)/r, computed stably for large c
    return (c / r) * math.exp(r - c) * (1 - math.exp(-2 * r)) / (1 - math.exp(-2 * c)) if r > 1e-8 else c / math.sinh(c)
@lru_cache(maxsize=None)
def zsum(c, T):
    eps = c / T
    s = sum(ell(c, eps, g) / g for g in G10)
    f = lambda t: ell(c, eps, float(t)) / float(t) * math.log(float(t) / (2 * math.pi)) / (2 * math.pi)
    # integrate in log t
    a, b = math.log(50.0), math.log(T)
    n = 4000; h = (b - a) / n
    tot = 0.0
    for i in range(n + 1):
        t = math.exp(a + i * h); w = 1 if i in (0, n) else (4 if i % 2 else 2)
        tot += w * f(t) * t
    return s + tot * h / 3
def delta(b, c, T, a, log3c, x0_is_eb=False):
    """log-space version (valid for large b): x0 = e^{lx0}, with e^{alpha eps} x0 = e^b unless x0_is_eb."""
    eps = c / T; lx0 = b if x0_is_eb else b - a * eps
    m, n = munu(c, a)
    logB0 = math.log(eps) + lx0 - eps + math.log(n / (2 * m))
    if logB0 <= 0: return None
    ism = math.exp(-lx0 / 2); ix0 = math.exp(-lx0)
    E1 = math.exp(2 * eps) * (eps + lx0) * (2 * eps * n / logB0 + 2.01 * eps * ism
          + math.log(math.log(2) + 2 * lx0) * ix0 / 2) + math.expm1(a * eps)
    E2 = 0.16 * (1 + ix0) / math.sinh(c) * math.exp(0.71 * math.sqrt(c * eps)) * (math.log(3 * c) if log3c else 1) * math.log(T)
    E3 = 2 * ism * zsum(c, T) + 2 * ix0
    return math.exp(a * eps) * (E1 + E2 + E3), (E1, E2, E3)
def optimise(b, log3c, Tmax, c0=25.0, T0=3e9, a0=0.12):
    best = None
    cs = [c0 + 0.5 * i for i in range(-16, 17)]
    Ts = [min(Tmax, T0 * 10 ** (0.05 * k)) for k in range(-20, 21)]
    As = [a0 * f for f in (0.5, 0.75, 1.0, 1.25, 1.5)]
    for rnd in range(3):
        for c in cs:
            for T in sorted(set(Ts)):
                for a in As:
                    r = delta(b, c, T, a, log3c)
                    if r and (best is None or r[0] < best[0]): best = (r[0], c, T, a, r[1])
        _, c0, T0, a0, _ = best
        step = 0.5 / 4 ** (rnd + 1)
        cs = [c0 + step * i for i in range(-8, 9)]
        Ts = [min(Tmax, T0 * 10 ** (0.05 / 4 ** (rnd + 1) * k)) for k in range(-8, 9)]
        As = [a0 * (1 + 0.25 / 4 ** (rnd + 1) * k) for k in range(-4, 5)]
    return best
