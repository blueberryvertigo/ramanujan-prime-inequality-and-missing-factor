"""Cross-check of fks_eps.py (float64) against a 40-digit mpmath evaluation of the same formulas."""
import math
from mpmath import mp, mpf, gammainc, log, exp, sqrt, pi
import fks_eps as F
mp.dps = 40
def B1m(U, V):
    U, V = mpf(U), mpf(V); lu = log(U)
    Rz = F.B1c * lu + F.B2c * log(lu) + F.B3c
    return ((1 / (2 * pi) + (F.B1c * lu + F.B2c) / (U * lu * log(U / (2 * pi))))
            * log(V / U) * log(sqrt(V * U) / (2 * pi)) + 2 * Rz / U)
def B0m(s, U, V):
    c1, c2 = F.c12(s); p, q = mpf(8) * (1 - s) / 3, 5 - 2 * mpf(s)
    lV, lU = log(V), log(U)
    return (c1 * lV**q / mpf(V)**(1 - p) + c2 * lV**2 / V
            + c1 / (1 - p)**(q + 1) * (gammainc(q + 1, (1 - p) * lU) - gammainc(q + 1, (1 - p) * lV))
            + c2 * (gammainc(3, lU) - gammainc(3, lV)))
def Ntm(s, t):
    key = round(s, 3); c1, c2 = F.TABLE7[key]
    p, q = mpf(8) * (1 - s) / 3, 5 - 2 * mpf(s); lt = log(t)
    return c1 * exp(p * lt) * lt**q + c2 * lt**2
def eps_m(L, s2, R=5.5666305, c=3, N=400, K=4000, s1=0.9, H0v=3e12):
    L, R = mpf(L), mpf(R)
    T = exp(c * sqrt(L / R))
    e1 = 2 * L**2 / T
    e2 = 2 * exp(-L / 2) * (F.S0 + B1m(F.T0, T)) + (exp((s1 - 1) * L) - exp(-L / 2)) * B1m(H0v, T)
    Hs = lambda s: max(mpf(H0v), exp(1 / (R * (1 - s))))
    d = (mpf(s2) - s1) / N
    e3 = 2 * exp(-(1 - s1) * L + d * L) * B0m(s1, Hs(s1), T)
    acc = mpf(0)
    for n in range(1, N):
        sn = s1 + n * d; acc += B0m(float(sn), Hs(sn), T) * exp(d * (n + 1) * L)
    e3 += 2 * exp(-(1 - s1) * L) * (1 - exp(-d * L)) * acc
    t0 = max(Hs(mpf(s2)), exp(sqrt(L / R)))
    lam = (T / t0) ** (mpf(1) / K)
    tk = [t0 * lam**k for k in range(K + 1)]
    w = [exp(-L / (R * log(t))) / t for t in tk]
    Nk = [Ntm(s2, t) for t in tk]
    e4 = 2 * sum(w[k] * (Nk[k + 1] - Nk[k]) for k in range(1, K)) + 2 * Nk[1] * w[0]
    return e1 + e2 + e3 + e4
for L, s2, R in [(3000, 0.993, 5.5666305), (2600, 0.993, 4.896)]:
    fl = sum(F.eps_parts(L, s2, N=400, K=4000, R=R))
    mm = eps_m(L, s2, R=R)
    print(f"L={L} s2={s2} R={R}: float {fl:.12e}  mp {mp.nstr(mm, 13)}  rel.diff {float(abs(fl - mm) / mm):.2e}")
