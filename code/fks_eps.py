"""
Re-implementation of Fiori-Kadiri-Swidinsky, Sharper bounds for the Chebyshev function psi(x),
J. Math. Anal. Appl. 527 (2023) 127426 (arXiv:2204.02588), Theorem 1.1:

  for log x0 > 1000:  E_psi(x) <= eps(x0) := eps1(x0,T) + eps2(x0,0.9,T0,T) + eps3(x0,0.9,s2,N,T)
                                             + eps4(x0,s2,K,T)   for all x >= x0,
  T = exp(c sqrt(log x0 / R)), c = 3, H0 = 3e12, T0 = 30 610 046 000, R = 5.5666305.

Ingredients (equation numbers of FKS):
  eps1 = 2 (log x)^2 / T                                                           (3.4)
  eps2 = 2 x^{-1/2} (S0 + B1(T0,T)) + (x^{s1-1} - x^{-1/2}) B1(H0,T)                (3.8)
  eps3 = 2 x^{-(1-s1)+(s2-s1)/N} B0(s1,H_{s1},T)
         + 2 x^{-(1-s1)} (1 - x^{-(s2-s1)/N}) sum_{n=1}^{N-1} B0(s(n),H(n),T) x^{(s2-s1)(n+1)/N}   (3.17)
  eps4 = 2 sum_{k=1}^{K-1} x^{-1/(R log t_k)}/t_k (Nt(s2,t_{k+1}) - Nt(s2,t_k))
         + 2 Nt(s2,t_1) x^{-1/(R log t_0)}/t_0                                     (3.22)
  t0 = max(H_{s2}, exp(sqrt(log x / R))), t_k = t0 (T/t0)^{k/K},  H_s = max(H0, exp(1/(R(1-s))))
  B1, B0 as in (2.1), (2.10); Rosser's b1 = 0.137, b2 = 0.443, b3 = 1.588; S0 = 39.5797647802.
  H0v: the RH verification height used for the zero sums (default 3e12).  A larger verified
  height (Theorem D) may be inserted here; the density constants of Table 7/8 were computed with
  H0 = 3e12, which remains a verification height, so they stay valid.
  Zero density (Cor. 2.9, Table 8, subconvexity a1 = 0.77):  N(s,T) <= Nt(s,T) =
         c1 T^{p} (log T)^{q} + c2 (log T)^2,  p = 8(1-s)/3, q = 5 - 2s, (c1,c2) piecewise in s.

Evaluated in float64 (scipy) with a relative safety factor; validated against FKS Table 5.
"""
import math
import numpy as np
from scipy.special import gammaincc, gamma as Gamma

H0 = 3e12
T0, S0 = 30610046000.0, 39.5797647802
B1c, B2c, B3c = 0.137, 0.443, 1.588
SAFETY = 1 + 1e-9

# FKS Table 8 (a1 = 0.77): (sigma_lo, sigma_hi, c1, c2)
TABLE8 = [(0.90, 0.91, 12.5220, 2.8915), (0.91, 0.92, 13.0287, 2.8112), (0.92, 0.93, 13.5468, 2.7309),
          (0.93, 0.94, 14.0757, 2.6506), (0.94, 0.95, 14.6145, 2.5702), (0.95, 0.96, 15.1623, 2.4897),
          (0.96, 0.97, 15.7181, 2.4093), (0.97, 0.98, 16.2805, 2.3287), (0.98, 0.99, 16.8481, 2.2481),
          (0.99, 1.00, 17.4194, 2.1675)]

# FKS Table 7 (a1 = 0.77): per-sigma optimised (c1, c2), valid at that sigma (Thm 2.7 with k = 1)
TABLE7 = {0.990: (16.8481, 2.1648), 0.991: (16.9051, 2.1567), 0.992: (16.9621, 2.1486),
          0.993: (17.0192, 2.1406), 0.994: (17.0762, 2.1325), 0.995: (17.1334, 2.1244),
          0.996: (17.1905, 2.1163), 0.997: (17.2477, 2.1083), 0.998: (17.3049, 2.1002),
          0.999: (17.3621, 2.0921)}

# FKS Table 7, second block (subconvexity a1 = 0.618, Hiary-Patel-Yang, J. Number Theory 2024)
TABLE7_618 = {0.990: (16.3799, 2.1437), 0.991: (16.4495, 2.1357), 0.992: (16.5192, 2.1277),
              0.993: (16.5890, 2.1197), 0.994: (16.6590, 2.1117), 0.995: (16.7292, 2.1037),
              0.996: (16.7994, 2.0958), 0.997: (16.8698, 2.0878), 0.998: (16.9404, 2.0798),
              0.999: (17.0110, 2.0718)}
# Bellotti, An explicit log-free zero density estimate, arXiv:2405.12545, Thm 1.1 / Table 2:
# N(sigma,T) <= C T^{B(1-sigma)} for T in (T0, T1], sigma >= alpha0   (rows: logT0, logT1, alpha0, C, B)
BELLOTTI = [(math.log(3e12), 29, 0.9927, 37341.72, 14.160), (29, 30, 0.9927, 37564.18, 14.084),
            (30, 31, 0.9925, 40014.59, 13.588), (31, 32, 0.9923, 42840.01, 13.108),
            (32, 33, 0.9921, 46015.02, 12.647), (33, 34, 0.9919, 49270.71, 12.226),
            (34, 35, 0.9918, 52249.66, 11.856), (35, 36, 0.9916, 55524.80, 11.505),
            (36, 37, 0.9914, 59715.72, 11.134), (37, 38, 0.9912, 64228.93, 10.782),
            (38, 39, 0.9910, 68446.65, 10.475), (39, 40, 0.9909, 71907.32, 10.217),
            (40, 41, 0.9907, 76267.80, 9.951), (41, 42, 0.9906, 80004.48, 9.720),
            (42, 43, 0.9905, 83811.09, 9.501), (43, 44, 0.9903, 89209.33, 9.259),
            (44, 45, 0.9901, 95233.68, 9.020), (45, 46.2, 0.9900, 99227.30, 8.843),
            (46.2, 50, 0.9899, 104761.04, 8.634), (50, 60, 0.9892, 128806.85, 7.936),
            (60, 70, 0.9876, 194323.75, 6.661), (70, 80, 0.9860, 278045.07, 5.809),
            (80, 90, 0.9850, 370655.73, 5.216), (90, 100, 0.9850, 425721.47, 4.831),
            (100, 110, 0.9850, 488901.14, 4.513), (110, 120, 0.9850, 545744.21, 4.264),
            (120, 130, 0.9850, 629490.27, 4.032), (130, 140, 0.9850, 694045.43, 3.855),
            (140, 150, 0.9850, 771373.78, 3.696), (150, 160, 0.9850, 825913.07, 3.572),
            (160, 170.2, 0.9850, 909966.07, 3.448), (170.2, 500, 0.9850, 1.12e6, 3.337)]
DENSITY = {"mode": "KLN"}   # "KLN" (a1 = 0.77), "KLN618" (a1 = 0.618), "BEL618" (min of KLN618 and Bellotti)

BUMP = 1e-4          # one unit in the last printed digit of FKS Tables 7/8 (printed to 4 decimals)
TABLE7 = {k: (a + BUMP, b + BUMP) for k, (a, b) in TABLE7.items()}
TABLE7_618 = {k: (a + BUMP, b + BUMP) for k, (a, b) in TABLE7_618.items()}
TABLE8 = [(lo, hi, c1 + BUMP, c2 + BUMP) for lo, hi, c1, c2 in TABLE8]
BELLOTTI = [(l0, l1, a0, C * (1 + 1e-6), B + 0.001) for l0, l1, a0, C, B in BELLOTTI]   # printed to 2-3 decimals

def bellotti(s, lt):
    """Bellotti bound at log T = lt (numpy array); inf where not applicable."""
    out = np.full_like(lt, np.inf, dtype=float)
    for l0, l1, a0, C, B in BELLOTTI:
        if s >= a0:
            m = (lt > l0) & (lt <= l1)
            out[m] = C * np.exp(B * (1 - s) * lt[m])
    return out

def c12(s):
    """density constants valid at sigma = s (if s is on a boundary, the larger pair)."""
    best = None
    for lo, hi, c1, c2 in TABLE8:
        if lo - 1e-15 <= s <= hi + 1e-15:
            cand = (c1, c2)
            best = cand if best is None else (max(best[0], c1), max(best[1], c2))
    assert best is not None, s
    return best

def Rz(U): return B1c * math.log(U) + B2c * math.log(math.log(U)) + B3c

def B1(U, V):
    lu = math.log(U)
    return ((1 / (2 * math.pi) + (B1c * lu + B2c) / (U * lu * math.log(U / (2 * math.pi))))
            * math.log(V / U) * math.log(math.sqrt(V * U) / (2 * math.pi)) + 2 * Rz(U) / U)

def uGamma(a, x):  # upper incomplete gamma
    return gammaincc(a, x) * Gamma(a)

def B0(s, U, V):
    c1, c2 = c12(s)
    p, q = 8 * (1 - s) / 3, 5 - 2 * s
    lV, lU = math.log(V), math.log(U)
    return (c1 * lV**q / V**(1 - p) + c2 * lV**2 / V
            + c1 / (1 - p)**(q + 1) * (uGamma(q + 1, (1 - p) * lU) - uGamma(q + 1, (1 - p) * lV))
            + c2 * (uGamma(3, lU) - uGamma(3, lV)))

def Nt(s, t):
    key = round(s, 3)
    on_grid = abs(s - key) < 1e-12 and key in TABLE7
    mode = DENSITY["mode"]
    if on_grid:
        c1, c2 = (TABLE7_618 if mode in ("KLN618", "BEL618") else TABLE7)[key]
    else:
        c1, c2 = c12(s)
    p, q = 8 * (1 - s) / 3, 5 - 2 * s
    lt = np.log(t)
    kln = c1 * np.exp(p * lt) * lt**q + c2 * lt**2
    if mode == "BEL618":
        return np.minimum(kln, bellotti(s, np.atleast_1d(lt).astype(float)))
    return kln

def Hs(s, R, H0v=H0): return max(H0v, math.exp(1 / (R * (1 - s)))) if s < 1 else math.inf

def eps_parts(L, s2, R=5.5666305, c=3, N=2000, K=20000, s1=0.9, H0v=H0, Lmax=None):
    """FKS eps(x0) components at log x0 = L for given sigma_2."""
    T = math.exp(c * math.sqrt(L / R))
    assert 3 * L < T and math.log(T) < L / 2 - math.log(3) and L > 50
    e1 = 2 * L**2 / T
    xm12 = math.exp(-L / 2)
    e2 = 2 * xm12 * (S0 + B1(T0, T)) + (math.exp((s1 - 1) * L) - xm12) * B1(H0v, T)
    # eps3
    d = (s2 - s1) / N
    e3 = 2 * math.exp(-(1 - s1) * L + d * L) * B0(s1, Hs(s1, R, H0v), T)
    acc = 0.0
    for n in range(1, N):
        sn = s1 + n * d
        acc += B0(sn, Hs(sn, R, H0v), T) * math.exp(d * (n + 1) * L)
    e3 += 2 * math.exp(-(1 - s1) * L) * (1 - math.exp(-d * L)) * acc
    # eps4
    # FKS Prop 3.11 needs t0 >= exp(sqrt(log x / R)) at every x where the row is used (x <= e^Lmax)
    t0 = max(Hs(s2, R, H0v), math.exp(math.sqrt((L if Lmax is None else max(L, Lmax)) / R)))
    assert T > t0
    k = np.arange(0, K + 1)
    tk = t0 * (T / t0) ** (k / K)
    w = np.exp(-L / (R * np.log(tk))) / tk
    Ntk = Nt(s2, tk)
    e4 = 2 * np.sum(w[1:K] * (Ntk[2:K + 1] - Ntk[1:K])) + 2 * Ntk[1] * w[0]
    return e1, e2, e3, float(e4)

def eps(L, s2, **kw):
    return SAFETY * sum(eps_parts(L, s2, **kw))

def best_eps(L, lo=0.980, hi=0.9985, **kw):
    """minimise eps over sigma_2 (golden section; eps is unimodal in sigma_2 in practice)."""
    g = (math.sqrt(5) - 1) / 2
    a, b = lo, hi
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = eps(L, x1, **kw), eps(L, x2, **kw)
    for _ in range(40):
        if f1 < f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - g * (b - a); f1 = eps(L, x1, **kw)
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + g * (b - a); f2 = eps(L, x2, **kw)
    return (f1, x1) if f1 < f2 else (f2, x2)

def best_eps_table7(L, **kw):
    """eps(x0) minimised over sigma_2 in {0.990, ..., 0.999} (Table 7 constants in eps4)."""
    out = []
    for s2 in sorted(TABLE7):
        try:
            out.append((eps(L, s2, **kw), s2))
        except AssertionError:          # T <= t0 for this sigma_2: not admissible
            pass
    return min(out)

if __name__ == "__main__":
    print("Validation against FKS Table 5 (their sigma_2, their eps):")
    for L, s2, printed in [(2000, 0.99221, 2.4746e-12), (2500, 0.99245, 1.0972e-13),
                           (2900, 0.99260, 9.2091e-15), (3000, 0.99263, 4.9678e-15)]:
        parts = eps_parts(L, s2, N=10000, K=100000)
        print(f"  log x0 = {L}: ours {sum(parts):.5e} (e1 {parts[0]:.4e}, e2 {parts[1]:.4e}, "
              f"e3 {parts[2]:.4e}, e4 {parts[3]:.4e})   printed {printed:.4e}")
    print("With Table 7 constants (sigma_2 on the 0.001 grid):")
    for L, printed in [(2000, 2.4746e-12), (2500, 1.0972e-13), (2900, 9.2091e-15), (3000, 4.9678e-15)]:
        e, s2 = best_eps_table7(L, N=10000, K=100000)
        print(f"  log x0 = {L}: ours {e:.5e} at sigma_2 = {s2}   printed {printed:.4e}")
