"""FKS2 (Fiori-Kadiri-Swidinsky, arXiv:2206.12557, Res. Number Theory 2023), Table 7: ranges
1 <= log x <= L* on which E_pi(x) <= B(x).  The upper ends L* are where B meets the numerical
bounds eps_pi,num of Table 4, which for 10^19 <= x < e^2073 rest on Buthe 2016 Thm 1 without log(3c).
(1) reproduce each L* from Table 4 (log-linear interpolation);  (2) add to eps_pi the increase of
eps_psi caused by restoring log(3c) -- at the uncorrected optimum parameters ('fixed') and after
re-optimising (c, alpha) ('reopt') -- and recompute L*.  Uses zeros_model.py (non-rigorous zero sum;
reproduces FKS2 Table 2 to <0.2%).  Also a uniform +2.6% and +8% perturbation."""
import math, sys
sys.path.insert(0, '.')
from zeros_model import delta, optimise
H = 3.0e12
T4 = [(44,1.7893e-8),(45,1.1449e-8),(46,7.2959e-9),(47,4.6388e-9),(48,2.9451e-9),(49,1.8680e-9),(50,1.1785e-9),
(51,7.4479e-10),(52,4.7046e-10),(53,2.9707e-10),(54,1.8753e-10),(55,1.1785e-10),(56,7.4191e-11),(57,4.6692e-11),
(58,2.9380e-11),(59,1.8774e-11),(60,1.2330e-11),(61,8.4134e-12),(62,6.0325e-12),(63,4.5827e-12),(64,3.6978e-12),
(65,3.1557e-12),(66,2.8216e-12),(67,2.6138e-12),(68,2.4828e-12),(69,2.3985e-12),(70,2.3427e-12),(71,2.3043e-12),
(72,2.2766e-12),(73,2.2555e-12),(74,2.2387e-12),(75,2.2244e-12),(76,2.2120e-12),(77,2.2006e-12),(78,2.1901e-12),
(79,2.1802e-12),(80,2.1708e-12),(81,2.1617e-12),(82,2.1530e-12),(83,2.1446e-12),(84,2.1364e-12),(85,2.1284e-12),
(90,2.0919e-12),(95,2.0598e-12),(100,2.0339e-12),(110,1.9853e-12),(120,1.9457e-12),(130,1.9126e-12),(140,1.8847e-12),
(150,1.8608e-12),(200,1.7789e-12),(250,1.7311e-12),(260,1.7238e-12),(270,1.7171e-12),(280,1.7108e-12),(300,1.6997e-12),
(400,1.6609e-12),(500,1.6380e-12),(1000,1.5925e-12),(1300,1.5822e-12),(1325,1.5815e-12),(1350,1.5809e-12),
(1375,1.5803e-12),(1400,1.5797e-12),(1500,1.5776e-12),(2000,1.5701e-12),(2100,1.3254e-12),(2200,7.1013e-13),
(3000,4.9698e-15),(3700,6.8015e-17),(3800,3.7027e-17),(4000,1.1024e-17)]
def eps_pi(L):
    for (a, ea), (b, eb) in zip(T4, T4[1:]):
        if a <= L <= b:
            return math.exp(math.log(ea) + (L - a) / (b - a) * (math.log(eb) - math.log(ea)))
    raise ValueError(L)
ROWS = [("2 log x x^-1/2", lambda L: 2 * L * math.exp(-L / 2), 57),
        ("log^1.5 x x^-1/2", lambda L: L**1.5 * math.exp(-L / 2), 65.65),
        ("(1/8pi) log^2 x x^-1/2", lambda L: L**2 * math.exp(-L / 2) / (8 * math.pi), 60.8),
        ("log^2 x x^-1/2", lambda L: L**2 * math.exp(-L / 2), 70.6),
        ("log^3 x x^-1/2", lambda L: L**3 * math.exp(-L / 2), 80),
        ("x^-1/3", lambda L: math.exp(-L / 3), 80.55),
        ("x^-1/4", lambda L: math.exp(-L / 4), 107.6),
        ("x^-1/5", lambda L: math.exp(-L / 5), 134.8),
        ("x^-1/10", lambda L: math.exp(-L / 10), 270.8),
        ("x^-1/50", lambda L: math.exp(-L / 50), 1358.6),
        ("x^-1/100", lambda L: math.exp(-L / 100), 3757.6)]
def crossing(B, ef, lo=44.0, hi=3999.0):
    """largest L in [lo, hi] such that B(l) >= ef(l) for all l in [lo, L] (scan + bisection)."""
    L = lo; step = 0.01
    if B(L) < ef(L): return None
    while L + step <= hi and B(L + step) >= ef(L + step): L += step
    a, b = L, L + step
    for _ in range(40):
        m = (a + b) / 2
        if B(m) >= ef(m): a = m
        else: b = m
    return a
_cache = {}
def deltas(L):
    """additive increase of eps_psi at log x = L: (fixed-parameter, re-optimised)"""
    Lr = round(L, 1)
    if Lr in _cache: return _cache[Lr]
    r0 = optimise(Lr, False, H, c0=33.3, T0=H, a0=max(0.002, min(0.12, 6.0 / Lr)))
    d0, c, T, a, _ = r0
    d_fixed = delta(Lr, c, T, a, True)[0] - d0
    r1 = optimise(Lr, True, H, c0=c + 1.4, T0=T, a0=a)
    _cache[Lr] = (d_fixed, r1[0] - d0, d0, r1[0], c, r1[1])
    return _cache[Lr]
if __name__ == "__main__":
    print(f"{'B(x)':24s} {'stated':>8s} {'repro':>9s} | {'fixed':>9s} {'reopt':>9s} | {'+2.6%':>9s} {'+8%':>9s} | dpsi_fixed/eps  dpsi_reopt/eps")
    for name, B, stated in ROWS:
        Ls = crossing(B, eps_pi)
        dfx, dre, d0, d1, c0, c1 = deltas(Ls)
        Lf = crossing(B, lambda l: eps_pi(l) + dfx)
        Lr = crossing(B, lambda l: eps_pi(l) + dre)
        L26 = crossing(B, lambda l: 1.026 * eps_pi(l))
        L8 = crossing(B, lambda l: 1.08 * eps_pi(l))
        print(f"{name:24s} {stated:8.2f} {Ls:9.3f} | {Lf:9.3f} {Lr:9.3f} | {L26:9.3f} {L8:9.3f} | "
              f"{dfx/eps_pi(Ls):+.2%} {dre/eps_pi(Ls):+.2%}  (model eps_psi {d0:.4e} -> {d1:.4e}; c {c0:.2f}->{c1:.2f})", flush=True)
