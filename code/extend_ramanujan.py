"""
Extensions of verify_ramanujan.py (same criterion (C), same interval arithmetic).

Mode F (unconditional, large L).  |psi(x)-x| <= eps x with eps the least of
  * BKLNW Table 8, rows b >= 2500 only (rows b <= 2000 come from Buthe 2016, Theorem 1 as
    printed, whose E2 lacks a factor log(3c); see buthe_bound.py);
  * FKS = Fiori-Kadiri-Swidinsky, J. Math. Anal. Appl. 527 (2023) 127426, Table 5:
        E_psi(x) <= eps(x0) for all x >= x0   (their Theorem 1.1);
  * FKS Corollary 1.4 (published version), valid for all x > 2:
        E_psi(x) < 9.22022 (log x)^{3/2} exp(-0.8476836 sqrt(log x)).
  All three are nonincreasing in x on the ranges used (Table 8 and Table 5 are decreasing
  columns; Cor. 1.4 is decreasing for log x > 12.6), so the value attached to
  [e^k, e^(k+1)] is valid on all of [e^k, x].  Nothing here uses Dudek-Platt's Lemma 2.1.

Tail.  For L >= L1 every term of the error, multiplied by the power of L it carries in the
  criterion, is nonincreasing (checked below), so one evaluation at L1 covers [L1, oo).

Gap map.  Best available unconditional eps against the eps the criterion needs, on the open range.
"""
import sys
from mpmath import iv, mp, mpf
import verify_ramanujan as V

I = V.I
A_FKS = I('9.22022')
C_FKS = I('0.8476836')

def classical(l):
    """FKS Cor. 1.4 bound for E_psi at log x = l (interval in, interval out)."""
    return A_FKS * l**I('1.5') * iv.exp(-C_FKS * iv.sqrt(l))

# FKS Table 5: (log x0, eps(x0)), printed to 5 significant digits -> bumped one unit.
FKS5 = [(1000, '1.3812e-9'), (2000, '2.4746e-12'), (2100, '1.3246e-12'), (2200, '7.0976e-13'),
        (2300, '3.8059e-13'), (2400, '2.0426e-13'), (2500, '1.0972e-13'), (2600, '5.8978e-14'),
        (2700, '3.1729e-14'), (2800, '1.7087e-14'), (2900, '9.2091e-15'), (3000, '4.9678e-15'),
        (4000, '1.1021e-17'), (5000, '2.9942e-20'), (6000, '1.2976e-22'), (7000, '8.5162e-25'),
        (8000, '7.7860e-27'), (9000, '9.2262e-29'), (10000, '1.3688e-30')]
FKS5 = [(b, V.bump(e)) for b, e in FKS5]

KMAX = 10050
def eps_seg(k):
    """upper bound for |psi(t)-t|/t valid for all t in [e^k, e^KMAX] (k integer >= 20)."""
    cands = [classical(I(k)).b]
    for a, b, e in V.SEG:                       # BKLNW Table 8 (decreasing, so valid to e^3500);
        if a >= 2500 and a <= k < b:            # only the rows computed by the Platt-Trudgian
            cands.append(e.b)                   # method: rows b <= 2000 use Buthe 2016 Thm 1 as printed
    rows = [e for b, e in FKS5 if b <= k]
    if rows: cands.append(rows[-1].b)
    return I(min(cands))

EPS = {k: eps_seg(k) for k in range(20, KMAX + 1)}
# C[n] = (1/e^n) * sum_{20<=k<n} EPS[k] * e^{k+1} / ((k+1)^2 (1-2/k))   (bound for int_{e^20}^{e^n})
C = {20: I(0)}
E_INV = iv.exp(I(-1))
for k in range(20, KMAX):
    C[k + 1] = C[k] * E_INV + EPS[k] / (I(k + 1)**2 * (1 - I(2) / k))

def d_F(l):
    N = int(mpf(l.a))
    assert 20 <= N < KMAX
    pt, it = V.psitheta_terms(l)
    point_eps = I(min(EPS[N].b, classical(l).b))
    point = (point_eps + pt) / l
    integ = V.K0 * iv.exp(-l) + C[N] * iv.exp(N - l) + EPS[N] / (l**2 * (1 - I(2) / N))
    return point + integ + it

# ---------------- tail: one evaluation at L1 covers all L >= L1 ----------------
def d_T(l):
    """majorant of |pi - J|/x using only the classical bound, split at t = e^{l/2}."""
    pt, it = V.psitheta_terms(l)
    point = (classical(l) + pt) / l
    c20 = classical(I(20))                      # Cor. 1.4 bound is decreasing for log t > 12.6
    integ = (V.K0 * iv.exp(-l)
             + c20 * iv.exp(-l / 2) / ((l / 2)**2 * (1 - I(2) / 20))
             + classical(l / 2) / (l**2 * (1 - I(4) / l)))
    return point + integ + it

def tail_margin(L1):
    """rigorous lower bound for L^6*(C) valid for every L >= L1, given the monotonicity facts."""
    L = I(L1)
    u_all = I([0, (1 / L).b])
    main_inf = V.horner(V.N6, u_all) / (1 - u_all)**V.n          # enclosure over all L >= L1
    r1, r0 = V.rho(L - 1), V.rho(L)
    Pu = V.Pf(1 / L)
    corr = (r1 / L - 2 * Pu * r0 - r0**2) * L**6
    d1, d0 = d_T(L - 1), d_T(L)
    err = (d1 / L + 2 * (Pu + r0) * d0 + d0**2) * L**6
    return (main_inf.a - abs(corr).b - err.b).a, main_inf, corr, err

def monotone_checks(L1):
    """the factors whose monotonicity tail_margin relies on, with the threshold beyond which each
    is nonincreasing:  L^5 classical(L): L > (13/C)^2;  L^4 classical(L/2): L > 2 (11/C)^2;
    L^k e^{-cL} pieces: L > k/c (at most 6/(1/2) = 12);  L^6 rho(L) ~ 12!/L^7."""
    c = mpf('0.8476836')
    need = max((13 / c)**2, 2 * (11 / c)**2, 13)
    assert L1 > need, (L1, need)
    # numeric sanity: L^6 d_T(L) decreasing along a geometric sequence
    vals = [(V.I(L)**6 * d_T(V.I(L))).b for L in [L1, 2 * L1, 4 * L1, 8 * L1, 64 * L1]]
    assert all(a > b for a, b in zip(vals, vals[1:])), vals
    return need

def margin_adaptive(a, b, depth=0):
    """True if the lower bound for L^6*(C) is positive on [a, b], bisecting up to 12 times."""
    if V.margin(I([a, b]), d_F)[0].a > 0: return True, 1
    if depth >= 12: return False, 1
    m = (a + b) / 2
    ok1, n1 = margin_adaptive(a, m, depth + 1)
    if not ok1: return False, n1
    ok2, n2 = margin_adaptive(m, b, depth + 1)
    return ok2, n1 + n2

def run_F(lo, hi, step):
    """like V.run, but with d_F (V.run only knows d_U and d_Q), and bisecting failed intervals."""
    last_ok = lo; nint = 0; nsplit = 0
    for a, b in V.grid(lo, hi, step, 'F'):
        ok, k = margin_adaptive(mpf(a), mpf(b))
        nint += k; nsplit += (k > 1)
        if not ok:
            print(f"  [F] FAIL on L in [{a}, {b}]")
            return last_ok, nint, nsplit
        last_ok = b
    return last_ok, nint, nsplit

# ---------------- gap map ----------------
def eps_need(L):
    """eps for which 3 eps L^4 = 1 + 14/L: the constant-eps approximation of (C).  Slightly
    optimistic: with eps(600,700) it gives 643.28, against the rigorous U-mode endpoint 642.84."""
    L = mpf(L)
    return (1 + 14 / L) / (3 * L**4)

if __name__ == "__main__":
    print(f"classical(20) = {float(classical(I(20)).b):.6g};  EPS[3000] = {float(EPS[3000].b):.6g}"
          f";  EPS[3158] = {float(EPS[3158].b):.6g}")
    # 1. where does Mode F start to work?  scan down from 3200, then refine
    last_fail = None
    for a, b in V.grid(3000, 3200, 0.5, 'F'):
        if V.margin(I([a, b]), d_F)[0].a <= 0: last_fail = b
    print(f"[F] coarse: last failing interval ends at L = {last_fail}")
    lo = float(last_fail) - 3
    Lstar = lo
    for a, b in V.grid(lo, float(last_fail) + 0.5, 0.002, 'F'):
        if V.margin(I([a, b]), d_F)[0].a <= 0: Lstar = b
    print(f"[F] fine: last failing interval ends at L = {Lstar}")
    # 2. contiguous verification from L* to L1
    L1 = 10000
    ok, cnt, nsplit = run_F(Lstar, L1, 0.5)
    print(f"[F] verified contiguously for L in [{Lstar}, {ok}]  ({cnt} intervals, {nsplit} of the "
          f"0.5-wide ones bisected)")
    # 3. tail
    need = monotone_checks(L1)
    tm, main_inf, corr, err = tail_margin(L1)
    print(f"[tail] monotonicity holds for L > {mp.nstr(need, 6)}; for all L >= {L1}: margin >= "
          f"{float(mpf(tm)):.10f}  (main >= {float(mpf(main_inf.a)):.10f}, |corr| <= "
          f"{float(mpf(abs(corr).b)):.3g}, err <= {float(mpf(err.b)):.3g})")
    # 4. gap map on the open range
    print("[gap] L      best eps        needed eps     shortfall")
    for L in [643, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3100, 3150]:
        best = float(mpf(min(EPS[L].b, classical(I(L)).b)))
        need_L = float(eps_need(L))
        print(f"      {L:5d}  {best:12.4g}  {need_L:12.4g}  {best / need_L:8.3g}x")
