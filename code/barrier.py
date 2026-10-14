"""
Barrier for Ramanujan's inequality from spectral data up to height T (sharpness of Theorem D).

Class A_T: sequences a_n >= 0 whose Dirichlet series A(s) = sum a_n n^{-s} converges for Re s > 1,
continues meromorphically to C, has a simple pole at s = 1 with residue 1 and no other pole with
|Im s| <= T.  [CH, Thm 1.1] applies uniformly to A_T.  Put psi_a(y) = sum_{n<=y} a_n and
pi_a(y) = sum_{2<=n<=y} a_n / log n.

Construction (Chirre-Helfgott, Prop. 3.1): a_n = F_K(T' log n), F_K the Fejer kernel, T' > T, so
A(s) = sum_{|k|<=K} (1-|k|/(K+1)) zeta(s - ikT') is in A_T.  With u = pi/T', at
log y = (2 pi N + delta)/T' (just after a peak), psi_a(y)/y -> u(coth u + 1) e^{-delta/T'} as K -> oo,
and at log y = (2 pi N - delta)/T' (just before a peak) psi_a(y)/y -> u(coth u - 1) e^{delta/T'}.
Taking T' = 2 pi M + 2 delta puts x just after a peak and x/e just before one.

Partial summation: pi_a(y) = J(y) + (S(y) - 1) y / log y + O*(R(y)), J(y) = y/log y + int_2^y dt/log^2 t,
with S(y) = psi_a(y)/y up to O(K^2 T' log y + K)/y and |R(y)| <= (2 pi^2/3) y / (T'^2 log^2 y) + o(y)
(the oscillating part of int_2^y psi_a(t)/(t log^2 t) dt; see the paper).

Then, with L = log x and j(l) = J(e^l)/e^l, the fake violates (R) at x as soon as
   Main(L) + L^6 [ s_minus / (L (L-1)) - 2 j(L) s_plus / L - s_plus^2 / L^2 ] + L^6 * 3 R/x^2 < 0,
Main(L) = L^6 (j(L-1)/L - j(L)^2) (computed exactly as in verify_ramanujan.py),
s_plus = u(coth u + 1) - 1, s_minus = u(coth u - 1) - 1.

T_bar(L) = sup of T' for which this holds.  T_suf(L) = least T for which [CH, Prop 9.1] + (C)
proves (R) on a neighbourhood of L (computed with ch_ramanujan.py's machinery).
"""
from mpmath import mp, mpf, iv, coth, findroot
import verify_ramanujan as V
import ch_ramanujan as CH

mp.dps = 40
I = V.I

def main_exact(L):
    """Main(L) = L^6 (j(L-1)/L - j(L)^2), enclosed exactly (main + corr in V.margin)."""
    Li = I(L)
    u = 1 / Li
    main = V.horner(V.N6, u) / (1 - u)**V.n
    r1, r0 = V.rho(Li - 1), V.rho(Li)
    Pu = V.Pf(u)
    corr = (r1 / Li - 2 * Pu * r0 - r0**2) * Li**6
    return main + corr, (V.Pf(u) + r0)          # also j(L) enclosure

def fake_margin(L, Tp):
    """L^6 * (R-margin)/x^2 for the Fejer fake with parameter T' = Tp (K -> oo, delta -> 0)."""
    L, Tp = mpf(L), mpf(Tp)
    M, jL = main_exact(L)
    M = mpf(M.b); jL = mpf(jL.a)               # upper end of Main, lower end of j: an upper bound
                                               # for the margin (the j-term enters with a minus sign)
    u = mp.pi / Tp
    sp = u * (coth(u) + 1) - 1
    sm = u * (coth(u) - 1) - 1
    R = 2 * mp.pi**2 / 3 / (Tp**2 * L**2)      # >= 2 x the bound for |Osc(y)|/y + |Err(y)|/y (see proof)
    return M + L**6 * (sm / (L * (L - 1)) - 2 * jL * sp / L - sp**2 / L**2) + L**6 * 3 * R / L

def T_bar(L):
    f = lambda t: fake_margin(L, mp.e**t)
    lo, hi = mp.log(mp.pi * L**4), mp.log(10 * mp.pi * L**4)
    for _ in range(200):
        mid = (lo + hi) / 2
        if f(mid) < 0: lo = mid
        else: hi = mid
    return mp.e**lo

def T_suf(L, width=0.25):
    """least T (geometric bisection) for which the CH-based criterion holds on [L - width, L]."""
    def ok(T):
        f, lmin = CH.ch_prop91(T)
        B = CH.FnBounds([(CH.ch_cor13, 0.0), (f, lmin)], kmax=int(L) + 3)
        return CH.ok_adaptive(B, mpf(L) - width, mpf(L))
    lo, hi = 1e10, 1e17
    assert ok(hi)
    for _ in range(40):
        mid = (lo * hi) ** 0.5
        if ok(mid): hi = mid
        else: lo = mid
    return hi

def L_bar(T, lo=500, hi=4000):
    """largest L (bisection) with T_bar(L) < T: beyond it, A_T contains counterexamples."""
    lo, hi, T = mpf(lo), mpf(hi), mpf(T)
    for _ in range(60):
        mid = (lo + hi) / 2
        if T_bar(mid) < T: lo = mid
        else: hi = mid
    return lo

if __name__ == "__main__":
    T_RH = mpf(3000175332800)
    print(f"[E1] T = 3,000,175,332,800: A_T has counterexamples at every L >= {mp.nstr(L_bar(T_RH), 8)} "
          f"(Theorem A reaches 753.985)")
    tb = T_bar(mpf('3157.93'))
    print(f"[E2] T_bar(3157.93) = {mp.nstr(tb, 10)}; Theorem D uses T* = 9.33566e14, ratio "
          f"{mp.nstr(mpf('9.33566e14') / tb, 8)}")
    print("   L        T_bar(L)          T_suf(L)        T_suf/T_bar    3 pi L^4/(1+14/L)")
    for L in [800, 1000, 1500, 2000, 2500, 3000, 3157.93]:
        tb = T_bar(L)
        ts = T_suf(L)
        print(f"{L:9.2f}  {mp.nstr(tb, 8):>16}  {ts:16.8g}  {float(ts / tb):12.6f}   "
              f"{mp.nstr(3 * mp.pi * mpf(L)**4 / (1 + 14 / mpf(L)), 8)}", flush=True)
