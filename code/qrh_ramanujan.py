"""
Theorem F: if zeta(s) != 0 for Re s > 7/8 (formalized in Lean 4 as OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re,
openai/math), then (R) holds for all x >= 38,358,837,683.

psi-bound under the hypothesis.  For x >= e^40 and max{51, log x} < T < (sqrt x - 2)/5 there is T* in [T, 2.45 T] with
    |psi(x) - x + sum_{|gamma| <= T*} x^rho / rho| <= 2.091 x log x / T          [Cully-Hugill--Johnston, Thm 1.2;
                                                                                  FKS (JMAA 2023) Thm 3.2]
We use 100 * 2 (log x)^2 / T instead (>= the above for log x >= 2), a safety factor far beyond any known correction.
Zeros: beta = 1/2 for 0 < gamma <= H = 3,000,175,332,800 (Platt--Trudgian); beta <= 7/8 for all zeros (hypothesis);
|x^{rho-1}/rho| <= x^{beta-1}/|gamma|; zeros symmetric under gamma -> -gamma; and for t >= 5000,
    sum_{0 < gamma <= t} 1/gamma <= log^2(t/(2 pi)) / (4 pi)                       [Buthe 2016, Lemma 3].
Hence, with T = e^{L/2 - 2} (L = log x >= 110; then T > max{51, L} and T < (sqrt x - 2)/5):
    |psi(x) - x|/x <= 200 L^2/T + 2 x^{-1/2} log^2(H/2pi)/(4pi) + 2 x^{-1/8} log^2(2.45 T/2pi)/(4pi) =: eps_F(L).
Each term is decreasing in L for L >= 110, so eps_F(a) bounds |psi(x)-x|/x for all x >= e^a.
(C) is then checked on [L1, 3035.5] with the machinery of the paper (ch_ramanujan.FnBounds / run), using only eps_F
and the default classical/FKS-Table-5 bounds, i.e. no Chirre--Helfgott input.  Lemma base (x <= 10^19),
Proposition A' (refereed; up to e^672.84 with Buthe's Theorem 1 as corrected) and Theorem C (x >= e^3035.01) cover the rest.
"""
import sys, time
from mpmath import iv, mpf
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
H = I('3000175332800')
PI = iv.pi
LOGH2 = iv.log(H / (2 * PI))**2

def eps_F(L):
    L = I(L)
    lT = L / 2 - 2
    t1 = 200 * L**2 * iv.exp(-lT)
    t2 = 2 * iv.exp(-L / 2) * LOGH2 / (4 * PI)
    t3 = 2 * iv.exp(-L / 8) * (lT + iv.log(I('2.45') / (2 * PI)))**2 / (4 * PI)
    return t1 + t2 + t3

def f(l):                       # bound valid for all x >= e^{l.a}
    return I(eps_F(I(l.a)).b)

if __name__ == "__main__":
    t0 = time.time()
    for L in [110, 200, 250, 300, 400, 500, 672, 1000, 3000]:
        print(f"   eps_F({L}) <= {float(eps_F(L).b):.4e}", flush=True)
    LMIN = 110.0
    B = CH.FnBounds([(f, LMIN)], kmax=3100)
    lo = float(sys.argv[1]) if len(sys.argv) > 1 else 200.0
    # find the first L from which (C) holds contiguously up to 3035.5
    last_fail = None
    for a, b in V.grid(lo, 3036.0, 0.5, 'B'):
        if not CH.ok_adaptive(B, mpf(a), mpf(b)):
            last_fail = b
    start = float(last_fail) if last_fail is not None else lo
    end, fail = CH.run(B, start, 3036.0, 0.5)
    print(f"[F] under Re s > 7/8 zero-free: (R) verified contiguously for L in [{start}, {end}] (failure: {fail})   "
          f"({time.time() - t0:.0f} s)", flush=True)
