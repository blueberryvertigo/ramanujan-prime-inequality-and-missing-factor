"""
Ramanujan's inequality from the Chirre-Helfgott bounds (arXiv:2512.15709v1).

  [CH, Cor. 1.3]  for every x >= 1:  |psi(x) - x| <= pi/(3e12) x + 113.67 sqrt(x)
                  (uses only the verification of RH to height 3,000,175,332,800, Platt-Trudgian).
  [CH, Prop. 9.1 with sigma = 0, and (9.2)]  if RH holds up to height T >= 1e7, then for
                  x > max(T, 1e9):  |psi(x) - x| <= pi/(T - 1 - pi/3) x + C_T sqrt(x),
                  C_T = log(T/2pi)^2/(2pi) - log(T/2pi)/(6pi).
  [FJ, Thm 1.4]   if every zero has beta <= theta with 7/8 <= theta <= 1, then for x > x0
                  E_psi(x) < A(x0) log(x)^2 x^{-1/4} + C(x0) x^{theta-1}  (column [7/8, 1]).

Modes:
  A+  unconditional: CH Cor. 1.3 (and FKS above e^1000 where smaller), from L = 43.749 (Axler 2018 below).
  D+  RH to height T: CH Prop. 9.1; bisection for the least T closing [L_A+ - 1, 3158.5].
  B+  zero-free half-plane Re s > theta: FJ Thm 1.4 + CH Cor. 1.3; bisection for the largest theta.
Theorem C (extend_ramanujan.py) covers L >= 3157.932 unconditionally.
"""
import math, time
from mpmath import iv, mp, mpf
import verify_ramanujan as V
import extend_ramanujan as X

I = V.I
PI = iv.pi
T_RH = 3000175332800
B19 = 43.749

def ch_cor13(l):
    """[CH, Cor 1.3] at log x = l (interval)."""
    return PI / I('3e12') + I('113.67') * iv.exp(-l / 2)

def ch_prop91(T):
    T = I(T)
    CT = iv.log(T / (2 * PI))**2 / (2 * PI) - iv.log(T / (2 * PI)) / (6 * PI)
    lmin = max(math.log(float(T.b)), math.log(1e9)) + 1e-9
    def f(l):
        assert l.a > lmin
        return PI / (T - 1 - PI / 3) + CT * iv.exp(-l / 2)
    return f, lmin

FJ_ROWS = [(5, '0.0198944', '0.02636425'), (150, '0.0198944', '0.02636425'), (200, '0.0198944', '0.02506448'),
           (300, '0.0198944', '0.01840514'), (500, '0.0198944', '0.01505842'), (1000, '0.0198944', '0.01320991'),
           (2000, '0.0198944', '0.01243612'), (5000, '0.0198944', '0.01201123')]

def fj_theta(theta):
    th = I(theta)
    assert 0.875 <= theta < 1
    def f(l):
        rows = [(b, I(a), I(c)) for b, a, c in FJ_ROWS if b < l.a]
        b, A, C = rows[-1]
        return A * l**2 * iv.exp(-l / 4) + C * iv.exp((th - 1) * l)
    return f, 5.0 + 1e-9

class FnBounds:
    """psi-bounds given as a list of (function of l, valid for l > lmin); each function is
    nonincreasing in l on the range used, so its value at the left end of a unit segment is valid
    on the whole segment and beyond.  Below 10^19 Buthe 2018 (0.94 sqrt x) is used."""
    def __init__(self, fns, kmax=3400):
        self.fns = fns
        self.EPS = {}
        for k in range(20, kmax + 1):
            if k + 1 <= B19:
                self.EPS[k] = I(0.94) * iv.exp(I(-k) / 2)
                continue
            cands = [X.classical(I(k)).b]
            fr = [e for b, e in X.FKS5 if b <= k]
            if fr: cands.append(fr[-1].b)
            for f, lmin in fns:
                if k > lmin: cands.append(f(I(k)).b)
            e = I(min(cands))
            if k < B19:
                e = I(max(e.b, (I(0.94) * iv.exp(I(-k) / 2)).b))
            self.EPS[k] = e
        self.C = {20: I(0)}
        einv = iv.exp(I(-1))
        for k in range(20, kmax):
            self.C[k + 1] = self.C[k] * einv + self.EPS[k] / (I(k + 1)**2 * (1 - I(2) / k))
        self.kmax = kmax

    def point(self, l):
        if l.b <= B19:
            return I(0.94) * iv.exp(-l.a / 2)
        N = int(mpf(l.a))
        cands = [self.EPS[N].b, X.classical(l).b]
        for f, lmin in self.fns:
            if l.a > lmin: cands.append(f(l).b)
        return I(min(cands))

    def d(self, l):
        N = int(mpf(l.a))
        assert 21 <= N < self.kmax
        pt, it = V.psitheta_terms(l)
        point = (self.point(l) + pt) / l
        integ = V.K0 * iv.exp(-l) + self.C[N] * iv.exp(N - l) + self.EPS[N] / (l**2 * (1 - I(2) / N))
        return point + integ + it

def ok_adaptive(B, a, b, depth=0):
    if V.margin(I([a, b]), B.d)[0].a > 0: return True
    if depth >= 12: return False
    m = (a + b) / 2
    return ok_adaptive(B, a, m, depth + 1) and ok_adaptive(B, m, b, depth + 1)

def run(B, lo, hi, step):
    last = lo
    for a, b in V.grid(lo, hi, step, 'B'):
        if not ok_adaptive(B, mpf(a), mpf(b)):
            return last, (a, b)
        last = b
    return last, None

if __name__ == "__main__":
    t0 = time.time()
    # ---- A+ ----
    BA = FnBounds([(ch_cor13, 0.0)])
    end, fail = run(BA, B19, 900, 0.25)
    end2, fail2 = run(BA, float(end), float(end) + 2, 0.002)
    LA = float(end2)
    print(f"[A+] (R) verified for L in [{B19}, {end2}]; first failure {fail2}   ({time.time()-t0:.0f} s)", flush=True)
    # ---- D+ ----
    def closes_D(T):
        f, lmin = ch_prop91(T)
        B = FnBounds([(ch_cor13, 0.0), (f, lmin)])
        e, fl = run(B, LA - 1, 3158.5, 0.5)
        return fl is None, fl
    lo_T, hi_T = float(T_RH), 1e16
    assert closes_D(hi_T)[0]
    for _ in range(18):
        mid = (lo_T * hi_T) ** 0.5
        ok, fl = closes_D(mid)
        print(f"   [D+] T = {mid:.6g}: {'closes' if ok else 'fails at ' + str(fl)}", flush=True)
        if ok: hi_T = mid
        else: lo_T = mid
    print(f"[D+] RH to height {hi_T:.6g} suffices; method fails at {lo_T:.6g}   ({time.time()-t0:.0f} s)", flush=True)
    # ---- B+ ----
    def closes_B(theta):
        f, lmin = fj_theta(theta)
        B = FnBounds([(ch_cor13, 0.0), (f, lmin)])
        e, fl = run(B, LA - 1, 3158.5, 0.5)
        return fl is None, fl
    lo_t, hi_t = 0.875, 0.999
    assert closes_B(lo_t)[0]
    for _ in range(16):
        mid = (lo_t + hi_t) / 2
        ok, fl = closes_B(mid)
        print(f"   [B+] theta = {mid:.6f}: {'closes' if ok else 'fails at ' + str(fl)}", flush=True)
        if ok: lo_t = mid
        else: hi_t = mid
    print(f"[B+] zero-free half-plane Re s > {lo_t:.6f} suffices; fails at {hi_t:.6f}   ({time.time()-t0:.0f} s)")
