"""
Theorem A' (unconditional) and Theorem D (conditional on a finite RH verification), using
psi bounds that do not depend on BKLNW Table 8:

  * x <= 10^19:  |psi(x) - x| <= 0.94 sqrt(x)      (Buthe, Math. Comp. 87 (2018), Thm 2, (1.5));
  * x >= 10^19:  the least of
      - Buthe 2016 Theorem 1 with E2 corrected (buthe_bound.py), RH to height
        T_RH = 3,000,175,332,800 (Platt-Trudgian 2021), parameters optimised per row;
      - FKS Table 5 rows and FKS Cor. 1.4 (as in extend_ramanujan.py).

Criterion (C), the exact main term and the interval arithmetic are those of verify_ramanujan.py.
"""
import sys, time
from mpmath import iv, mp, mpf
import verify_ramanujan as V
import extend_ramanujan as X
import buthe_bound as BB

I = V.I
T_RH = 3000175332800
B19 = 43.749            # < log(10^19) = 43.7491...

def best_row(b, T_list, c_list, a_list):
    best = None
    for T in T_list:
        for c in c_list:
            for a in a_list:
                try:
                    d = BB.delta(b, c, T, a, log3c=True).b
                except AssertionError:
                    continue
                if best is None or d < best[0]: best = (d, c, T, a)
    return best

def build_rows(T_top=T_RH, bmin=B19, bmax=3300, verbose=False):
    """rows (b, delta): |psi(x)-x| <= delta x for x >= e^{b + 1e-6}; RH up to T_top assumed."""
    small_T = [T_top] + [t for t in [1e12, 3e11, 1e11, 3e10, 1e10, 3e9, 1e9, 3e8] if t < T_top]
    bs = [B19] + list(range(44, 121)) + list(range(125, 801, 5)) + list(range(825, bmax + 1, 25))
    bs = [b for b in bs if b >= bmin]
    rows = []
    for b in bs:
        if b < 70:
            r = best_row(b, small_T, [28, 30, 32, 33, 34, 35], [1e-3, 3e-3, 1e-2])
        else:
            import math
            c0 = 35 + max(0.0, math.log(T_top / T_RH))         # optimum c ~ log T + 6.3
            if T_top == T_RH:
                cs, als = [c0 + d for d in (-1.5, -0.75, 0, 0.75, 1.5)], [3e-4, 1e-3, 3e-3]
            else:
                cs, als = [c0 + d for d in (-1, 0, 1)], [3e-4, 1e-3]
            r = best_row(b, [T_top], cs, als)
        rows.append((b, I(r[0])))
        if verbose and b in (bmin, 50, 60, 100, 300, 600, 700, 1000, 3000):
            print(f"   row b={b}: delta <= {mp.nstr(r[0], 6)}  (c={r[1]}, T={r[2]:.4g}, alpha={r[3]})")
    return rows

class Bounds:
    """unit-segment and pointwise psi bounds from Buthe 2018 / Buthe 2016 rows / FKS."""
    def __init__(self, rows, kmax=3400):
        self.rows = rows
        self.EPS = {}
        for k in range(20, kmax + 1):
            if k + 1 <= B19:
                e = I(0.94) * iv.exp(I(-k) / 2)
            elif k < B19:                      # [k, k+1] straddles 10^19: Buthe 2018 below, row B19 above
                e = I(max((I(0.94) * iv.exp(I(-k) / 2)).b, rows[0][1].b))
                assert rows[0][0] == B19
            else:
                cands = [X.classical(I(k)).b]
                rb = [d for b, d in rows if b + 1e-6 <= k]
                if rb: cands.append(rb[-1].b)
                fr = [e for b, e in X.FKS5 if b <= k]
                if fr: cands.append(fr[-1].b)
                e = I(min(cands))
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
        rb = [d for b, d in self.rows if b + 1e-6 <= l.a]
        if rb: cands.append(rb[-1].b)
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

def run(B, lo, hi, step, stop=True):
    """returns the right end of the last verified interval (contiguous from lo)."""
    last = lo
    for a, b in V.grid(lo, hi, step, 'B'):
        if not ok_adaptive(B, mpf(a), mpf(b)):
            return last, (a, b)
        last = b
    return last, None

def merge(rows_a, rows_b):
    """pointwise minimum of two row lists on the same b grid (rows_b may start later)."""
    mb = {b: d for b, d in rows_b}
    return [(b, I(min(d.b, mb[b].b)) if b in mb else d) for b, d in rows_a]

def eps_need(L):
    L = mpf(L)
    return (1 + 14 / L) / (3 * L**4)

D_TOP = float(__import__("os").environ.get("D_TOP", "3158.5"))   # top of the range closed by Theorem D

if __name__ == "__main__":
    t0 = time.time()
    print("Theorem A': rows at T_RH = 3,000,175,332,800")
    rows = build_rows(verbose=True)
    B = Bounds(rows)
    print(f"   built {len(rows)} rows in {time.time() - t0:.0f} s")
    end, fail = run(B, B19, 800, 0.25)
    print(f"[A'] coarse: verified contiguously for L in [{B19}, {end}]; first failure {fail}")
    end2, fail2 = run(B, float(end), float(end) + 2, 0.002)
    print(f"[A'] fine:   verified contiguously for L in [{B19}, {end2}]; first failure {fail2}")
    LA = float(end2)
    print("[gap] L      best eps (T_RH)   needed eps     shortfall")
    for L in [674, 700, 800, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3100, 3150]:
        best = float(mpf(B.point(I(L)).b))
        print(f"      {L:5d}  {best:14.4g}  {float(eps_need(L)):12.4g}  {best / float(eps_need(L)):8.3g}x")
    # Theorem D: smallest verification height T* (on a geometric bisection) closing [L_A, D_TOP]
    print(f"Theorem D: RH to height T* => (R) on [e^L_A, e^{D_TOP}]")
    def closes(T):
        rT = build_rows(T_top=T, bmin=600)
        BT = Bounds(merge(rows, rT))
        e, f = run(BT, LA - 1, D_TOP, 0.5)
        return f is None, e, f
    lo_T, hi_T = float(T_RH), 1e17
    ok, e, f = closes(hi_T)
    print(f"   T = {hi_T:.4g}: closes = {ok}", flush=True)
    assert ok
    for _ in range(16):
        mid = (lo_T * hi_T) ** 0.5
        ok, e, f = closes(mid)
        print(f"   T = {mid:.5g}: closes = {ok}" + ("" if ok else f" (first failure at L in {f})"), flush=True)
        if ok: hi_T = mid
        else: lo_T = mid
    print(f"[D] RH verified up to height {hi_T:.5g} suffices for (R) on [e^{LA - 1:.2f}, e^{D_TOP}]; "
          f"the method does not reach it at {lo_T:.5g}")
    print(f"total {time.time() - t0:.0f} s")
