"""
Trade-off between the RH verification height T and the zero-free-region constant R.
Hypothesis H(T, R): every zero with 0 < gamma <= T has beta = 1/2, and zeta(beta + i gamma) != 0 for |gamma| > T,
beta >= 1 - 1/(R log|gamma|).  Then on [L_A - 1, L_TOP]
  |psi(x) - x| <= x * min( [CH, Cor 1.3],  [CH, Prop 9.1](T),  FKS Thm 1.1 with H0 = T and constant R ),
and Theorems A (L <= L_A) and C (L >= L_TOP = 3035.002) cover the rest.  For each R, bisection (log scale) for the
least T in [T_PT, 6e14] for which (C) holds on [L_A - 1, L_TOP + 0.5].  T_PT = 3,000,175,332,800 (Platt--Trudgian).
FKS rows below log x0 = 1000: as in zfr_B.py (ingredients need x > e^50, 3 log x < T' < sqrt(x)/3, T' > t0; t0 is
chosen via Lmax so that each row is valid on the whole L-range where it is used).
argv: R tag [density mode, default KLN618]
"""
import sys, time
import fks_eps as F
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
R = float(sys.argv[1]); tag = sys.argv[2]
F.DENSITY['mode'] = sys.argv[3] if len(sys.argv) > 3 else 'KLN618'
LA, L_TOP = 753.985, 3035.002
T_PT = 3000175332800.0
GRID = list(range(700, 1200, 3)) + list(range(1200, int(L_TOP) + 5, 4))

def rows_fn(T):
    rows = []
    for i, L in enumerate(GRID):
        Lnext = GRID[i + 1] if i + 1 < len(GRID) else L + 4
        try:
            e, s2 = F.best_eps_table7(L, N=4000, K=40000, R=R, H0v=T, Lmax=Lnext)
        except (ValueError, AssertionError):
            continue
        rows.append((L, Lnext, I(repr(float(e)))))
    def f(l):
        cand = [v for L, Ln, v in rows if L <= l.a and l.b <= Ln]
        return cand[-1] if cand else I(1)
    return f, float(GRID[0])

def closes(T):
    g, gmin = CH.ch_prop91(T)
    f, fmin = rows_fn(T)
    B = CH.FnBounds([(CH.ch_cor13, 0.0), (g, gmin), (f, fmin)])
    e, fl = CH.run(B, LA - 1, L_TOP + 0.5, 0.5)
    return fl is None, fl

if __name__ == "__main__":
    t0 = time.time()
    ok, fl = closes(T_PT)
    print(f"   [tradeoff R={R}] T = {T_PT:.6g} (current verification): {'closes' if ok else 'fails at ' + str(fl)}", flush=True)
    if ok:
        print(f"[tradeoff R={R}, {tag}] current verification height suffices   ({time.time() - t0:.0f} s)")
        sys.exit()
    lo_T, hi_T = T_PT, 6e14
    ok, fl = closes(hi_T)
    assert ok, f"T = 6e14 fails at {fl}"
    for _ in range(12):
        mid = (lo_T * hi_T) ** 0.5
        ok, fl = closes(mid)
        print(f"   [tradeoff R={R}] T = {mid:.6g}: {'closes' if ok else 'fails at ' + str(fl)}   ({time.time()-t0:.0f} s)", flush=True)
        if ok: hi_T = mid
        else: lo_T = mid
    print(f"[tradeoff R={R}, {tag}] RH to height {hi_T:.6g} suffices; fails at {lo_T:.6g}   ({time.time() - t0:.0f} s)")
