"""fks_fine.py with a different zero-free-region constant R (argv[1]) and row range (argv[2], argv[3]):
  R = 5.558691  Mossinghoff-Trudgian-Yang, Res. Number Theory (2023): no zeros for sigma >= 1 - 1/(R log|t|), |t| >= 2;
  R = 4.896     Bellotti-Trudgian-Yang, arXiv:2603.21490 (2026): no zeros for t >= 3, sigma >= 1 - 1/(R log t).
FKS's Theorem 1.1 uses the zero-free region only through R (eps3: H_sigma; eps4: t0 and x^{-1/(R log t)})."""
import json, sys, time
from mpmath import mpf
import fks_eps as F
import ch_ramanujan as CH
import verify_ramanujan as V
I = V.I
R = float(sys.argv[1]); LO, HI = int(sys.argv[2]), int(sys.argv[3]); tag = sys.argv[4]
if len(sys.argv) > 5: F.DENSITY['mode'] = sys.argv[5]   # KLN | KLN618 | BEL618
t0 = time.time()
import os
if os.path.exists(f"fks_rows_{tag}.json"):
    rows = [tuple(r) for r in json.load(open(f"fks_rows_{tag}.json"))]
else:
    rows = []
    for L in range(LO, HI + 1):
        e, s2 = F.best_eps_table7(L, N=10000, K=100000, R=R, Lmax=L + 1)
        rows.append((L, float(e), s2))
    json.dump(rows, open(f"fks_rows_{tag}.json", "w"))
for L, e, s2 in rows[::25]:
    print(f"  row log x0 = {L}: eps <= {e:.5e} (sigma_2 = {s2})")
print(f"  R = {R}: {len(rows)} rows in {time.time() - t0:.0f} s", flush=True)
vals = [(L, I(repr(float(e)))) for L, e, s2 in rows]
def f(l):
    if l.b > HI + 1: return I(1)        # rows are certified only up to their own L + 1
    cand = [v for L, v in vals if L <= l.a]
    return cand[-1] if cand else I(1)
B = CH.FnBounds([(f, float(LO))], kmax=10050)
last_fail = None
for a, b in V.grid(LO + 20, HI + 40, 0.5, 'B'):
    if not CH.ok_adaptive(B, mpf(a), mpf(b)): last_fail = b
lo = float(last_fail) - 1
LC = lo
for a, b in V.grid(lo, float(last_fail) + 0.5, 0.002, 'B'):
    if not CH.ok_adaptive(B, mpf(a), mpf(b)): LC = b
end, fail = CH.run(B, float(LC), 10000, 0.5)
print(f"[C, R={R}] (R) verified contiguously for L in [{LC}, {end}] (failure: {fail})   ({time.time() - t0:.0f} s)", flush=True)
LA = 753.985
def closes(T):
    g, gmin = CH.ch_prop91(T)
    BT = CH.FnBounds([(CH.ch_cor13, 0.0), (g, gmin), (f, float(LO))])
    e, fl = CH.run(BT, LA - 1, float(LC) + 0.5, 0.5)
    return fl is None, fl
lo_T, hi_T = 3000175332800.0, 1e16
for _ in range(18):
    mid = (lo_T * hi_T) ** 0.5
    ok, fl = closes(mid)
    print(f"   [D, R={R}] T = {mid:.6g}: {'closes' if ok else 'fails at ' + str(fl)}", flush=True)
    if ok: hi_T = mid
    else: lo_T = mid
print(f"[D, R={R}] RH to height {hi_T:.6g} suffices; fails at {lo_T:.6g}   ({time.time() - t0:.0f} s)")
