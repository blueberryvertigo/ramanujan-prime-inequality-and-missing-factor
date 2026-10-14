"""Remark (Buthe's method, refereed inputs) at T* = 1.3277e15 with Buthe's Theorem 1 as corrected in Bhattacharjee's thesis:
if RH holds for 0 < gamma <= T*, then (C) holds on [e^{671.85}, e^{3035.502}]; with Proposition A' and Theorem C this
gives (R) for every x >= 38,358,837,683."""
import time
import buthe_bound as BB
BB.E2MODE = 'thesis'
import buthe_ramanujan as BR
T_STAR = 1.3277e15
t0 = time.time()
rows = BR.build_rows()
rT = BR.build_rows(T_top=T_STAR, bmin=600)
B = BR.Bounds(BR.merge(rows, rT))
for b in (600, 1000, 2000, 3000):
    d = dict(rT)[b]
    print(f"row b={b}: delta <= {float(d.b):.5g} under RH to {T_STAR:.5g}")
end, fail = BR.run(B, 671.85, 3035.502, 0.5)
print(f"[D] RH to {T_STAR:.5g}: (C) verified contiguously for L in [671.85, {end}]; failure: {fail}")
print(f"total {time.time() - t0:.0f} s")
