"""Robustness of Proposition A' on [10^19, e^247.5]: inflate every row of Buthe's Theorem 1 (thesis-corrected)
by a factor K and report how far (C) still holds.  argv: K"""
import sys, time
import buthe_bound as BB
BB.E2MODE = 'thesis'
K = float(sys.argv[1])
_d = BB.delta
BB.delta = lambda *a, **k: _d(*a, **k) * K
import buthe_ramanujan as BR
t0 = time.time()
rows = BR.build_rows()
B = BR.Bounds(rows)
end, fail = BR.run(B, BR.B19, 700, 0.25)
print(f"[A' x{K}] (C) verified contiguously for L in [{BR.B19}, {end}]; first failure {fail}   ({time.time()-t0:.0f} s)", flush=True)
