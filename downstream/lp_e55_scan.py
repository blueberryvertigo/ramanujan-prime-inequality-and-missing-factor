"""Broad grid check of the least corrected Buthe Thm 1 value at x >= e^55 with RH height 3.061e10
(Buthe 2016 Table 1 setting; E3 by zeros_model): c in [20, 40], alpha in [0.01, 0.4], T in [1e9, 3.061e10]."""
import sys
sys.path.insert(0, '.')
from zeros_model import delta
H = 3.061e10
best = None
for i in range(41):
    c = 20 + 0.5 * i
    for j in range(40):
        a = 0.01 * (j + 1)
        for k in range(16):
            T = min(H, 1e9 * 10 ** (0.1 * k))
            r = delta(55, c, T, a, True)
            if r and (best is None or r[0] < best[0]): best = (r[0], c, a, T)
print(f"coarse best with log3c: {best[0]:.6e} at c={best[1]}, alpha={best[2]}, T={best[3]:.4e}")
b0 = best
for i in range(41):
    c = b0[1] - 0.5 + 0.025 * i
    for j in range(41):
        a = max(0.001, b0[2] - 0.01 + 0.0005 * j)
        for T in [H, H * 0.98, H * 0.95, H * 0.9]:
            r = delta(55, c, T, a, True)
            if r and r[0] < best[0]: best = (r[0], c, a, T)
print(f"fine best with log3c: {best[0]:.6e} at c={best[1]:.3f}, alpha={best[2]:.4f}, T={best[3]:.4e}")
print(f"theta bound: {best[0] + 1.02 * 2.718281828459045 ** (-27.5) + 3 * 2.718281828459045 ** (-110 / 3):.6e} vs 2.9e-10")
