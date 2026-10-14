"""Least (re-optimised) value of Buthe's bound for the rows that matter downstream, under the three E2
versions of e2_variants.py; and the excess over the printed value."""
import sys, math
sys.path.insert(0, '.')
from e2_variants import optimise
sys.stdout.reconfigure(encoding="utf-8")
H1, H2, H3, H4 = 3.061e10, 2.445e12, 2445999556030.0, 3000175332800.0
ROWS = [("Buthe T1 e^55", 55, H1, 2.88434e-10, 30, 0.1), ("Buthe T1 e^500", 500, H1, 1.47067e-10, 30, 0.01),
        ("Buthe T2 e^1000", 1000, H2, 1.94751e-12, 34, 0.005), ("Buthe T2 e^4000", 4000, H2, 1.90866e-12, 34, 0.001),
        ("BKLNW T8 b=100", 100, H3, 2.45299e-12, 34, 0.06), ("BKLNW T8 b=200", 200, H3, 2.18154e-12, 34, 0.03),
        ("BKLNW T8 b=600", 600, H3, 1.98894e-12, 34, 0.01), ("BKLNW T8 b=1000", 1000, H3, 1.94751e-12, 34, 0.005),
        ("BKLNW T8 b=2000", 2000, H3, 1.92279e-12, 34, 0.003), ("BKLNW ext. b=2300", 2300, H3, 1.91944e-12, 34, 0.0025),
        ("JY Lemma 2.3 b=2000", 2000, H4, 1.570e-12, 34, 0.003), ("FKS2 T3 b=1350", 1350, 3.0e12, 1.5795e-12, 34, 0.003)]
print(f"{'row':22s} {'printed':>11s} | {'opt printed-E2':>14s} | {'opt log3c':>11s} {'excess':>7s} | {'opt thesis':>11s} {'excess':>7s}")
for name, b, H, p, c0, a0 in ROWS:
    r0 = optimise(b, 'printed', H, c0=c0 - 1, a0=a0)
    r1 = optimise(b, 'log3c', H, c0=c0 + 1, a0=a0)
    r2 = optimise(b, 'thesis', H, c0=c0 + 1, a0=a0)
    print(f"{name:22s} {p:11.5e} | {r0[0]:14.5e} | {r1[0]:11.5e} {r1[0]/p-1:+7.2%} | {r2[0]:11.5e} {r2[0]/p-1:+7.2%}   (c: {r0[1]:.2f}/{r1[1]:.2f}/{r2[1]:.2f})", flush=True)
