"""How large is the log(3c) correction once Buthe's parameters (c, alpha) are re-optimised?
For each printed row: least corrected value over c in [20, 48] (step 0.25) and alpha on a log grid,
at the row's own height T; and whether any (c, alpha) reaches the printed value.
Then: the reach of criterion (C) with corrected rows at Buthe's/BKLNW's height 2,445,999,556,030."""
import math
from mpmath import mp
import buthe_bound as BB
ROWS = [("Buthe Table 2, e^1000", 1000, 2.445e12, '1.94751e-12'), ("Buthe Table 2, e^4000", 4000, 2.445e12, '1.90866e-12'),
        ("BKLNW Table 8, b=200", 200, 2.445999556030e12, '2.18154e-12'), ("BKLNW Table 8, b=600", 600, 2.445999556030e12, '1.98894e-12'),
        ("Buthe Table 1, e^500", 500, 3.061e10, '1.47067e-10')]
CS = [20 + 0.25 * i for i in range(113)]
ALPHAS = [10 ** (-5 + 0.25 * j) for j in range(19)]          # 1e-5 ... 0.3
for name, b, T, printed in ROWS:
    best = None
    for c in CS:
        for a in ALPHAS:
            try:
                d = BB.delta(b, c, T, a, log3c=True).b
            except AssertionError:
                continue
            if best is None or d < best[0]: best = (d, c, a)
    p = float(printed)
    print(f"{name}: printed {printed}; least corrected value {mp.nstr(best[0], 6)} at c={best[1]}, alpha={best[2]:.3g} "
          f"(+{100 * (float(best[0]) / p - 1):.2f}%); printed value reachable: {float(best[0]) <= p}", flush=True)
import buthe_ramanujan as BR
H = 2445999556030
rows = BR.build_rows(T_top=H)
B = BR.Bounds(rows)
end, fail = BR.run(B, BR.B19, 800, 0.25)
end2, fail2 = BR.run(B, float(end), float(end) + 2, 0.002)
print(f"Reach of (C) with corrected Buthe rows at height {H}: verified for L in [{BR.B19}, {end2}]; first failure {fail2}")
