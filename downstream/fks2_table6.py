"""FKS2 Table 6 (Corollary 23) and Corollary 22: E_pi(x) <= A (L/R)^B exp(-C sqrt(L/R)), L = log x >= log x0,
R = 5.5666305.  On the affected range 44 <= L < 2073 compare the bound with eps_pi,num of FKS2 Table 4
(interpolated): the minimum ratio bound/eps is the room available there."""
import math, sys
sys.path.insert(0, '.')
from fks2_table7 import eps_pi
R = 5.5666305
ROWS = [(0.000120, 0.25, 1.00, 22.955), (0.826, 0.25, 1.00, 1.0), (1.41, 0.5, 1.5, 2.0), (1.76, 1.0, 1.5, 3.0),
        (2.22, 1.5, 1.5, 3.0), (12.4, 1.5, 1.9, 1.0), (38.8, 1.5, 1.95, 1.0), (121.107, 1.5, 2.0, 1.0), (6.60, 2.0, 2.0, 3.0)]
for A, B, C, lx0 in ROWS:
    best = None
    L = 44.0
    while L <= 2072.0:
        r = A * (L / R)**B * math.exp(-C * math.sqrt(L / R)) / eps_pi(L)
        if best is None or r < best[0]: best = (r, L)
        L += 0.25
    print(f"A={A}, B={B}, C={C}, log x0={lx0}: min over 44<=L<2073 of bound/eps_pi = {best[0]:.3f} at L={best[1]:.2f}")
