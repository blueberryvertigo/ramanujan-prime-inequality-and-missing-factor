"""BKLNW Table 15 / Corollary 11.1 (M_k for k = 1..5 on the middle range).  Table 10 gives
B_k(b, b') ~ eps(b) b'^k on a 25-step grid; near b = 2300 the max is attained:
  M1 = 4.4627e-9 = B1(2300,2325), M2 = 1.0376e-5 = B2(2300,2325)   [eps(2300) from Buthe Thm 1]
  M3 = 2.4334e-2 = B3(2325,2350), M4 = 57.184 = B4(2325,2350), M5 = 1.3441e5 = B5(2350,2375) [Platt-Trudgian method]
(1) reproduce eps(2300) with Buthe Thm 1 as printed (H = 2,445,999,556,030), (2) re-optimise with log(3c),
(3) recompute B_k(2300,2325) and compare with the printed maxima."""
import sys, math
sys.path.insert(0, '.')
from zeros_model import delta, optimise
H = 2445999556030.0
eps_printed = 4.4627e-9 / 2325          # = B1(2300,2325)/2325
r0 = optimise(2300, False, H, c0=33.3, T0=H, a0=0.003)
r1 = optimise(2300, True, H, c0=34.8, T0=H, a0=0.003)
fixed = delta(2300, r0[1], H, r0[3], True)[0]
print(f"eps(2300) from Table 10: {eps_printed:.5e};  Thm 1 as printed, optimised: {r0[0]:.5e} (c={r0[1]:.3f}, alpha={r0[3]:.4g})")
print(f"with log(3c) at those parameters: {fixed:.5e} (+{fixed/r0[0]-1:.2%});  re-optimised with log(3c): {r1[0]:.5e} (c={r1[1]:.3f}, alpha={r1[3]:.4g}) (+{r1[0]/r0[0]-1:.2%})")
printed_max = {1: 4.4627e-9, 2: 1.0376e-5, 3: 2.4334e-2, 4: 57.184, 5: 1.3441e5}
B2300 = {1: 4.4627e-9, 2: 1.0376e-5, 3: 2.4124e-2, 4: 56.088, 5: 1.3040e5}   # Table 10 row b = 2300
for k in range(1, 6):
    for name, e in [("re-opt", r1[0]), ("fixed", fixed), ("+2.6%", 1.026 * r0[0]), ("+8%", 1.08 * r0[0])]:
        newB = B2300[k] * e / eps_printed
        flag = "EXCEEDS printed M_k" if newB > printed_max[k] else "ok"
        print(f"k={k} {name:6s}: B_k(2300,2325) -> {newB:.5g}  vs printed M_k {printed_max[k]:.5g}  [{flag}]")
# fine-grid estimate: PT-method eps from Table 10 B1 rows (b = 2325..2400), log-linear in b
pts = [(2325, 4.4062e-9 / 2350), (2350, 4.2245e-9 / 2375), (2375, 4.0498e-9 / 2400), (2400, 3.8820e-9 / 2425)]
rate = math.log(pts[0][1] / pts[-1][1]) / (pts[-1][0] - pts[0][0])
ept = lambda b: pts[0][1] * math.exp(rate * (pts[0][0] - b))
print(f"PT-method eps: {', '.join(f'{b}:{e:.4e}' for b, e in pts)}; decay {rate*25:.4f} per 25")
for name, e in [("re-opt", r1[0]), ("fixed", fixed)]:
    bc = pts[0][0] - math.log(e / pts[0][1]) / rate       # crossover eps_PT(bc) = corrected Buthe
    for k in (3, 4, 5):
        m = max(min(e, ept(b)) * (b + 1.0)**k for b in [2200 + 0.5 * i for i in range(600)])
        print(f"  fine grid (step 1), {name}: crossover b_c={bc:.1f}; max_b min(eps_Bu,eps_PT)(b)*(b+1)^{k} = {m:.5g} vs printed {printed_max[k]:.5g}")
