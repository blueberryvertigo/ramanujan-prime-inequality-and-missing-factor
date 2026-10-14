"""Johnston-Leong-Tudzi (arXiv:2408.04143v3), Appendix A: Propositions A.6/A.7 (El Marraki-Schoenfeld iteration).
Reproduce Table 1 rows that use beta = 4 with B = B4 = 75.1090 (Lemma A.5; = eps(2000,2500)*2500^4 from
BKLNW Table 8 row b = 2000), then redo them with B4 corrected (eps(2000) re-optimised with log(3c), or with
Bhattacharjee's E2).  (c0, x0) = (1/4345, 2160535) as in the paper; lambda1 = log xM of the paper."""
import math
PI2 = math.pi ** 2
def prop(alpha, A, beta, B, lam1, c0=1 / 4345, x0=2160535.0):
    D = (3 * B * beta / (PI2 * A)) ** (1 / (beta + 1))
    f = lambda q: 2 * A * D / (1 - q) ** alpha + 6 * B / (PI2 * D ** beta) * (1 - q)
    q0 = D * lam1 ** ((alpha - beta) / (beta + 1))
    cN = max(f(0), f(q0)) if q0 < 1 else f(0)
    a = (alpha * beta - 1) / (beta + 1); k = a + 1
    lx0 = math.log(x0)
    cM = (cN + (c0 + cN / lx0 ** k) * x0 * lam1 ** k / math.exp(lam1)
          + cN * lam1 ** k / (lx0 ** (k + 1) * math.exp(lam1 / 10)) + cN / ((0.9) ** (k + 1) * lam1))
    lam_min = (0.8 * B / (A * D ** beta * math.log(D))) ** ((beta + 1) / (beta - alpha))
    return k, cN, cM, lam_min
B4 = 7.51090e1
for label, Bnew in [("printed B4", B4), ("B4 x 1.0205 (log3c)", B4 * 1.96220 / 1.92279), ("B4 x 1.0226 (thesis E2)", B4 * 1.96625 / 1.92279)]:
    k1, cN1, cM1, lm1 = prop(0, 1 / 4345, 4, Bnew, 43)
    k2, cN2, cM2, lm2 = prop(4 / 5, cM1, 3, 3.95239e-1, 93)
    A5 = cM2 * 385 ** (2 - 27 / 20)
    print(f"{label:24s} B4={Bnew:.4f}: row1 k={k1:.3f} cM={cM1:.6e} (paper 7.80973e-3; lambda_min {lm1:.2f}) | "
          f"row2 k={k2:.4f} cM={cM2:.6e} (paper 6.09073e-2) -> (A.5) const {A5:.5f} (paper 2.91890)")
# last row: alpha = 31/16, A = 4.30429e-1, beta = 4, B = B4, log xM = 385
for label, Bnew in [("printed B4", B4), ("B4 x 1.0205 (log3c)", B4 * 1.96220 / 1.92279), ("B4 x 1.0226 (thesis E2)", B4 * 1.96625 / 1.92279)]:
    k7, cN7, cM7, lm7 = prop(31 / 16, 4.30429e-1, 4, Bnew, 385)
    print(f"{label:24s} row7 k={k7:.4f} cN={cN7:.6e} cM={cM7:.6e} (paper 4.88346; lambda_min {lm7:.1f})")
