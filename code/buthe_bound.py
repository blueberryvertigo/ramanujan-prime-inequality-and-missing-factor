"""
Rigorous evaluation of Buthe's Theorem 1 (Math. Comp. 85 (2016) 2483-2498; arXiv:1410.7015v4):

  if beta = 1/2 for every zero with 0 < gamma <= c/eps, then for x >= e^{alpha eps} x0
      |psi(x) - x| <= x e^{alpha eps} (E1 + E2 + E3),
  E1 = e^{2eps} log(e^eps x0) [ 2 eps |nu_c(alpha)| / log B0 + 2.01 eps / sqrt(x0)
                                + log log(2 x0^2) / (2 x0) ] + (e^{alpha eps} - 1),
  B0 = eps x0 e^{-eps} |nu_c(alpha)| / (2 mu_c(alpha))       (must exceed 1),
  E2 = 0.16 (1 + 1/x0) / sinh(c) * e^{0.71 sqrt(c eps)} * LOG3C * log(c/eps),
  E3 = (2/sqrt(x0)) sum_{0<gamma<=c/eps} ell(gamma)/gamma + 2/x0.

LOG3C.  Theorem 1 as printed has no log(3c) in E2.  But E2 is (1/x) times the bound (4.1) of the
paper's own Proposition 3, which carries the factor log(3c); so does the lemma it rests on (the
paper's Lemma 2, from Franke-Kleinjung-Buthe-Jost, Lemma 2.4).  `table2_check` reproduces the
paper's Table 2 with and without the factor: the printed values match the version WITHOUT it.  We
use the version with it.

E3 is bounded with ell <= 1 and the paper's Lemma 3:  sum_{0<gamma<t} 1/gamma <= log(t/2pi)^2/(4pi)
for t >= 5000.

mu_c, nu_c (paper, Section 5):  eta(t) = c I0(c sqrt(1-t^2)) / (2 sinh c) on [-1,1];
  mu_c(a) = 1/2 - int_0^a eta  in  [1/2 - a eta(0), 1/2 - a eta(a)]           (0 < a < 1),
  |nu_c(0)| = I1(c) / (2 sinh c),   |nu_c(a)| = |nu_c(0)| - int_0^a mu_c  in  [|nu_c(0)| - a/2, |nu_c(0)| - a mu_c(a)].
"""
from mpmath import iv, mp, mpf

iv.prec = 320
I = iv.mpf
PI = iv.pi

def sinh(z):
    return (iv.exp(z) - iv.exp(-z)) / 2

def bessel_I(nu, z, terms=400):
    """I_nu(z) for integer nu >= 0, z > 0 an interval: power series plus a geometric tail bound."""
    half = z / 2
    term = half**nu / I(mp.factorial(nu))
    s = term
    for n in range(1, terms):
        term = term * half**2 / (n * (n + nu))
        s += term
    # tail: ratio of consecutive terms <= (z/2)^2 / (terms (terms+nu)) =: r < 1
    r = (half**2 / (terms * (terms + nu))).b
    assert r < 0.5
    tail = (term * half**2 / (terms * (terms + nu))).b / (1 - r)
    return s + I([0, tail])

def mu_nu(c, a):
    """enclosures of mu_c(a) and |nu_c(a)| for 0 < a < 1."""
    sh = sinh(c)
    eta0 = c * bessel_I(0, c) / (2 * sh)
    eta_a = c * bessel_I(0, c * iv.sqrt(1 - a**2)) / (2 * sh)
    mu = I([(0.5 - a * eta0).a, (0.5 - a * eta_a).b])
    nu0 = bessel_I(1, c) / (2 * sh)
    nu = I([(nu0 - a / 2).a, (nu0 - a * mu.a).b])
    assert mu.a > 0 and nu.a > 0
    return mu, nu

E2MODE = 'log3c'   # 'log3c' | 'thesis'

def delta(b, c, T, alpha, log3c=True, T_RH=None):
    """upper bound delta with |psi(x)-x| <= delta x for all x >= e^{b + alpha eps} (x0 = e^b);
    requires RH for 0 < gamma <= T (= c/eps).  Returns an interval; use .b."""
    c, T, alpha, b = I(c), I(T), I(alpha), I(b)
    if T_RH is not None: assert T.b <= T_RH
    eps = c / T
    assert eps.b < 1e-3 and c.a >= 3 and b.a >= mp.log(100)
    x0 = iv.exp(b)
    mu, nu = mu_nu(c, alpha)
    B0 = eps * x0 * iv.exp(-eps) * nu / (2 * mu)
    B0lo = (eps * x0 * iv.exp(-eps) * nu.a / (2 * mu.b)).a
    assert B0lo > 1
    logB0 = iv.log(I(B0lo))
    E1 = (iv.exp(2 * eps) * (eps + b) * (2 * eps * I(nu.b) / logB0 + 2.01 * eps / iv.sqrt(x0)
                                        + iv.log(iv.log(2) + 2 * b) / (2 * x0))
          + (iv.exp(alpha * eps) - 1))
    if E2MODE == 'thesis' and log3c:
        # Bhattacharjee (MSc thesis, Lethbridge 2024), Theorem 7, (2.15)-(2.16): corrected Lemma 2 / Prop. 3 / Prop. 4
        Lt = iv.log(T / PI)                                   # log(c/(pi eps))
        E2 = (I('0.16') * (1 + 3 / x0) / sinh(c) * iv.exp(I('0.71') * iv.sqrt(c * eps))
              * (Lt * iv.log(c) + 2 * Lt + I('18.75') * eps * iv.log(T) + 1))
        assert T.a >= 5000
        E3 = 2 * (1 / iv.sqrt(x0) + 1 / x0) * iv.log(T / (2 * PI))**2 / (4 * PI) + I('1.92') / x0
    else:
        E2 = (I('0.16') * (1 + 1 / x0) / sinh(c) * iv.exp(I('0.71') * iv.sqrt(c * eps))
              * (iv.log(3 * c) if log3c else 1) * iv.log(T))
        assert T.a >= 5000
        E3 = 2 / iv.sqrt(x0) * iv.log(T / (2 * PI))**2 / (4 * PI) + 2 / x0
    return iv.exp(alpha * eps) * (E1 + E2 + E3)

# Buthe 2016, Table 2:  (log of e^{alpha eps} x0, c, T, alpha, printed delta_0)
TABLE2 = [(500, 33, 2.445e12, 0.012, '1.99986e-12'), (1000, 33, 2.445e12, 0.005, '1.94751e-12'),
          (2000, 33, 2.445e12, 0.003, '1.92155e-12'), (3000, 33, 2.445e12, 0.001, '1.91298e-12'),
          (4000, 33, 2.445e12, 0.001, '1.90866e-12'), (100, 33, 2.445e12, 0.06, '2.45229e-12')]

def table2_check():
    rows = []
    for b, c, T, a, printed in TABLE2:
        without = delta(b, c, T, a, log3c=False)
        with_ = delta(b, c, T, a, log3c=True)
        rows.append((b, printed, without, with_))
    return rows

if __name__ == "__main__":
    print("Buthe 2016 Table 2 reproduced (T = 2.445e12, c = 33):")
    print("   b     printed      ours, no log(3c)    ours, with log(3c)")
    for b, printed, wo, w in table2_check():
        print(f"{b:5d}  {printed:>12}   [{mp.nstr(wo.a, 6)}, {mp.nstr(wo.b, 6)}]   "
              f"[{mp.nstr(w.a, 6)}, {mp.nstr(w.b, 6)}]")
