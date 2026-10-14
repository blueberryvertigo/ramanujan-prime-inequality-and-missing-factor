"""Independent (non-interval) check of buthe_bound.py: mu_c, nu_c by high-precision quadrature,
Buthe 2016 Tables 1 and 2 reproduced with and without log(3c), and the Theorem A rows."""
from mpmath import mp, mpf, besseli, quad, sinh, log, exp, sqrt, pi
import buthe_bound as BB
mp.dps = 40
def eta(c, t): return c * besseli(0, c * sqrt(1 - t * t)) / (2 * sinh(c))
def mu(c, a): return mpf(1) / 2 - quad(lambda t: eta(c, t), [0, a])
def nu(c, a): return besseli(1, c) / (2 * sinh(c)) - quad(lambda s: mu(c, s), [0, a])
def delta_quad(b, c, T, a, log3c=True):
    eps = mpf(c) / T; x0 = exp(b); m = mu(c, a); n = nu(c, a)
    B0 = eps * x0 * exp(-eps) * n / (2 * m)
    E1 = (exp(2 * eps) * (eps + b) * (2 * eps * n / log(B0) + mpf('2.01') * eps / sqrt(x0)
          + log(log(2) + 2 * b) / (2 * x0)) + (exp(a * eps) - 1))
    E2 = (mpf('0.16') * (1 + 1 / x0) / sinh(c) * exp(mpf('0.71') * sqrt(c * eps))
          * (log(3 * c) if log3c else 1) * log(mpf(c) / eps))
    E3 = 2 / sqrt(x0) * log(T / (2 * pi))**2 / (4 * pi) + 2 / x0
    return exp(a * eps) * (E1 + E2 + E3)
if __name__ == "__main__":
    for a in [1e-3, 3e-3, 1e-2]:
        m_iv, n_iv = BB.mu_nu(BB.I(35), BB.I(a))
        m, n = mu(35, a), nu(35, a)
        print(f"c=35 a={a}: mu {mp.nstr(m, 10)} in enclosure: {m_iv.a <= m <= m_iv.b};  "
              f"|nu| {mp.nstr(n, 10)} in enclosure: {n_iv.a <= n <= n_iv.b}")
    for name, rows, c, T in [("Table 1", [(500, .01, '1.47067e-10'), (1000, .005, '1.43770e-10'), (3000, .001, '1.41594e-10')], 29, 3.061e10),
                             ("Table 2", [(1000, .005, '1.94751e-12'), (4000, .001, '1.90866e-12')], 33, 2.445e12)]:
        for b, a, printed in rows:
            print(f"{name} e^{b}: printed {printed}  no-log3c {mp.nstr(delta_quad(b, c, T, a, False), 6)}  "
                  f"with-log3c {mp.nstr(delta_quad(b, c, T, a, True), 6)}")
    for b in [600, 700, 1000]:
        q = delta_quad(b, 35, 3000175332800, 0.003); u = BB.delta(b, 35, 3000175332800, 0.003).b
        print(f"Theorem A row b={b}: quadrature {mp.nstr(q, 8)} <= interval bound {mp.nstr(u, 8)}: {q <= u}")
