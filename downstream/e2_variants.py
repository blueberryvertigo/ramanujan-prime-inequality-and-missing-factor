"""Buthe 2016 Thm 1 with three versions of E2 (E3 by the zero-sum model of zeros_model.py):
  'printed' : 0.16 (1+1/x0)/sinh c * e^{0.71 sqrt(c eps)} * log(c/eps)                (Thm 1 as printed)
  'log3c'   : the same times log(3c)                                                  (Prop. 3 / FKBJ Lemma 2.4)
  'thesis'  : 0.16 (1+3/x0)/sinh c * e^{0.71 sqrt(c eps)} *
              [log(c/(pi eps)) log c + 2 log(c/(pi eps)) + 18.75 eps log(c/eps) + 1]   (Bhattacharjee 2024, eq. (2.15))
For 'thesis' E3 also gets the extra 2/x0 * sum (negligible here).  delta(...) returns the bound;
optimise(...) minimises over (c, T <= H, alpha)."""
import math, sys
sys.path.insert(0, '.')
from zeros_model import munu, zsum
def e2fac(c, eps, x0inv, variant):
    if variant == 'printed': return (1 + x0inv), math.log(c / eps)
    if variant == 'log3c': return (1 + x0inv), math.log(3 * c) * math.log(c / eps)
    if variant == 'thesis':
        L = math.log(c / (math.pi * eps))
        return (1 + 3 * x0inv), L * math.log(c) + 2 * L + 18.75 * eps * math.log(c / eps) + 1
    raise ValueError(variant)
def delta(b, c, T, a, variant):
    eps = c / T; lx0 = b - a * eps
    m, n = munu(c, a)
    logB0 = math.log(eps) + lx0 - eps + math.log(n / (2 * m))
    if logB0 <= 0: return None
    ism = math.exp(-lx0 / 2); ix0 = math.exp(-lx0)
    E1 = math.exp(2 * eps) * (eps + lx0) * (2 * eps * n / logB0 + 2.01 * eps * ism
          + math.log(math.log(2) + 2 * lx0) * ix0 / 2) + math.expm1(a * eps)
    pre, fac = e2fac(c, eps, ix0, variant)
    E2 = 0.16 * pre / math.sinh(c) * math.exp(0.71 * math.sqrt(c * eps)) * fac
    E3 = 2 * (ism + (ix0 if variant == 'thesis' else 0)) * zsum(c, T) + 2 * ix0
    return math.exp(a * eps) * (E1 + E2 + E3), (E1, E2, E3)
def optimise(b, variant, H, c0=33.5, T0=None, a0=0.01):
    T0 = H if T0 is None else T0
    best = None
    cs = [c0 + 0.5 * i for i in range(-12, 13)]
    Ts = sorted(set(min(H, T0 * 10 ** (0.05 * k)) for k in range(-20, 21)))
    As = sorted(set([a0 * f for f in (0.25, 0.5, 0.75, 1.0, 1.5, 2, 3)]))
    for rnd in range(4):
        for c in cs:
            for T in Ts:
                for a in As:
                    r = delta(b, c, T, a, variant)
                    if r and (best is None or r[0] < best[0]): best = (r[0], c, T, a, r[1])
        _, c0, T0, a0, _ = best
        s = 4 ** (rnd + 1)
        cs = [c0 + 0.5 / s * i for i in range(-8, 9)]
        Ts = sorted(set(min(H, T0 * 10 ** (0.05 / s * k)) for k in range(-8, 9)))
        As = [a0 * (1 + 0.5 / s * k) for k in range(-4, 5)]
    return best
