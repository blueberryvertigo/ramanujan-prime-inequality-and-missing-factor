"""
Rigorous (interval-arithmetic) verification of Ramanujan's prime counting inequality

        pi(x)^2 < (e x / log x) * pi(x/e)                                  (R)

on ranges of L = log x, from explicit bounds on psi(x) - x.

Mode "U" (superseded, kept for reference; not run): uses Table 8 of Broadbent-Kadiri-Lumley-Ng-Wilk,
    Math. Comp. 90 (2021) 2281-2315 (arXiv:2002.11068v2), |psi(x)-x| < eps(b,b') x on
    [e^b, e^b'], computed with RH verified to H0 = 2,445,999,556,030; that height is
    covered by Platt-Trudgian's rigorous verification to 3,000,175,332,800.
Mode "Q" (conditional): assumes zeta(s) != 0 for Re(s) > 7/8 and uses
    Fiori-Jaskari, arXiv:2609.23222v2, Theorem 1.4 with sigma = 7/8, x0 = e^5:
    |psi(x)-x| < x (0.0198944 log(x)^2 x^(-1/4) + 0.02636425 x^(-1/8))  for x > e^5.

Proof skeleton (all inequalities are implemented below with outward rounding):
  * pi(x) = theta(x)/L + int_2^x theta(t)/(t log^2 t) dt  (x >= 2).
  * With E = theta - t:  pi(x) = J(x) + E(x)/L + int_2^x E(t)/(t log^2 t) dt,
    J(x) = x/L + int_2^x dt/log^2 t.  So |pi(x) - J(x)| <= Delta(x).
  * 0 <= psi - theta <= K (sqrt x + x^(1/3) log_2 x), K = 1.01624,
    from theta(y) < 1.01624 y for all y > 0 (Rosser-Schoenfeld 1962, Thm 9).
  * |E(t)| <= t trivially on [2, e^20] (0 <= theta(t) < 1.01624 t).
  * (R) holds at x if  (j(L-1) - d(L-1))/L - (j(L) + d(L))^2 > 0,
    where j(l) = J(e^l)/e^l and d(l) = Delta(e^l)/e^l.
  * j is expanded exactly:  j(l) = P(1/l) + rho(l),  P(u) = sum_{k<n} k! u^{k+1},
    and P(1/(L-1))/L - P(1/L)^2 = u^6 N6(u)/(1-u)^n with an exact integer polynomial N6,
    so the x^2/L^6 main term of (R) is extracted without cancellation.
"""
import sys
from math import factorial
from fractions import Fraction as Fr
from mpmath import iv, mpf, mp

iv.prec = 320
mp.prec = 320
I = iv.mpf
LOG2 = iv.log(I(2))
E1 = iv.exp(I(1))

# ---------------- exact polynomial part ----------------
n = 12          # expansion length; m = n+1 = 13 <= s0
s0 = 20         # J split point: t = e^20
theta_split = I(4)/5

def padd(a, b):
    r = [0]*max(len(a), len(b))
    for i, x in enumerate(a): r[i] += x
    for i, x in enumerate(b): r[i] += x
    return r
def pmul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] += x*y
    return r
def ppow(a, e):
    r = [1]
    for _ in range(e): r = pmul(r, a)
    return r

P = [0]*(n+1)
for k in range(n): P[k+1] = factorial(k)
omu = [1, -1]
term1 = [0]
for k in range(n):
    t = [0]*(k+2) + [factorial(k)]        # k! u^{k+2}
    term1 = padd(term1, pmul(t, ppow(omu, n-1-k)))
Npoly = padd(term1, [-c for c in pmul(ppow(omu, n), pmul(P, P))])
assert all(c == 0 for c in Npoly[:6]), Npoly[:8]
N6 = Npoly[6:]
# sanity: series of N/(1-u)^n should reproduce Hassani's -d_k, k=4..n-2
def hassani_d(k):
    return sum(factorial(j)*(factorial(k-j) - (factorial(k)//(factorial(j)*factorial(k-j)))) for j in range(k+1))
ser = [Fr(0)]*(n)
inv = [Fr(comb_) for comb_ in [1]]  # series of (1-u)^(-n)
from math import comb
inv = [Fr(comb(n-1+i, i)) for i in range(n)]
for i in range(n):
    ser[i] = sum(Fr(N6[j])*inv[i-j] for j in range(min(i+1, len(N6))))
for k in range(4, n-2):
    assert ser[k-4] == -hassani_d(k), (k, ser[k-4], -hassani_d(k))

def horner(coeffs, u):
    acc = I(0)
    for c in reversed(coeffs):
        acc = acc*u + c
    return acc

Pf = lambda u: horner(P, u)

# ---------------- rho(l): rigorous enclosure of j(l) - P(1/l) ----------------
# K0 = int_2^{e^20} dt/log^2 t  in [0, e^10/log^2 2 + e^20/100]
K0 = I([0, (iv.exp(I(10))/LOG2**2 + iv.exp(I(20))/100).b])
B0 = iv.exp(I(s0))*sum(I(factorial(k))/I(s0)**(k+1) for k in range(1, n))
m = n+1
NF = I(factorial(n))
def rho(l):
    # n! e^{-l} int_{s0}^{l} e^s s^{-m} ds, g increasing on [s0,l] since m <= s0
    tl = theta_split*l
    q = I(m)/tl
    assert q.b < 1 and tl.a > s0
    gl = l**(-m)                                  # e^{-l} g(l)
    gt = iv.exp(-(1-theta_split)*l)*tl**(-m)      # e^{-l} g(theta l)
    lo = gl - gt
    hi = gl/(1-q) + (tl - s0)*gt
    rem = NF*I([lo.a, hi.b])
    return rem + iv.exp(-l)*(K0 - B0)

# ---------------- d(l) : bound for |pi - J| / x ----------------
K = I('1.01624')                 # Rosser-Schoenfeld: theta(y) < 1.01624 y for all y > 0
T8 = [(20,'4.26760e-5'),(21,'2.58843e-5'),(22,'1.56996e-5'),(23,'9.52229e-6'),(24,'5.77556e-6'),
      (25,'3.50306e-6'),(30,'2.87549e-7'),(35,'2.36034e-8'),(40,'1.93378e-8'),(45,'1.09073e-8'),
      (50,'1.11990e-9'),(100,'2.45299e-12'),(200,'2.18154e-12'),(300,'2.09022e-12'),(400,'2.03981e-12'),
      (500,'1.99986e-12'),(600,'1.98894e-12'),(700,'1.97643e-12'),(800,'1.96710e-12'),(900,'1.95987e-12'),
      (1000,'1.94751e-12'),(1500,'1.93677e-12'),(2000,'1.92279e-12'),(2500,'9.06304e-13'),(3000,'4.59972e-14')]
def bump(s):  # add one unit in the last printed digit (table values are rounded)
    mant, ex = s.split('e'); dec = len(mant.split('.')[1])
    return I(s) + I(10)**(int(ex)-dec)
SEG = [(T8[i][0], T8[i+1][0] if i+1 < len(T8) else 3500, bump(T8[i][1])) for i in range(len(T8))]
def seg_of(l):
    for a, b, e in SEG:
        if a <= l.a and l.b <= b: return a, b, e
    raise ValueError(f"interval {l} straddles a table boundary")

def psitheta_terms(l):
    # pointwise: (psi-theta)/x <= K(e^{-l/2} + e^{-2l/3} l/log2)
    pt = K*(iv.exp(-l/2) + iv.exp(-2*l/3)*l/LOG2)
    # (1/x) int_2^x (psi-theta)(t)/(t log^2 t) dt
    #   <= K/x [ 2x^{1/4}/log^2 2 + 8 sqrt(x)/l^2 + (1/log2)(3x^{1/6}/log2 + 6x^{1/3}/l) ]
    it = K*(2*iv.exp(-3*l/4)/LOG2**2 + 8*iv.exp(-l/2)/l**2
            + (3*iv.exp(-5*l/6)/LOG2 + 6*iv.exp(-2*l/3)/l)/LOG2)
    return pt, it

def d_U(l):
    a, b, e = seg_of(l)
    pt, it = psitheta_terms(l)
    point = (e + pt)/l
    integ = K0*iv.exp(-l)                                     # t in [2, e^20]
    for (aj, bj, ej) in SEG:                                  # completed earlier segments
        if bj <= a:
            integ += ej*iv.exp(I(bj)-l)/(I(bj)**2*(1-I(2)/aj))
    integ += e/(l**2*(1-I(2)/a))                              # current segment up to x
    return point + integ + it

A_FJ, C_FJ = I('0.0198944'), I('0.02636425')
def d_Q(l):
    pt, it = psitheta_terms(l)
    point = (A_FJ*l**2*iv.exp(-l/4) + C_FJ*iv.exp(-l/8) + pt)/l
    integ = K0*iv.exp(-l)
    integ += A_FJ*I(4)/3*iv.exp(-l/4)
    integ += C_FJ*(iv.exp(-l/8)/(l**2*(I(7)/8 - 4/l)) + I(8)/7/400*iv.exp(-9*l/16))
    return point + integ + it

def margin(L, dfun):
    """rigorous lower bound for L^6 * [ (j(L-1)-d(L-1))/L - (j(L)+d(L))^2 ] over interval L"""
    u = 1/L
    main = N6val = horner(N6, u)/(1-u)**n                     # = L^6 * (P(1/(L-1))/L - P(1/L)^2)
    r1, r0 = rho(L-1), rho(L)
    Pu = Pf(u)
    corr = (r1/L - 2*Pu*r0 - r0**2)*L**6
    d1, d0 = dfun(L-1), dfun(L)
    err = (d1/L + 2*(Pu + r0)*d0 + d0**2)*L**6
    return main + corr - err, N6val, corr, err

def grid(lo, hi, step, mode):
    cuts = set()
    if mode == 'U':
        for a, b, e in SEG: cuts |= {a, a+1, b, b+1}
    pts = [lo]
    x = lo
    while x < hi:
        nxt = min(hi, x+step)
        for c in sorted(cuts):
            if x < c < nxt: nxt = c; break
        pts.append(nxt); x = nxt
    return list(zip(pts[:-1], pts[1:]))

def run(mode, lo, hi, step, stop_on_fail=True, report_every=None):
    dfun = d_U if mode == 'U' else d_Q
    worst = None; last_ok = lo; nint = 0
    for a, b in grid(lo, hi, step, mode):
        L = I([a, b])
        mg, main, corr, err = margin(L, dfun)
        nint += 1
        if worst is None or mg.a < worst[0]: worst = (mg.a, a, b, main, corr, err)
        if mg.a <= 0:
            print(f"  [{mode}] FAIL on L in [{a}, {b}]: margin lower bound {mp.nstr(mg.a, 6)} "
                  f"(main {mp.nstr(main.a,6)}, corr {mp.nstr(corr.a,3)}, err {mp.nstr(err.b,6)})")
            if stop_on_fail: return last_ok, worst, nint
        else:
            last_ok = b
        if report_every and nint % report_every == 0:
            print(f"  [{mode}] L in [{a:.2f},{b:.2f}]  margin >= {mp.nstr(mg.a, 6)}  (main {mp.nstr(main.a,6)}, err <= {mp.nstr(err.b,4)})")
    return last_ok, worst, nint

if __name__ == "__main__":
    print("exact N6 polynomial: leading coeffs", N6[:4], "(series matches Hassani d_k for k=4..%d)" % (n-3))
    print("K0 <=", mp.nstr(K0.b, 6), "  B0 =", mp.nstr(B0.a, 10))
    # Mode U (BKLNW Table 8 as printed) is no longer run: rows b <= 2000 of that table come from
    # Buthe 2016, Thm 1 as printed, whose E2 lacks a factor log(3c) (see buthe_bound.py).
    # Theorem A is proved in buthe_ramanujan.py.
    # ---- conditional on Re(rho) <= 7/8: where does the F-J bound alone start to suffice? ----
    last_fail = None
    for a, b in grid(100, 200, 0.25, 'Q'):
        if margin(I([a, b]), d_Q)[0].a <= 0: last_fail = b
    print(f"[Q] F-J bound alone fails somewhere below L = {last_fail}; passes on every interval above it in [100,200]")
    qstart = 200
    okq, worstq, cntq = run('Q', qstart, 3158.5, 0.5, report_every=1500)
    print(f"[Q] verified contiguously for L in [{qstart}, {okq}]  ({cntq} intervals); "
          f"smallest margin {mp.nstr(worstq[0],6)} on [{worstq[1]},{worstq[2]}]")
