# Ramanujan's prime counting inequality holds for all x ≥ 38,358,837,683

**Paper:** [`paper/ramanujan_inequality.pdf`](paper/ramanujan_inequality.pdf) (28 pages; LaTeX source alongside)
**Author:** liqsweep ([@liqsweep](https://x.com/liqsweep)) · October 2026 · computer-assisted, not yet peer reviewed

In his notebooks Ramanujan showed that

$$\pi(x)^2 < \frac{e\,x}{\log x}\,\pi\!\left(\frac{x}{e}\right)$$

holds for all sufficiently large $x$. Let $N_R$ be the least integer such that it holds for every real $x \ge N_R$. In 2015 Dudek and Platt found the last counterexample below $10^{11}$, at $x = 38{,}358{,}837{,}682$, and showed that **the Riemann hypothesis implies $N_R = 38{,}358{,}837{,}683$**.

## Main result

> **Corollary F.** $N_R = 38{,}358{,}837{,}683$. The inequality holds for every real $x \ge 38{,}358{,}837{,}683$.

To our knowledge, this is the first time $N_R$ has been determined without assuming the Riemann hypothesis.

The proof combines two ingredients.

1. **A zero-free half-plane.** OpenAI's preprint of 30 September 2026 shows that $\zeta(s) \ne 0$ for $\operatorname{Re} s > 7/8$.
   - Its statement for $\zeta$ is formalized in Lean 4 as `riemannZeta_ne_zero_of_seven_eighths_lt_re`. It uses Mathlib's `riemannZeta` and only the standard axioms.
   - The formalization has been [rebuilt independently](https://github.com/tomoto0/quasi-riemann-hypothesis-7-8-verification).
   - The written proof has **not** been peer reviewed. We did not rebuild the Lean project ourselves.
2. **Theorem B (this paper).** If $\zeta(s) \ne 0$ for $\operatorname{Re} s > \theta$ with $\theta \le 0.945092$, then $N_R = 38{,}358{,}837{,}683$.
   - The proof uses only refereed inputs, together with Bhattacharjee's correction of a theorem of Büthe.
   - At $\theta = 7/8$ the bound we need is far smaller than required: at $\log x = 672$ we have $|\psi(x)-x|/x \le 5.8\cdot10^{-33}$, where about $10^{-12}$ is needed.

Four ranges cover every $x \ge 38{,}358{,}837{,}683$, with overlaps:

| Range | Proved by |
|---|---|
| $[N_R,\ 10^{19}]$ | Dudek–Platt computation; Axler (2018) |
| $[10^{19},\ e^{672.84}]$ | Proposition A′: Büthe's method with RH verified to $3\cdot10^{12}$, refereed inputs. Independently, Theorem A (Chirre–Helfgott) covers this range up to $e^{753.98}$. |
| $[e^{247.5},\ e^{3036}]$ | Theorem B with $\theta = 7/8$: the explicit formula of Cully-Hugill–Johnston plus the half-plane |
| $[e^{3035.01},\ \infty)$ | Theorem C |

## Results that use only the peer-reviewed literature (and preprints where stated)

| | Result |
|---|---|
| **A** | The inequality holds for $N_R \le x \le e^{753.98}$. The previous range ended at $e^{103}$. With refereed inputs only, the range is $e^{672.84}$. |
| **C** | It holds for $x \ge e^{3035.01}$, replacing Axler's $e^{3158.442}$, whose proof uses a defective case of a lemma. With two recent preprints the threshold is $e^{2591.01}$. |
| **B\*, B′** | Either of the following suffices: the half-plane $\operatorname{Re} s > 0.96897$ (using preprints), or a classical zero-free region $\beta \le 1 - 1/(1.7373\log\lvert\gamma\rvert)$ above $3\cdot10^{12}$. |
| **D** | RH verified up to $4.1602\cdot10^{14}$, about 139 times the current height, suffices. With refereed inputs only, $1.3277\cdot10^{15}$ suffices. |
| **E** | **Spectral barrier.** Arguments that use only positivity, analytic continuation and pole data up to height $T$ cannot go past $e^{754.40}$ when $T = 3\cdot10^{12}$. So A is optimal among such arguments to within $0.41$ in $\log x$. |
| **Trade-off** | For each zero-free-region constant $R$, the RH verification height needed to settle the inequality: for example $R = 2 \Rightarrow 5.4\cdot10^{12}$ and $R = 4 \Rightarrow 9.7\cdot10^{13}$ (Table 7). |

## The missing factor (erratum)

S. Bhattacharjee's MSc thesis (Lethbridge, 2024) found that the constant $E_2$ in Theorem 1 of Büthe (*Math. Comp.* 2016) lacks a factor $\log(3c)$, along with further corrections. Our additions are below.

- **BKLNW's table inherits the error.** The rows $40 \le b \le 2000$ of Table 8 in Broadbent–Kadiri–Lumley–Ng–Wilk (*Math. Comp.* 2021) reproduce from the printed formula, so they share the error.
  - After re-optimizing the parameters, rows 100–2000 are understated by 1.3–2.2%.
  - The rows below about $e^{55}$ are essentially unaffected.
- **Consequences in the literature.** No published theorem we examined becomes false. Several proofs need a patch or a newer reference:
  - **Johnston–Yang (2023).** Their Lemma 2.3 already leaves $(e^{2487.37}, e^{2488}]$ uncovered as printed. Corrected, it covers only $\log x \le 2484.52$. Their theorem still follows from Fiori–Kadiri–Swidinsky's Cor. 1.4.
  - **Lichtman–Pomerance (2018).** Their Prop. 2.1(iii) had 0.14% to spare, and the correction overshoots it by 0.73%.
  - **PrimeNumberTheoremAnd (Lean).** The project states Büthe's Theorem 1 without the factor, with its proof left as `sorry`.
  - **Model-based shifts.** These come from a computation with a modelled zero sum, not checked independently: Fiori–Kadiri–Swidinsky (*Res. Number Theory* 2023, Table 7), BKLNW's $M_3$ and $M_4$, and Johnston–Leong–Tudzi (2025).

  Details are in [`downstream/DOWNSTREAM_REPORT.md`](downstream/DOWNSTREAM_REPORT.md).

## Repository layout

```
paper/        ramanujan_inequality.pdf, ramanujan_inequality.tex
code/         every script used in the proofs (Python 3, mpmath; fks_eps.py and its users also need NumPy and SciPy)
outputs/      the output files that every number in the paper is taken from
downstream/   the literature audit of the Büthe/BKLNW correction: report, scripts and outputs
              (the papers themselves are not redistributed; fetch_arxiv.sh downloads them)
```

## Reproducing the results

All proofs use outward-rounded interval arithmetic (`mpmath.iv`). Run the scripts from `code/`; each writes its result to standard output. The reference outputs are in `outputs/`.

```bash
cd code
pip install mpmath numpy scipy

python qrh_ramanujan.py 200          # Theorem B at theta = 7/8: (C) on [247.5, 3036]        -> qrh_ramanujan_output.txt  (~30 s)
python qrh_theta.py 672.5            # Theorem B: largest theta = 0.945092                    -> qrh_theta_output.txt      (~2 min)
D_TOP=3035.502 python buthe_thesis_run.py   # Proposition A' (e^672.84) and refereed Theorem D -> buthe_thesis_output.txt  (~15 min)
python ch_ramanujan.py               # Theorem A (e^753.98) and B*                             -> ch_output.txt
bash regen_final.sh                  # Theorems C, C*, B', D and variants (parallel, ~30 min)  -> final_*.txt
python barrier.py                    # Theorem E / Corollary E                                 -> barrier_output.txt
python tradeoff.py 2.0 TR2.0         # one row of the trade-off table (~1 h per row)           -> tradeoff_R*.txt
python crosscheck_buthe.py           # the missing log(3c): quadrature reproduction            -> crosscheck_buthe_output.txt
```

The paper's *Data availability* section maps each output file to the theorem it supports.

## Status and caveats

- **Unreviewed.** This is a computer-assisted manuscript. Neither the paper nor the code has been independently reviewed.
- **What Corollary F rests on.** It depends on the 7/8 half-plane, whose written proof is unrefereed but whose statement for $\zeta$ is machine-checked, and on Bhattacharjee's (examined, unpublished) correction of Büthe.
  - Every other input is refereed.
  - The middle range survives even if Büthe's corrected bound were twice as weak.
  - Theorem A, which uses a different method, covers that range independently.
- **Without the half-plane.** Relying only on peer-reviewed results, the range $(e^{753.98}, e^{3035.01})$ remains open. Theorems B, B′ and D say exactly what would close it.

## Acknowledgments

The computations and the drafting were carried out with assistance from an AI system (Claude, Anthropic).

## Citation

```bibtex
@misc{liqsweep2026ramanujan,
  author = {liqsweep},
  title  = {Ramanujan's prime counting inequality holds for all $x \ge 38{,}358{,}837{,}683$},
  year   = {2026},
  note   = {Computer-assisted manuscript},
  url    = {https://github.com/blueberryvertigo/ramanujan-prime-inequality-and-missing-factor}
}
```
