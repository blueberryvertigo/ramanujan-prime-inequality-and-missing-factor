# Downstream effects of the missing log(3c) in Büthe's Theorem 1 (E₂)

*13 October 2026. This is read-only literature research plus computation. All numbers come from the scripts in this folder (listed at the end).*

**Title note.**
- arXiv:1410.7015 is J. Büthe, *Estimating π(x) and related functions under partial RH assumptions*, Math. Comp. 85 (2016) 2483–2498. Theorem 1 and Tables 1–2 are in this paper.
- *An analytic method for bounding ψ(x)* is a different paper: arXiv:1511.02032, Math. Comp. 87 (2018). It is the source of the 0.94√x bound and of BKLNW Table 8 rows b ≤ 35. Neither of those is affected.

## 0. Findings that change the framing

1. **Prior art.** S. Bhattacharjee's MSc thesis (*A survey of Büthe's method for estimating prime counting functions*, U. Lethbridge, 2024) already reports that Büthe's E₂ is missing log(3c). The supervisor was N. Ng; the committee included H. Kadiri and A. Fiori.
   - The thesis also claims a further error in Lemma 2, which comes from FKBJ Lemma 2.4. It says an integral is computed as log(2c/2π)/(2πc) where it should be log(2c/(2πε))/(2πc).
   - It replaces log(3c)·log(c/ε) by log(c/πε)·log c + 2 log(c/πε) + 18.75 ε log(c/ε) + 1, and (1 + 1/x₀) by (1 + 3/x₀). See its eq. (2.15).
   - At c = 33 this makes E₂ about 15% larger than the log(3c)-only version.
   - Its Table 2.7 recomputes Büthe's Table 2 at Büthe's parameters. Example: e^1000 gives 2.13205e-12, against 1.94751e-12 printed.
   - My implementation of that formula (`e2_variants.py`, variant `thesis`) reproduces Table 2.7 to within 0.06%.
   - **I have not checked the Lemma 2 claim independently.** So every margin below is computed twice: with log(3c) only ("log3c"), and with the thesis E₂ ("thesis").
   - The PrimeNumberTheoremAnd blueprint cites this thesis as a resource. Yet it still formalizes the printed Theorem 1.
2. **Where the omission matters.** It only matters where the optimal truncation T = c/ε reaches the RH height H. Below that, c can be raised at almost no cost, so E₂ becomes negligible. As a result:
   - BKLNW Table 8 rows b = 40 (= ε(19 log 10)), 45 and 50 are essentially unaffected: at most +0.03% after re-optimization.
   - FKS2 Table 2 rows log x₀ = 44…57 are essentially unaffected, as is Büthe's Table 1 row e^45.
   - The effect sets in around log x ≈ 55–60, and from about 65 upward it is +1.3% to +2.5% after re-optimization.
   - So BKLNW Cor. 2.1, θ(x) ≤ (1 + 1.93378·10⁻⁸)x, is **not** affected. Re-optimized with log(3c), ε(10¹⁹) is 1.91948e-8 < 1.93378e-8. This is the value cited by about 17 of the papers checked.
3. **Two more independent uses of the printed formula** (outside Büthe's and BKLNW's tables):
   - FKS2 (arXiv:2206.12557), Tables 1–4. I reproduced their Table 2 rows from the (α, c, T) they print, to within 0.2%, without log(3c).
   - Johnston–Yang (arXiv:2204.01980), Lemma 2.3.
4. **Least corrected values** (`variants_table_output.txt`):

| row | printed | log3c, re-opt. | thesis E₂, re-opt. |
|---|---|---|---|
| Büthe T1 e^55 (H = 3.061e10) | 2.88434e-10 | 2.90948e-10 (+0.87%) | 2.91436e-10 (+1.04%) |
| Büthe T1 e^500 | 1.47067e-10 | 1.50799e-10 (+2.54%) | 1.51156e-10 (+2.78%) |
| Büthe T2 e^1000 | 1.94751e-12 | 1.98965e-12 (+2.16%) | 1.99376e-12 (+2.37%) |
| BKLNW T8 b=100 | 2.45299e-12 | 2.50684e-12 (+2.20%) | 2.51203e-12 (+2.41%) |
| BKLNW T8 b=200 | 2.18154e-12 | 2.20905e-12 (+1.26%) | 2.21361e-12 (+1.47%) |
| BKLNW T8 b=2000 | 1.92279e-12 | 1.96220e-12 (+2.05%) | 1.96625e-12 (+2.26%) |
| BKLNW extended ε(2300) | 1.91944e-12 | 1.95873e-12 (+2.05%) | 1.96278e-12 (+2.26%) |
| JY Lemma 2.3, b=2000, H=3.0002e12 | 1.570e-12 | 1.60480e-12 (+2.22%) | 1.60810e-12 (+2.43%) |
| FKS2 T3/T4 region, b=1350, H=3e12 | 1.5795e-12 | 1.61545e-12 (+2.28%) | 1.61877e-12 (+2.49%) |
| BKLNW T8 b=40 (ε(10¹⁹)) | 1.93378e-8 | 1.91948e-8 (still below) | 1.91948e-8 |

At the published (fixed) parameters, the excess with log(3c) is +5.2% to +5.5% for log x ≳ 70. With the thesis E₂ it is +6.2% to +6.5%, and +9.2% to +9.7% at Büthe's c = 33.

## 1. Main table

Verdicts:
- **a**: survives with room.
- **b**: survives, but a stated constant or threshold shifts slightly.
- **c**: the conclusion as stated loses its proof.
- **unverified**: possibly affected, not checked.

"Re-opt" means the corrected bound re-optimized over (c, α, T ≤ H).

| Paper | Affected values used | Where | Margin analysis | Verdict | Conf. |
|---|---|---|---|---|---|
| **Johnston–Yang 2023**, JMAA (arXiv:2204.01980v2) | Own recomputation of Büthe Thm 1 (printed E₂, α = 0, c = 33.6, H = 3.0002e12): 1.570e-12 for x ≥ e^2000. Also BKLNW T8 rows 50–2000. | Lemma 2.3; §3.1 (covers e^2000 < x ≤ e^2488 for Thm 1.1); §3.2 starts at x₀ = e^2488 | Thm 1.1 needs ε ≤ f(L) = 9.39 L^1.515 e^{−0.8274√L}, and f(2488) = 1.5624e-12. **As printed, 1.570e-12 only covers L ≤ 2487.37**, which is already a small gap. Corrected and re-optimized: 1.6048e-12 (log3c) or 1.6081e-12 (thesis), covering L ≤ 2484.5 or 2484.2. At JY's own parameters: 1.6543e-12, covering L ≤ 2480.6. Uniform +2.6% → 2484.0; +8% → 2477.4. Uses of BKLNW T8 have ≥ 40× room. | **c** (proof gap on (e^2484.5, e^2488]). Statement still true: FKS1 Cor. 1.4 implies Thm 1.1 for all x ≥ 2. | high |
| **Lichtman–Pomerance 2018**, JNT (arXiv:1705.02442) | Büthe Table 1 (H = 3.061e10) rows e^45, e^50, e^55: 1.118e-8, 1.165e-9, 2.885e-10 | Prop. 2.1; used in the W* ranges, §5 | (iii) \|θ−x\|/x < 2.9e-10 for x > e^55: least corrected ε_ψ is 2.9095e-10 (log3c) or 2.9144e-10 (thesis). Adding ψ−θ ≤ 1.02√x + 3x^{1/3} gives 2.921e-10 or 2.926e-10, so it **fails by 0.73% / 0.90%**. (ii) 1.2e-9 at e^50: ε_θ ≤ 1.1812e-9, room 1.6% (fails only under a uniform +2.6%). (ii′) 1.2e-8 at e^45: room 7.7%. (i) 2.3e-8 for x > 10¹⁹: large room. | **c** for Prop. 2.1(iii). Statement true via FKS2 Table 2 (H = 3e12): ε_θ(e^55) = 1.137e-10, a row with T < H, so unaffected. | high |
| **BKLNW 2021**, Math. Comp. (arXiv:2002.11068); its own derived results | Table 8 b ≥ 100; extended Table 8 near b = 2300 | Tables 8, 10, 14, 15; Cor. 11.1 | **Unaffected:** Cor. 2.1 (1.93378e-8), Table 8 b = 40/45/50, Cor. 5.1. **Affected, b:** Table 8 b = 100–2000 (+1.3% to +2.5%); Table 14 M₀/m₀ for log X = 60–2000. **Table 15 middle-range maxima** come from B_k(2300, 2325) = ε(2300)·2325^k (Büthe-based) or B_k(2325, 2350) (Platt–Trudgian-based). Corrected ε(2300) is ×1.0205 (log3c) or ×1.0226 (thesis). This gives M₁ 4.4627e-9 → 4.554e-9 / 4.564e-9 and M₂ 1.0376e-5 → 1.0588e-5 / 1.0610e-5. **M₃ = 0.024334 (Cor. 11.1) and M₄ = 57.184 are exceeded on BKLNW's own 25-step grid** (0.024618 / 0.024669; 57.236 / 57.354). M₅ is fine. On a fine grid the max moves to the Büthe/PT crossover: ≈ 0.02399 and 55.3 (estimate; extrapolates their PT values). | **b** (M₁, M₂, M₀ rows shift about 2%). Derivation of M₃ and M₄ fails as printed but looks recoverable. | high / medium (fine-grid part) |
| **Fiori–Kadiri–Swidinsky 2023b**, Res. Number Theory (arXiv:2206.12557v2) ("FKS2") | Own recomputation of Büthe Thm 1 (printed E₂, H = 3e12) for 10¹⁹ ≤ x₀ < e^2073. Reproduced from their printed (α, c, T) to within 0.2%. | Tables 2–4 (and Table 5 rows ≤ 2000); Cor. 22–26; Table 7 | Rows 44–57 unaffected (≤ 0.1%). From about 60 upward, +0.5% to +2.5% (re-opt) or up to +6.5% (fixed). **Table 7 range ends** (corrected re-opt, log3c/thesis), shown as stated → corrected: (log x)^{3/2}x^{−1/2}: 65.65 → 65.62; (log x)²x^{−1/2}: 70.6 → 70.56; (log x)³x^{−1/2}: 80 → 79.95; x^{−1/3}: 80.55 → 80.50; x^{−1/10}: 270.8 → 270.66–270.68; **x^{−1/50}: 1358.6 → 1357.4–1357.5**. Survive with thin room: 2 log x·x^{−1/2} ≤ 57 (58.11); (1/8π)log²x·x^{−1/2} ≤ 60.8 (60.92); x^{−1/4} ≤ 107.6 (107.66); x^{−1/5} ≤ 134.8 (134.83). x^{−1/100} ≤ 3757.6 is unaffected. Table 6 / Cor. 23 have room ≥ 58% (worst row A = 1.2e-4); Cor. 22 (9.2211) has room ×10.9; Cor. 26 (0.4298) is trivial. | **b** (6 of the 11 Table 7 ranges shrink by 0.04–1.2 in log x; tables understated by ≤ 2.5%) | high (shifts) / medium (exact ends) |
| **Johnston–Leong–Tudzi 2025**, Exp. Math. (arXiv:2408.04143v3) | BKLNW T8: B₄ = 75.1090 = ε(2000, 2500)·2500⁴ (the b = 2000 row). B₁–B₃ come from b = 20 (unaffected). | Lemma A.5 → Table 1 → Thm A.1 → Thm A.2 → Thm 1.2 | Corrected B₄ = 76.65 / 76.81. Re-running their Prop. A.6/A.7 (my run reproduces their cM to ≤ 0.01%): row 1 cM 7.80973e-3 → 7.8415e-3. **(A.5)** 2.91890 → 2.9276 / 2.9285 (+0.3%). **(A.6)** 4.88346 → 4.9143 / 4.9175 (+0.6%). Effect on m(x) (4.591, 7.397) and on \|W(x)\| < 2260x not recomputed. | **b** (A.5, A.6 shift). 2260: unverified, small sensitivity. Using FKS2 Table 3 instead gives B₄ ≈ 30, so all constants are recoverable. | high (A.5/A.6) |
| **Bordignon–Starichkova 2024**, Ramanujan J. (arXiv:2211.08844) | BKLNW Table 15 M₁ = 4.4627e-9 for x > e^1000 | §2, bound on Σ1/p | 4.47e-9 becomes 4.56e-9 / 4.57e-9. Immaterial to the final GRH–Chen results (not traced in detail). | **b** (intermediate constant); final results a | medium |
| **Axler 2022a** (arXiv:2203.05917) | BKLNW Table 15: 0.024334 (k = 3, x ≥ e^29), 57.184 (k = 4) | Prop. 1.1 (N = 1,757,126,630,797), Lemma 3.1, Thms 1.3–1.4 | The inputs keep their values only after a fine-grid recomputation or with FKS-era inputs (M₃ ≈ 0.015, M₄ ≈ 30). N is computed below e^29 and is unaffected. | **a** (input needs re-justification) | medium |
| Axler 2022b (arXiv:2203.09290); Axler 2022c, Robin (arXiv:2110.13478); Axler 2024 (arXiv:2406.04018) | Same: 0.024334 and/or 57.184 | Various lemmas | As above | **a** (same caveat) | medium |
| **Fiori–Kadiri–Swidinsky 2023a**, JMAA (arXiv:2204.02588) ("FKS1") | BKLNW T8: ε(e^40) = 1.93378e-8, ε(e^1000) = 1.94751e-12 | Lemma 5.2 → Cor. 1.4 (9.22022) | 1.93378e-8 is unaffected, with ×7.5 room against 1.443e-7. ε(e^1000) corrected to 1.989e-12 / 1.993e-12, still below 2.593e-12 (room 23%). Survives +8%. | **a** | high |
| **PrimeNumberTheoremAnd / IEANTN** (Lean, GitHub, Oct 2026) | `BKLNW_app.bklnw_thm_16` states Büthe Thm 1 with the printed E₂ (`sorry`); `theorem_2`/`table_8_ε`; FKS2 Cor. 24 "trusted" windows | BKLNW_app.lean; FKS2Cor24TrustedNumerics.lean | The trusted windows at the Table 7 ends lose their numerical support on: row 2 (65.62, 65.65]; row 4 (70.56, 70.6]; row 5 (79.95, 80]; row 6 (80.50, 80.55]; row 9 (270.66, 270.8]; and **row 10, the whole window [e^1358, e^1358.6]**. Ramanujan.lean uses only BKLNW T15 row 43, k = 2 (unaffected). | **c** (for these trusted facts as justified by FKS2) | high |
| Ramaré–Zuniga Alterman 2025, Funct. Approx. (arXiv:2408.05969) | BKLNW T8 rows 21 and 40; T15 row J | (13), (14), (16) | Unaffected rows | a | high |
| Johnston–Ramaré–Trudgian 2023, RNT (arXiv:2303.13785) | T8 row 40 (as 1.994e-8); T15 row 10⁵ | Lemmas 6, 8 | Unaffected | a | high |
| Hathi–Johnston 2024, JNT (arXiv:2109.11883) | T15 rows J and 10⁹ | Lemmas 5.1, 5.2 | Unaffected | a | high |
| Johnston 2022, Ramanujan J. (arXiv:2109.02249) | Uses Büthe **Prop. 3 including log(3c)** (Lemma 2.1), not Thm 1; BKLNW Cor. 5.1 | §§2–5 | Independent confirmation of the factor. Unaffected by log(3c). Would be touched only if Bhattacharjee's Lemma 2 claim holds (see §4). | a | high |
| Axler 2022, Math. Inequal. Appl. (arXiv:2207.02486) | FKS2 Cor. 14 (from FKS1; survives) | Lemma 2.1 | x ≥ e^3158, outside the affected range | a | high |
| Visser 2022 ×2 (arXiv:2211.00840, 2211.06469), Visser 2025 (arXiv:2504.14458) | FKS1/FKS2/JY asymptotic forms (9.22, 9.2211, 9.39/9.40) | throughout | The statements used are still true (JY Thm 1.1 follows from FKS1 Cor. 1.4) | a | medium |
| Cully-Hugill–Johnston 2021/2024 (arXiv:2111.10001, 2402.04272); Cully-Hugill–Lee 2022 (arXiv:2103.05986); Cully-Hugill 2023 (arXiv:2107.14468); Cully-Hugill–Dudek (arXiv:1902.05065); Lee 2023/2024/2026 (arXiv:2301.13457, 2403.11814, 2602.14340); Lee–Nosal (arXiv:2312.05628); Johnston–Starichkova (arXiv:2208.01229) | BKLNW Cor. 5.1 / 2.1: a₁ = 1 + 1.93378e-8 (or 1 + ε(500) = 1 + 1.99986e-12) and a₂ | ψ − θ bounds | 1.93378e-8 is unaffected. The ε(500) change moves a₁ by 4e-14. | a | high |
| Garcia–Lee ×2 (arXiv:2007.10313, 2012.03102); Francis–Lee (arXiv:2003.08083); Lee–O'Clarey (arXiv:2506.11814); Teo 2025 (arXiv:2508.06764); Bordellès (arXiv:2402.12333); Yamada ×2 (arXiv:2103.16936, 2303.16853); Fan (arXiv:2306.03347); Fan–Pomerance (arXiv:2306.03339) | Cor. 2.1, Thm 1 at small X, T15 small-X rows, Cor. 11.2 | various | Unaffected values, or room of orders of magnitude | a | high |
| Nguyen–Pomerance 2019, Math. Comp.; Pomerance–Spicer 2019, AMM; Lichtman 2018 (arXiv:1801.01925); Lichtman–Pomerance 2018 (arXiv:1609.05569) | Only the 2.3e-8 (x > 10¹⁹) part of LP Prop. 2.1 / Büthe | — | Large room | a | high |
| Platt–Trudgian 2021 Math. Comp. (arXiv:1809.03134); Platt–Trudgian 2021 BLMS (arXiv:2004.09765); Morrill–Platt (arXiv:1809.10813, uses Dusart's Table 1) | Büthe Thm 2 or a remark only | — | No dependence on Thm 1 or Tables 1–2 | a (not dependent) | high |
| Axler–Nicolas 2023, Acta Arith. 209; Nicolas 2021, Ramanujan J.; Chalker 2016 (ANU thesis); Bërdëllima 2025, NNTDM; Ralaivaosaona–Razakarinoro 2026, JNT; Fan–Kobayashi–Molnar 2026 | Cite Büthe 2016 and/or BKLNW/LP; full text not obtained | — | — | unverified | — |

Of the other works in the Semantic Scholar citation lists, about 25 cite only for background, with no numerical use found. They are listed in §5.

## 2. Details of the c cases

### 2.1 Johnston–Yang, Lemma 2.3 and Theorem 1.1 (`jy_lemma23.py`, `summary_calcs.py` §A)

- Their Theorem 1.1 states |ψ(x) − x| ≤ 9.39 x (log x)^1.515 exp(−0.8274 √log x) for all x ≥ 2.
- §3.1 covers e^2000 < x ≤ e^2488 by Lemma 2.3 (1.570e-12). §3.2 handles x ≥ x₀ = e^2488. So §3.1 needs 1.570e-12 ≤ f(L) for 2000 < L ≤ 2488, and f is decreasing there.

| Value used for Lemma 2.3 | ε | f(L) ≥ ε for L ≤ |
|---|---|---|
| as printed | 1.570e-12 (computed 1.56936e-12) | 2487.37 (2487.43) |
| log(3c), at JY's c = 33.6, α = 0 | 1.65430e-12 | 2480.57 |
| log(3c), re-optimized | 1.60480e-12 | 2484.52 |
| thesis E₂, re-optimized | 1.60810e-12 | 2484.25 |
| uniform +2.6% / +8% | 1.6108e-12 / 1.6956e-12 | 2484.03 / 2477.36 |

- f(2488) = 1.56239e-12, so the printed argument already leaves (e^2487.4, e^2488] uncovered. The correction widens the gap to about 3.5–3.8 units of log x.
- Bridging the gap from their §3.2 would need A′ recomputed at x₀ ≈ e^2484, which the paper does not give. To keep the constants B and C, A would have to be 9.645 at L = 2488 with the corrected re-optimized value (9.436 already with the printed value).
- **Repair.**
  - FKS1 Cor. 1.4 (9.22022 (log x)^{3/2} exp(−0.8476836 √log x), x > 2) is at most the JY bound for every L = log x ≥ log 2: the ratio is 0.9819·L^{−0.015}·e^{−0.0203√L} < 1 there.
  - FKS1 Cor. 1.4 itself survives (§1, FKS1 row). So Theorem 1.1 and Corollaries 1.2–1.3 remain true; only JY's own proof has the gap.
- This is the arXiv v2 text; I did not see the journal version.
- Lemma 2.3 as a statement: if quoted as "1.570e-12 for x ≥ e^2000", the provable value is ≥ 1.605e-12 (with H = 3e12, by Büthe's method).

### 2.2 Lichtman–Pomerance, Prop. 2.1(iii) (`lp_prop21.py`, `lp_e55_scan.py`, `summary_calcs.py` §B)

- Prop. 2.1 claims |θ(x) − x|/x < 2.9e-10 for x > e^55. The proof takes 2.885e-10, which is Büthe Table 1, e^55 (2.88434e-10, H = 3.061e10). It adds ψ − θ ≤ 1.02 x^{1/2} + 3 x^{1/3}, which gives +1.163e-12.
- My zero-sum model (`zeros_model.py`) reproduces Büthe's printed value at his parameters (c = 30, α = 0.1) to 0.00%.
- A broad scan (c ∈ [20, 40], α ∈ [0.01, 0.4], T ≤ 3.061e10) gives a least corrected ε_ψ of 2.90948e-10 (log3c) or 2.91436e-10 (thesis). Then ε_θ ≤ 2.9211e-10 or 2.9260e-10, which is **above 2.9e-10 by 0.73% / 0.90%**.
- The model's uncertainty is the S(t) term of the zero sum, about 0.1% of the total. That is far below this margin.
- **Repair:** use H = 3·10¹² (Platt–Trudgian 2021). FKS2 Table 2, log x₀ = 55, has T = 8.9e11 < H, so the E₂ change there is negligible, and gives ε_θ = 1.1368e-10.
- I did not trace the paper's final Ψ(x, y) theorems. They use (iii) only through w > e^55, and (iii) is true, so I expect them to stand.
- The other users of Prop. 2.1 (Lichtman 2018, Nguyen–Pomerance, Pomerance–Spicer, Fan–Pomerance) cite only the 2.3e-8 part.

### 2.3 PrimeNumberTheoremAnd / IEANTN

- `PrimeNumberTheoremAnd/IEANTN/BKLNW/BKLNW_app.lean`: `bklnw_thm_16` gives E₂ = 0.16(1 + x₀⁻¹)/sinh c · e^{0.71√(cε)} · log(c/ε). There is no log(3c), and the proof is `sorry`. `theorem_2` (Table 8) is also `sorry`. Its blueprint note points to Bhattacharjee's thesis.
- `FKS2Cor24TrustedNumerics.lean` has 11 `sorry`'d "sliver/band" facts at the Table 7 thresholds. These rest on FKS2's (understated) Büthe-based numbers. With the corrected inputs the following parts are unsupported:
  - band_row2: (e^65.62, e^65.65]
  - band_row4: (e^70.56, e^70.6]
  - the row-5 window above e^79.95
  - sliver_row6: (e^80.50, e^80.55]
  - sliver_row9: (e^270.66, e^270.8]
  - **sliver_row10: all of [e^1358, e^1358.6]**, since the corrected end is ≈ 1357.4–1357.5
- These statements are very likely true, but the cited computation no longer proves them.
- With refereed inputs, Büthe's method is still the best available at L ≈ 1358, and FKS1's Theorem 1.1 is worse there. The unrefereed Chirre–Helfgott bound (arXiv:2512.15709, ε ≈ π/(3·10¹²) ≈ 1.05e-12) would cover row 10.
- I found no PNT+ issue mentioning log(3c) or the Büthe E₂ (GitHub search, 13 Oct 2026).

## 3. Details of the b cases

### 3.1 FKS2 Table 7 (`fks2_table7_v2.py` → `fks2_table7_v2_output.txt`)

Method:
- ε_π,num is interpolated log-linearly from FKS2 Table 4.
- The crossing L* where B(e^L) = ε_π is found by scan and bisection.
- The corrected ε_ψ increase at L* is added to ε_π, but only for L < 2073 (above that FKS2 use FKS1).
- My reproduction of each stated end is within +0.004 to +0.19 of the stated value. So the shifts, not the absolute ends, are the robust quantity.

| B(x) | stated end | reproduced | log3c fixed / re-opt | thesis fixed / re-opt | +2.6% | +8% |
|---|---|---|---|---|---|---|
| 2 log x·x^{−1/2} | 57 | 58.12 | 58.10 / 58.11 | 58.09 / 58.11 | 56.90* | 54.13* |
| (log x)^{3/2} x^{−1/2} | 65.65 | 65.667 | 65.557 / **65.619** | 65.536 / **65.615** | 65.596 | 65.456 |
| (1/8π)(log x)² x^{−1/2} | 60.8 | 60.986 | 60.852 / 60.927 | 60.827 / 60.922 | 60.683* | 60.077* |
| (log x)² x^{−1/2} | 70.6 | 70.608 | 70.498 / **70.561** | 70.476 / **70.556** | 70.552 | 70.439 |
| (log x)³ x^{−1/2} | 80 | 80.004 | 79.892 / **79.956** | 79.870 / **79.952** | 79.948 | 79.836 |
| x^{−1/3} | 80.55 | 80.575 | 80.418 / **80.508** | 80.387 / **80.502** | 80.497 | 80.341 |
| x^{−1/4} | 107.6 | 107.759 | 107.549 / 107.670 | 107.509 / 107.661 | 107.656 | 107.448 |
| x^{−1/5} | 134.8 | 134.949 | 134.686 / 134.837 | 134.636 / 134.827 | 134.820 | 134.561 |
| x^{−1/10} | 270.8 | 270.907 | 270.379 / **270.682** | 270.277 / **270.662** | 270.650 | 270.135 |
| x^{−1/50} | 1358.6 | 1358.658 | 1356.03 / **1357.53** | 1355.53 / **1357.43** | 1357.37 | 1354.81 |
| x^{−1/100} | 3757.6 | 3757.87 | unaffected (end lies in the FKS1 range) | | | |

\* The uniform perturbations overstate the effect for L ≲ 61. There the real increase is +0.04% to +0.6%, because T < H.

### 3.2 BKLNW Table 15 / Cor. 11.1 (`bklnw_t15.py`, `summary_calcs.py` §C)

- Table 10 shows the middle-range maxima sit at b = 2300–2350, at the switch from Büthe's method (ε(2300) = 1.9194e-12) to the Platt–Trudgian method (ε(2325) = 1.8750e-12).
- My Thm 1 optimum with the printed formula at b = 2300 is 1.91484e-12, slightly below BKLNW's value. With log(3c) it is 1.95873e-12; with the thesis E₂, 1.96278e-12.
- On BKLNW's 25-step grid, B_k(2300, 2325) becomes:
  - k = 1: 4.554e-9 (> M₁ = 4.4627e-9)
  - k = 2: 1.0588e-5 (> 1.0376e-5)
  - k = 3: 0.024618 (> 0.024334)
  - k = 4: 57.236 (> 57.184)
  - k = 5: 1.3307e5 (≤ 1.3441e5)
- Under the thesis E₂: 0.024669 and 57.354.
- On a step-1 grid, with PT values extrapolated log-linearly from the b = 2325–2400 rows, the maxima drop to 0.02399 (k = 3) and 55.3 (k = 4). So the printed M₃ and M₄ look recoverable. M₁ and M₂ are not recoverable with BKLNW's inputs: the fine-grid M₁ is ≈ 4.51e-9.
- With current inputs (FKS2 tables, H = 3e12, FKS1 beyond e^2073) all five are comfortably below the printed values.

### 3.3 Johnston–Leong–Tudzi (`jlt_iter.py`)

- Lemma A.5 takes B₄ = max over Table 8 intervals of ε(b, b′)·b′⁴ (+ b⁴e^{−b}). The maximum is at b = 2000: 1.92279e-12 × 2500⁴ = 75.109. The next largest is the b = 2500 PT row at 73.41.
- Corrected: 76.65 (log3c) or 76.81 (thesis).
- My implementation of Props. A.6–A.7 reproduces their row-1 cM (7.809722e-3 against 7.80973e-3) and row-7 cM (4.883451 against 4.88346).
- With corrected B₄:
  - row 1: 7.8415e-3
  - row 2: 6.1089e-2, so (A.5) = cM₂·385^{0.65} gives 2.9276 (printed 2.91890)
  - row 7, which is (A.6): 4.9143 (printed 4.88346)
- Thesis E₂: 2.9285 and 4.9175.
- cN scales like B^{1/(β+1)}, which is why the shifts are small.
- I have not traced the effect on Theorem A.2 (4.591, 7.397) and Theorem 1.2 (2260).

### 3.4 Bordignon–Starichkova

They use M₁(1 + 3/(2·1000)) ≤ 4.47e-9 for x > e^1000. Corrected, this is 4.561e-9 (log3c) or 4.570e-9 (thesis), when taken from BKLNW. It is a negligible secondary term.

## 4. Caveats and items not verified

- The E₃ zero sums for rows with b < 60 use a model: the first 10 zeros exactly, then the Riemann–von Mangoldt density. It is **not rigorous**, but reproduces FKS2 Table 2 rows (b = 44–75) to within 0.2% and Büthe's Table 1 rows to within 0.02%. For b ≥ 65 the E₃ term is negligible and the computation matches `crosscheck_buthe.delta_quad`.
- FKS2 Table 2 prints c = 23.93 at log x₀ = log 10¹⁹. The printed ε (1.9220e-8) reproduces with c = 28.93, so 23.93 is probably a typo. The value is unaffected either way.
- **Bhattacharjee's Lemma 2 claim** would also change Büthe's Prop. 3 itself. That would touch Büthe 2016 Thm 2 (via (7.6)) and Johnston 2022, Lemma 2.1, by about 15% in the tail term. Its effect on the constants 4.92 and 9.06 was **not** assessed: unverified.
- Not traced: LP's final Ψ(x, y) theorems, JLT's 2260, and Bordignon–Starichkova's final constants. In each case the affected input is either re-justifiable or enters with small weight.
- No full text obtained: Axler–Nicolas 2023, Nicolas 2021, Chalker 2016, Bërdëllima 2025, Ralaivaosaona–Razakarinoro 2026, Fan–Kobayashi–Molnar 2026, Hernández Cáceres–Rúa 2026. Dusart 2018 and the Fiori–Kadiri online bound calculator were not examined.

## 5. Search coverage

- **Semantic Scholar citations** (13 Oct 2026): BKLNW (59), Büthe 2016 (38), Johnston–Yang (33), FKS2 (17), Lichtman–Pomerance (12). JSON files are `cit_*.json`.
- **Full texts:** 72 arXiv texts were fetched (in `pdfs/`), plus Nguyen–Pomerance, Pomerance–Spicer and Bhattacharjee 2024.
- **How dependencies were found:** grep of citation contexts (`grep_ctx.py`, `grep_ctx2.py`, outputs `ctx*.txt`), then reading the relevant sections. PNT+ was checked through its repository tree, the IEANTN Lean files (`pnt/`) and a GitHub issue search.
- **Background-only citations (no numerical use found):** Bellotti–Keller–Trudgian; Mossinghoff–Trudgian–Yang; Chourasiya–Simonič; Creech et al.; Marques–Trojovský; Keliher–Lee; Bayless–Kinlaw–Lichtman; Brent–Platt–Trudgian ×2; Ernvall-Hytönen–Palojärvi; Staple; Büthe 2018 and the Mertens paper; Axler 2017 ×2 (Thm 2 only); Patel–Yang; Johnston 2024/2025 ×4; Chahal et al.; Harm; Starichkova; Dimitrov; Gorodetsky; Makdad–Sorenson; Bach–Sorenson; Fiori–Jaskari 2026; Kadiri 2005; and several off-topic citers.

## 6. Files

| File | What it does |
|---|---|
| `zeros_model.py` | Thm 1 with the zero-sum model; log-space. |
| `e2_variants.py` | Printed, log3c and thesis E₂, plus the optimizer. |
| `variants_table.py` | §0 table → `variants_table_output.txt`. |
| `jy_lemma23.py`, `summary_calcs.py` | §2.1, §2.2 and §3.2–3.4 → `jy_lemma23_output.txt`, `summary_calcs_output.txt`. |
| `lp_prop21.py`, `lp_e55_scan.py` | §2.2 → `lp_prop21_output.txt`. |
| `fks2_table7.py` (data, helpers), `fks2_table7_v2.py` | §3.1 → `fks2_table7_v2_output.txt`. This output is authoritative; running `fks2_table7.py` directly misapplies corrections beyond e^2073. |
| `fks2_table6.py` | Table 6 room → `fks2_table6_output.txt`. |
| `bklnw_t15.py` | §3.2 → `bklnw_t15_output.txt`. |
| `jlt_iter.py` | §3.3 → `jlt_iter_output.txt`. |
| `small_rows.py` | Early crude check (ℓ ≤ 1 E₃ bound); superseded by `zeros_model.py`. |
