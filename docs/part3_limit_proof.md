# Proposition 4: limit for certifying a threshold decision (RQ3)

Added 2026-11-19 (exploratory theory; not pre-registered). Numbers: `scripts/59_threshold_limit.py` →
`results/threshold_limit_v1/summary.json`; fast checks in `tests/test_limit.py`.

## Setting

A specimen has $M$ cells with labels $L\in\{0,1\}^M$ (1 = blast), $B=\sum_iL_i$ blasts and $p=B/M$. The examiner sees
side information $S$ (positions, images, classifier scores, anything) and has independent randomisation $U$. A
*procedure* verifies distinct cells $I_1,I_2,\dots$, where $I_t$ is a function of
$(S,U,I_1,L_{I_1},\dots,I_{t-1},L_{I_{t-1}})$; it stops at a stopping time $\tau$ of the filtration
$\mathcal F_t=\sigma(S,U,(I_s,L_{I_s})_{s\le t})$ and outputs a decision $D\in\{\text{certify},\text{not}\}$ that is
$\mathcal F_\tau$-measurable. It is *valid at level $\alpha$* for $H_0:p\le\theta$ if, for every fixed population
$(S,L)$ with $B\le\theta M$, $\Pr_U(D=\text{certify})\le\alpha$.

**Assumption 2 (no information about labels).** Under the model $\Pr_B$, $L$ is uniformly distributed over the
$\binom MB$ arrangements with $B$ blasts and independent of $(S,U)$: what the examiner sees does not tell blasts from
other cells. (Assumption 1 of Theorem 1 is the analogue for the absence problem.)

## Statement

Let $B_0=\lfloor\theta M\rfloor$ and $X_n\sim\mathrm{Hyp}(M,B,n)$ be the number of blasts in $n$ uniform draws without
replacement. Let $\beta^*_n(B)$ be the power at $B$ of the randomised level-$\alpha$ test of $B_0$ that rejects when
$X_n>c$ and with probability $\gamma$ when $X_n=c$, where $\Pr_{B_0}(X_n>c)+\gamma\Pr_{B_0}(X_n=c)=\alpha$.

**Proposition 4.** Under Assumption 2, for every procedure valid at level $\alpha$ for $H_0:p\le\theta$, every
$B>B_0$ and every $n\le M$,
$$\Pr_B(D=\text{certify},\ \tau\le n)\ \le\ \beta^*_n(B),$$
whatever the selection rule and the stopping rule. The bound is attained by the procedure that verifies $n$ cells
chosen uniformly and applies the randomised test. The same holds for certifying $p<\theta$ with
$B_0=\lceil\theta M\rceil$ and rejection for small $X_n$.

## Proof

*Step 1 (reduction).* Fix $t\le n$ and condition on $\mathcal F_{t-1}$ and on $I_t$, which is a function of the
history. Given $(S,U)$ and the labels observed so far, Assumption 2 makes the unobserved labels a uniformly random
arrangement of $B-X_{t-1}$ blasts among the $M-t+1$ unverified cells, independent of which of these cells the history
selects. Hence
$$\Pr_B\bigl(L_{I_t}=1\mid\mathcal F_{t-1},I_t\bigr)=\frac{B-X_{t-1}}{M-t+1}.$$
By induction, conditional on $(S,U)$ the verified labels $(L_{I_1},\dots,L_{I_n})$ have the law $Q_{B,n}$ of the
first $n$ labels of a uniformly random permutation of the population, which does not depend on $(S,U)$; the indices
themselves carry no information beyond $(S,U)$ and the labels. The probability of an observed label sequence
$\ell_{1:n}$ with $x$ blasts is $(B)_x(M-B)_{n-x}/(M)_n$ (falling factorials), so the likelihood ratio of
$\Pr_{B}$ to $\Pr_{B_0}$ for $(S,U,I_{1:n},L_{I_{1:n}})$ equals $h_B(x)/h_{B_0}(x)$, the ratio of hypergeometric
probabilities, a non-decreasing function of $x$ (monotone likelihood ratio; Lehmann and Romano, Testing Statistical
Hypotheses, 3rd ed., Springer 2005, Sec. 3.4).

*Step 2 (validity under random arrangement).* $\Pr_{B_0}$ averages over populations with $B_0\le\theta M$ blasts, for
each of which the procedure certifies with probability at most $\alpha$; hence $\Pr_{B_0}(D=\text{certify})\le\alpha$.

*Step 3 (Neyman–Pearson).* Extend the procedure after $\tau$ by uniform draws among the unverified cells (this
changes nothing it reports), so that $I_{1:n}$ is always defined. Let $\varphi=\mathbf 1\{D=\text{certify},\ \tau\le
n\}$. Because $\{\tau\le n\}$ and the decision on it are determined by the first $n$ verifications, $\varphi$ is a (possibly randomised) test based on
$(S,U,I_{1:n},L_{I_{1:n}})$, and $\mathbb E_{B_0}\varphi\le\Pr_{B_0}(D=\text{certify})\le\alpha$. For the simple null
$\Pr_{B_0}$ against the simple alternative $\Pr_B$, the Neyman–Pearson lemma gives, among all tests of size at most
$\alpha$ based on these observations, maximal power for the test that rejects when the likelihood ratio exceeds a
constant and randomises at equality; by Step 1 this is the stated test on $X_n$. Hence
$\mathbb E_B\varphi\le\beta^*_n(B)$. Attainment: the fixed-$n$ procedure is valid for every population with
$B\le B_0$ (its size is non-decreasing in $B$ by the monotone likelihood ratio) and its power is $\beta^*_n(B)$.
$\square$

## Fields

If whole fields are verified and the *fields*, not the cells, are exchangeable (the multiset of field compositions is
fixed and its arrangement is uniformly random and independent of what the examiner sees), Steps 1–3 hold with fields
as units: for every null field population $\Pi_0$ with $p\le\theta$, the probability of certifying within $n$
fields is at most the power of the most powerful level-$\alpha$ test of $\Pi_0$ against the actual field population
based on $n$ fields drawn uniformly without replacement. Each choice of $\Pi_0$ gives a valid bound. When every field
holds either no blast or $c$ blasts, taking $\Pi_0$ with $\lfloor\theta M/c\rfloor$ such fields gives a hypergeometric
bound in fields. In that family, the design effect of the arrangement is $(c-pm)/(1-p)$ for fields of $m$ cells
(Kish), and the verified cells needed for 80% power, relative to individually sampled cells, track it
(`results/threshold_limit_v1`, entries `fields`).

## What the proposition does and does not say

- It concerns the probability of certifying *within $n$ verified cells*; it holds for every $n$ simultaneously, but
  the fixed-$n$ test attains it only at its own $n$. Procedures that are valid at every stopping time, such as ALPHA,
  reached 0.29–0.88 of the bound at the Table II settings (`table2` entries).
- It assumes that the examiner cannot predict labels. Classifier-driven designs can do better through information in
  the scores (Section V-B), and lose the advantage when the scores lose information (Section V-C).
- It assumes perfect verification (as in Section III); with verifier error the bound concerns the verifier's count.
- KL relaxation: by data processing, $\mathrm{kl}(\Pr_B(\cdot),\alpha)\le\mathrm{KL}(\mathrm{Hyp}(M,B,n)\,\|\,
  \mathrm{Hyp}(M,B_0,n))$, but this is infinite as soon as $n>B_0$ (the sample can then contain more than $B_0$
  blasts), so the exact bound is the one to use.

## Checks

- Reduction: an adaptive policy that picks cells by uninformative features and past labels gives the hypergeometric
  count distribution (chi-square $p=0.53$, 20,000 runs).
- A valid sequential betting test (size 0.021–0.028) never certified by $n$ more often than the bound
  (largest excess $-0.22$ over four populations).
- At $M=2{,}000$ the bound's 80%-power ranges differ from those of the non-randomised exact test of
  `scripts/52_resolvable_margin.py` by at most 0.6 percentage points.
