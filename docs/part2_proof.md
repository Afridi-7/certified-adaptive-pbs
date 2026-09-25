# Part II — Rare-event certification: results and proofs (RQ2)

Status: draft v1 (2026-09-25). Numerical verification: `scripts/17_part2_verify.py` → `results/part2_v1/`
(all exact checks pass; see §7). This document gives the mathematics; the empirical check of Result 1 is H1d
(`docs/confirmatory_v1.md`, confirmed).

**Question (RQ2).** Why does certifying that rare, decision-changing evidence is *absent* require near-exhaustive
examination in the worst case, and when does this worst case stop describing the problem?

---

## 1. Setting and notation

- A specimen has $N$ fields, indexed $[N]=\{1,\dots,N\}$. (In our composites a field holds $m$ cells; a field is
  *critical* if it contains at least one blast. Everything below is at the field level.)
- $y_i\in\{0,1\}$ says whether field $i$ is critical; $C=\{i: y_i=1\}$ is the critical set and $k=|C|$.
- $x=(x_1,\dots,x_N)$ is everything the examiner can observe **before** verifying a field (e.g. an embedding or a
  low-resolution preview). It is treated as fixed and known.
- **Verification.** Inspecting field $i$ reveals $y_i$ without error (relaxed in Remark 4).
- **Examiner (policy) $\pi$.** At step $t=1,\dots,n$ it chooses a not-yet-inspected field $I_t$ as a function of
  $x$, the history $(I_s,y_{I_s})_{s<t}$, and independent internal randomness $U$. It may be adaptive and randomised,
  and may stop early, but inspects at most $n$ fields.
- **Miss event.** $M=\{C\cap\{I_1,\dots,I_n\}=\emptyset\}$: no critical field was inspected.
- **Probability.** $P_\pi(M\mid x,C)$ is over the examiner's randomness $U$ **only**. The specimen $(x,C)$ is fixed.
  All guarantees below are therefore *finite-population and design-based*: they hold for every specimen in the stated
  class, not on average over patients, and say nothing about a population of specimens or about diagnostic accuracy.
- $h(N,k,n) := \binom{N-k}{n}\big/\binom{N}{n}$ (defined as $0$ when $n>N-k$).
- $n^*(N,k,\alpha) := \min\{n: h(N,k,n)\le\alpha\}$.

**Assumption A (no location information).** For the given $x$, every placement $C\subseteq[N]$ with $|C|=k$ is
possible. In other words, the pre-verification observations do not restrict where the critical fields can be. The
worst case is taken over this class, $\mathcal S(x,k)=\{C:|C|=k\}$.

Assumption A is a **worst-case modelling assumption**, not an empirical claim about blood smears. §5 discusses when
it is violated, and what that does and does not buy.

---

## 2. Result 1 — one critical field, uniform sampling (equality)

**Proposition 1.** Let $\pi_U$ be simple random sampling without replacement (SRSWOR) of $n$ fields, and let $k=1$.
Then for every $C=\{c\}$:
$$P_{\pi_U}(M\mid C)=1-\frac{n}{N}\quad\text{(exact)},\qquad
P_{\pi_U}(M\mid C)\le\alpha \iff n\ge(1-\alpha)N,\qquad n^*(N,1,\alpha)=\lceil(1-\alpha)N\rceil .$$

*Proof.* Under SRSWOR every $n$-subset is equally likely. The number of $n$-subsets avoiding $c$ is
$\binom{N-1}{n}$, so $P(M)=\binom{N-1}{n}/\binom{N}{n}=\frac{(N-1)!\,n!\,(N-n)!}{n!\,(N-1-n)!\,N!}=\frac{N-n}{N}$.
Equivalently, by symmetry each field is inspected with probability $n/N$. The inequality follows by rearranging, and
$n^*$ is the smallest integer that satisfies it. $\square$

*Edge cases.* $n=0$ gives $1$ and $n=N$ gives $0$. For $\alpha=0.05$ and $N=200$: $n^*=190$. For $N=1000$: $950$.

---

## 3. Result 2 — $k$ critical fields, uniform sampling (equality and bounds)

**Proposition 2.** Under SRSWOR of $n$ fields and $|C|=k$, for $0\le n\le N-k$:
$$P_{\pi_U}(M\mid C)=h(N,k,n)=\frac{\binom{N-k}{n}}{\binom{N}{n}}=\frac{\binom{N-n}{k}}{\binom{N}{k}}
=\prod_{i=0}^{k-1}\frac{N-n-i}{N-i}\qquad\text{(exact).}$$
Moreover,
$$\Big(\frac{N-n-k+1}{N-k+1}\Big)^{k}\;\le\;h(N,k,n)\;\le\;\Big(1-\frac nN\Big)^{k}. \tag{B}$$

*Proof.* The first equality counts $n$-subsets avoiding all of $C$. For the second, both ratios equal
$\frac{(N-k)!\,(N-n)!}{N!\,(N-n-k)!}$. The product form follows by cancelling factorials. For (B), each factor is
$1-\frac{n}{N-i}$ with $0\le i\le k-1$, which lies between $1-\frac{n}{N-k+1}$ and $1-\frac nN$. $\square$

**Corollary 2 (sufficient sample size).** $n\ge N(1-\alpha^{1/k})$ implies $h\le\alpha$, by the upper bound in (B).
This closed form is **sufficient, not exact**. It never undershoots $n^*$, and it overshoots by at most 2 fields over
the grid in §7.

**Exact $n^*(N,k,\alpha)$ as a fraction of $N$** (from `results/part2_v1/nstar_table.csv`):

| $k$ | $N=200$, α=0.05 | $N=1000$, α=0.05 | $N=1000$, α=0.01 |
|---|---|---|---|
| 1 | 190 (0.950) | 950 (0.950) | 990 (0.990) |
| 2 | 155 (0.775) | 777 (0.777) | 900 (0.900) |
| 3 | 126 (0.630) | 631 (0.631) | 784 (0.784) |
| 5 | 90 (0.450) | 450 (0.450) | 601 (0.601) |
| 10 | 51 (0.255) | 258 (0.258) | 368 (0.368) |
| 20 | 27 (0.135) | 138 (0.138) | 204 (0.204) |
| 50 | 11 (0.055) | 57 (0.057) | 86 (0.086) |

The required *fraction* barely depends on $N$. It is governed by $k$ and $\alpha$, approximately $1-\alpha^{1/k}$.
Near-exhaustive examination is a property of small $k$, not of large specimens.

---

## 4. Result 3 — the lower bound holds for every strategy (minimax)

Results 1–2 describe uniform sampling. The thesis needs a stronger statement: under Assumption A, **no** strategy,
however adaptive or clever, can do better in the worst case.

**Theorem 3 (minimax lower bound).** Under Assumption A, for every policy $\pi$ that inspects at most $n$ fields:
$$\max_{C:\,|C|=k}P_\pi(M\mid x,C)\;\ge\;h(N,k,n).$$
SRSWOR attains this bound for every $C$, so it is minimax-optimal, and the minimax miss probability equals
$h(N,k,n)$ exactly. In particular:

1. To guarantee $P(M)\le\alpha$ for every placement, **any** policy needs $n\ge n^*(N,k,\alpha)$. For $k=1$ this
   means $n\ge\lceil(1-\alpha)N\rceil$.
2. Every **deterministic** policy, however adaptive, has worst-case miss probability $1$ whenever $n\le N-k$.

*Proof.*

*Step 1 (the all-negative path).* Fix $x$ and a realisation $u$ of the examiner's randomness. Consider the run in
which every inspected field turns out negative. Every quantity the policy conditions on is then determined by $u$, so
this run inspects a fixed sequence $a_1(u),\dots,a_{n'}(u)$ with $n'\le n$. Let $A(u)=\{a_1(u),\dots,a_{n'}(u)\}$.

*Step 2 (the miss event is determined by $A(u)$).* Take any $C$.
- If $C\cap A(u)=\emptyset$, every field on the all-negative path is negative under $C$. The actual run coincides with
  that path and misses.
- If $C\cap A(u)\neq\emptyset$, let $t$ be the first index with $a_t(u)\in C$. Up to step $t-1$ the actual
  observations are all negative, so the actual run coincides with the all-negative path through step $t$, inspects
  $a_t(u)\in C$, and does not miss.

Hence $P_\pi(M\mid x,C)=P_U\big(C\cap A(U)=\emptyset\big)$.

*Step 3 (average over placements).* Let $C$ be uniform over the $k$-subsets, independent of $U$. That this is allowed
is exactly Assumption A. For fixed $u$, the fraction of $k$-subsets disjoint from $A(u)$ is
$\binom{N-|A(u)|}{k}/\binom Nk\ge\binom{N-n}{k}/\binom Nk=h(N,k,n)$, using Proposition 2's identity and $|A(u)|\le n$.
Averaging over $U$ gives $\mathbb E_C\,P_\pi(M\mid x,C)\ge h(N,k,n)$.

*Step 4 (worst case $\ge$ average).* $\max_C P_\pi(M\mid x,C)\ge\mathbb E_C\,P_\pi(M\mid x,C)\ge h(N,k,n)$.
Proposition 2 shows SRSWOR achieves $h$ for every $C$.

*Item 2.* A deterministic policy has a single $A$ with $|A|\le n\le N-k$. The adversary places $C$ inside
$[N]\setminus A$, and the miss probability is $1$. $\square$

This is Yao's minimax principle specialised to a search problem. The randomised strategy that achieves the bound is
an *equaliser*: its miss probability does not depend on where the critical field is. That is also why the empirical
check H1d (uniform miss $=1-n/N$) is a check of the implementation, not of blood-smear data.

**Corollary 4 (unknown $k$: certifying absence).** Suppose the examiner declares "no critical evidence" whenever it
inspects no critical field, and we require $P(\text{declare absent}\mid C)\le\alpha$ for every $C\neq\emptyset$ with
$k$ unknown. Under Assumption A, every policy needs $n\ge\lceil(1-\alpha)N\rceil$.

*Proof.* If $C\subseteq C'$ then $M(C')\subseteq M(C)$: missing all of $C'$ implies missing all of $C$. This holds
pathwise, for any policy. So the supremum over non-empty $C$ is attained at a singleton, and Theorem 3 with $k=1$
applies. $\square$

This is the precise sense in which *absence* of rare evidence is expensive. The cost is not a property of uniform
sampling but a lower bound on every method, unless something beyond Assumption A is known.

**Remark 4 (imperfect verification).** Suppose verification detects an inspected critical field with probability
$1-\beta$, independently of the selection. Then, for $k=1$ under SRSWOR, $P(M)=1-(1-\beta)\,n/N$. At $n=N$ this is
$\beta$. So if $\beta>\alpha$, absence cannot be certified at level $\alpha$ even by a full census, and the certificate
can be no better than the verifier. This matters for the thesis because the "expert" in Part III is simulated by
ground truth, i.e. $\beta=0$.

---

## 5. When the worst case does not describe the problem

The bound is tight **only** under its assumptions. Each departure below changes the answer, and each corresponds to
part of the thesis.

| Departure | What changes | Where in the thesis |
|---|---|---|
| **Informative features** (Assumption A fails: critical fields look different before verification) | The adversary's class shrinks, so a relevance-ordered search can beat $h$. But the gain is a *model* property, not a guarantee: it holds only as long as the feature–label relationship holds. | Part I: KT1 relevance AUC 0.977 → 0.896 as dilution grows. RQ4: blast recall 0.906 → 0.340 under acquisition shift. |
| **Known $k\ge2$** (or a known lower bound on $k$) | Cost falls from $\approx(1-\alpha)N$ to $\approx N(1-\alpha^{1/k})$ (Prop. 2). | Rarity grid of Part I. |
| **Different estimand** (proportion vs threshold instead of presence) | Making "blast proportion $>\theta$" true requires $\ge\theta M$ blasts, so the adversary cannot hide the evidence in one field. Cost is then governed by the margin $\lvert p-\theta\rvert$, not by $k=1$. | Part III: 40–445 of 500 cells, depending on the margin. |
| **Imperfect verification** | The floor rises to $\beta$ (Remark 4). | Not tested; a stated limitation. |

**Connection to Part I (interpretation).** The confirmatory InfoDPP-PAC runs stopped at a median of 33 of 200 fields.
Under Assumption A, Theorem 3 implies that *any* policy inspecting 33 of 200 fields has worst-case miss
probability at least $1-33/200=0.835$ at $k=1$. InfoDPP-PAC's observed miss rates at $k=1$ (0.315, 0.681, 0.759,
0.805 at $m=1,5,10,20$ cells per field) are *below* this figure. That is possible only because Assumption A fails on
our composites: blasts are partly visible in the embedding. As dilution removes that visibility, the observed miss
rate approaches the worst-case figure.

Two caveats:
- These are averages over randomly placed critical fields, whereas the theorem bounds the worst case.
- The stop point varies by episode.

The comparison is therefore illustrative, not a test.

**What the theorem does not say.**
- It does not say that adaptive methods are useless. It says their advantage has to come from information about
  location, and that this advantage carries no guarantee unless that information is itself certified.
- It does not describe the $k>1$ or proportion settings as if they were the $k=1$ worst case (project rule §11).
- It gives no probability over patients or specimens.

---

## 6. Appendix — the α = 0.01 calibration claim in the proposal

The proposal states that α = 0.01 is not achievable at realistic calibration sizes because of a finite-sample
correction. This claim is reproduced **if** the correction is of conformal-risk-control type with loss bound $B=1$:
$$\frac{n}{n+1}\hat R_n+\frac{B}{n+1}\le\alpha\iff\hat R_n\le\alpha-\frac{1-\alpha}{n}\quad(B=1).$$
The usable budget is negative iff $n<(1-\alpha)/\alpha$. That is $n<99$ positive calibration episodes at
$\alpha=0.01$, and $n<19$ at $\alpha=0.05$ (verified in §7, check 7).

This is an inference about which correction the proposal means. The thesis text should name the procedure
explicitly. The claim is separate from Theorem 3 and concerns calibration, not search.

---

## 7. Numerical verification (`scripts/17_part2_verify.py`)

All checks except 6 are exact (rational arithmetic, or a linear program).

| # | Check | Scope | Result |
|---|---|---|---|
| 1 | $h(N,1,n)=1-n/N$; $n^*=\lceil(1-\alpha)N\rceil$ | $N\le300$; α ∈ {0.10, 0.05, 0.01} | pass |
| 2 | product form of $h$ and bounds (B) | all $k$, $n$, $N\le120$ | pass |
| 3 | $h$ non-increasing in $k$ | $N\le120$ | pass |
| 4 | closed form $\lceil N(1-\alpha^{1/k})\rceil\ge n^*$ | table grid | pass (overshoot ≤ 2) |
| 5 | **minimax game value by LP equals $h$**; every deterministic strategy has worst case 1 | all $N\le7$, $k$, $n$ | pass |
| 6 | Monte Carlo of SRSWOR, 200,000 reps per cell (seed 20260701) | 4 cells | $\lvert z\rvert\le1.58$ |
| 7 | CRC budget sign change | α = 0.05, 0.01 | $n=19$, $n=99$ |

Check 5 is the independent confirmation of Theorem 3. For small $N$ it solves the zero-sum game between an examiner
who mixes over all $n$-subsets and an adversary who picks $C$, and finds the value equals $h(N,k,n)$ exactly.
Without location information, adaptivity reduces to choosing a subset (Step 1 of the proof), so subsets are the full
strategy space here.

## 8. For the thesis text

- Present Proposition 1, Proposition 2, Theorem 3, Corollary 4 and Remark 4 as mathematical results that are
  proven here.
- Present §5 as interpretation, supported by the Part I and RQ4 data.
- Theorem 3 is a standard minimax (Yao) argument. Present it as an application to this problem, not as a new
  technique.
- Open item: a standard textbook reference for Yao's principle (Yao 1977, FOCS) should be verified before it is cited.
