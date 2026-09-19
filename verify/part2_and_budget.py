"""Numerical verification of proposal claims (Part II, alpha=0.01 budget, McNemar power)."""
import numpy as np
from math import comb, ceil
from scipy import stats
rng = np.random.default_rng(0)

print("== Part II, k=1: P(miss)=1-n/N, required n >= (1-a)N ==")
for N in [100, 500, 1000]:
    for a in [0.05, 0.01]:
        n = ceil((1-a)*N - 1e-9)
        # Monte Carlo check at that n
        sims = 200_000
        crit = rng.integers(0, N, sims)          # position of the critical field
        # uniform w/o replacement: prob. crit field in first n of a random permutation = n/N
        miss = np.mean(crit >= n)               # equivalent by exchangeability of permutation
        print(f"N={N:5d} a={a}: n*={n} ({n/N:.3f}N)  exact miss={1-n/N:.4f}  MC={miss:.4f}")

print("\n== k>1: P(all k missed)=C(N-k,n)/C(N,n); minimal n for alpha=0.05 ==")
def nmin(N,k,a):
    for n in range(N+1):
        if comb(N-k,n)/comb(N,n) <= a: return n
for N in [100,500,1000]:
    row=[]
    for k in [1,2,3,5,10,20,50]:
        if k<=N: row.append(f"k={k}:{nmin(N,k,0.05)}")
    print(f"N={N}: "+"  ".join(row))
# MC check for one cell
N,k,n=500,5,nmin(500,5,0.05); sims=100_000
miss=np.mean([not (rng.permutation(N)[:n] < k).any() for _ in range(sims)])
print(f"MC check N=500,k=5,n={n}: miss={miss:.4f} exact={comb(N-k,n)/comb(N,n):.4f}")
# approximation n ≈ N(1-a^{1/k})
print("approx N(1-a^(1/k)) for N=1000:",[round(1000*(1-0.05**(1/k))) for k in [1,2,3,5,10,20,50]])

print("\n== Minimax: any (possibly adaptive) strategy with no informative features ==")
print("Against uniformly random placement, any strategy examining n fields has P(miss)=1-n/N;")
print("so a randomised uniform sampler is minimax-optimal; deterministic orders have worst-case miss=1.")

print("\n== alpha budget: CRC-style finite-sample correction  n/(n+1)*R + B/(n+1) <= a ==")
print("Budget available for empirical risk: a - (B-a)/n  -> negative iff n < B/a - 1")
for a in [0.05, 0.01]:
    print(f"a={a}: need n >= {1/a-1:.0f} (B=1); budget at n=50: {a-(1-a)/50:+.4f}, n=100: {a-(1-a)/100:+.4f}, n=200: {a-(1-a)/200:+.4f}")
