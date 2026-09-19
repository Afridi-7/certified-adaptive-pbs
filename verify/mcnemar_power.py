"""Connor (1987) sample size for paired McNemar, two-sided alpha=0.05, power 0.8.
n = [z_{a/2} sqrt(psi) + z_b sqrt(psi - d^2)]^2 / d^2, psi = p10+p01 (discordance), d = p10-p01."""
from scipy.stats import norm, binom
import numpy as np
za, zb = norm.ppf(0.975), norm.ppf(0.8)
def n_connor(p10,p01):
    psi, d = p10+p01, p10-p01
    return (za*np.sqrt(psi)+zb*np.sqrt(psi-d**2))**2/d**2
print("p10 p01 -> n (episodes)")
for p10,p01 in [(0.10,0.03),(0.12,0.04),(0.08,0.02),(0.15,0.07),(0.10,0.04),(0.12,0.05),(0.20,0.10),(0.08,0.03)]:
    print(f"{p10:.2f} {p01:.2f} -> {np.ceil(n_connor(p10,p01)):.0f}")
# exact (conditional binomial) power check by simulation for one design
def exact_power(n,p10,p01,sims=20000,rng=np.random.default_rng(1)):
    from scipy.stats import binomtest
    hits=0
    cnt=rng.multinomial(n,[p10,p01,1-p10-p01],size=sims)
    for b,c,_ in cnt:
        if b+c==0: continue
        if binomtest(b,b+c,0.5).pvalue<0.05: hits+=1
    return hits/sims
for (p10,p01) in [(0.10,0.04),(0.15,0.07)]:
    n=int(np.ceil(n_connor(p10,p01)))
    print(f"exact-test power at n={n} for ({p10},{p01}): {exact_power(n,p10,p01,5000):.3f}")
