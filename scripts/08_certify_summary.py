"""Summarise Part III runs: certification rate (Wilson CI), median cells to certify, and fixed-convention calls."""
import sys, numpy as np, pandas as pd
def wilson(x,n,z=1.96):
    p=x/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return max(0,c-h),min(1,c+h)
for path in sys.argv[1:]:
    d=pd.read_csv(path); th=d.theta.iloc[0]
    print(f"\n===== {path}  theta={th:.2f}  M={d.M.iloc[0]} cells  episodes/cell={d.groupby('p_target').size().iloc[0]} =====")
    meths=[m for m in ['uniform','bet_uniform','bet_prop','bet_prop_cv','prio_naive'] if f'{m}_cert' in d]
    rows=[]
    for p,g in d.groupby('p_target'):
        r=dict(p=p, side='AT' if abs(p-th)<1e-9 else ('above' if p>th else 'below'))
        for m in meths:
            c=int(g[f'{m}_cert'].sum()); n=len(g); lo,hi=wilson(c,n)
            cert=g[g[f'{m}_cert']==1]
            r[m]=f"{c/n:.3f}"+(f" [{lo:.2f},{hi:.2f}]" if r['side']!='above' else "")
            r[m+'_cells']=int(cert[f'{m}_cells'].median()) if len(cert) and r['side']=='above' else ''
        r['conv100']=round(g.conv100_call.mean(),3)
        rows.append(r)
    pd.set_option('display.width',260); pd.set_option('display.max_columns',40)
    print(pd.DataFrame(rows).to_string(index=False))
