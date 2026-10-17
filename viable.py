"""Globally shell-crossing-free parabolic bang-time family for eta >= 0."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis import ray, linear_response, monotone_shape

OUT=Path(__file__).parent/'results'
OUT.mkdir(exist_ok=True)
ZG=np.linspace(.1,1.5,15)

def main():
    base=ray(0,.8,profile=monotone_shape)
    fig,ax=plt.subplots(figsize=(6.5,3.8))
    result={}
    for w in (.8,1.2):
        exact=linear_response(w,profile=monotone_shape)
        D0=base['at_z_D'](ZG)
        estimates=[]
        for h in (.002,.001,.0005):
            perturbed=ray(h,w,profile=monotone_shape)
            estimates.append((perturbed['at_D_drift'](D0)-base['at_D_drift'](D0))/h)
        # Second-order Richardson extrapolation of positive-eta secants.
        richardson=2*estimates[-1]-estimates[-2]
        err=float(np.max(np.abs(richardson-exact)))
        if err>1e-5:raise AssertionError('positive-amplitude derivative check failed')
        result[str(w)]={'z':ZG.tolist(),'fixed_D_variational':exact.tolist(),
                        'fixed_D_positive_step':richardson.tolist(),
                        'max_positive_step_error':err,'minimum':float(np.min(exact)),
                        'maximum':float(np.max(exact))}
        ax.plot(ZG,exact,marker='o',ms=2.5,label=f'w={w}')
    ax.axhline(0,color='black',lw=.7)
    ax.set(xlabel='Baseline redshift z (at equal $D_L$)',
           ylabel=r'$\partial_\eta\dot z|_{D_L}$',
           title=r'Monotone bang-time profiles ($\eta\geq0$)')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/'monotone_response.pdf')
    fig.savefig(OUT/'monotone_response.png',dpi=180)
    (OUT/'monotone.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:{f:v for f,v in x.items() if isinstance(v,float)}
                      for k,x in result.items()},indent=2))

if __name__=='__main__':main()
