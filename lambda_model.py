"""Exact flat Lambda-LTB bang-time family on an Omega_m=.3 background.

H0=1, Omega_m=.3, Omega_Lambda=.7, 2E=0, M=Omega_m*r^3/2.
R(t,r)=r*a(t-eta*f(r)), with a the flat Lambda-CDM scale factor.
For eta>=0 and monotone f'<=0, R'>0 for every post-bang shell.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.interpolate import PchipInterpolator, CubicSpline
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analysis import monotone_shape

OM=.3
OL=.7
AGE=2*np.arcsinh(np.sqrt(OL/OM))/(3*np.sqrt(OL))
OUT=Path(__file__).parent/'results'
OUT.mkdir(exist_ok=True)
ZG=np.linspace(.1,1.5,15)


def expansion(tau):
    if tau<=0:raise ValueError('Past the bang')
    b=1.5*np.sqrt(OL)*tau
    a=(OM/OL)**(1/3)*np.sinh(b)**(2/3)
    H=np.sqrt(OM/a**3+OL)
    adot=a*H
    addot=-OM/(2*a*a)+OL*a
    adddot=adot*(OM/a**3+OL)
    return a,adot,addot,adddot


def geometry(t,r,eta,width):
    f,fp=monotone_shape(r,width)
    a,adot,addot,adddot=expansion(t-eta*f)
    rp=a-r*eta*fp*adot
    drp=adot-r*eta*fp*addot
    ddrp=addot-r*eta*fp*adddot
    return r*a,rp,drp,ddrp


def ray(eta,width,t0=AGE):
    def rhs(r,y):
        R,rp,drp,ddrp=geometry(y[0],r,eta,width)
        return [-rp,(1+y[1])*drp,drp*y[2]+ddrp]
    def stop(r,y):return y[1]-1.6
    stop.terminal=True;stop.direction=1
    sol=solve_ivp(rhs,(0,1.8),[t0,0.,0.],events=stop,rtol=2e-11,
                  atol=2e-13,max_step=.003,dense_output=True)
    if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('Ray integration failed')
    r=np.linspace(0,sol.t_events[0][0],3001)
    t,z,Z=sol.sol(r)
    geo=np.array([geometry(ti,ri,eta,width) for ti,ri in zip(t,r)])
    if np.any(geo[:,1]<=0) or np.any(np.diff(z)<=0):raise AssertionError('Invalid cone')
    D=(1+z)**2*geo[:,0]
    if np.any(np.diff(D)<=0):raise AssertionError('Nonmonotone distance')
    return {'r':r,'z':z,'Z':Z,'D':D,
            'at_r_z':CubicSpline(r,z),'at_r_Z':CubicSpline(r,Z),
            'at_z_Z':PchipInterpolator(z,Z),
            'at_z_D':PchipInterpolator(z,D),
            'at_D_Z':PchipInterpolator(D,Z)}


def main():
    base=ray(0,.8)
    flrw=(1+ZG)-np.sqrt(OM*(1+ZG)**3+OL)
    baseline_error=float(np.max(abs(base['at_z_Z'](ZG)-flrw)))
    if baseline_error>1e-8:raise AssertionError('Lambda-CDM limit failed')
    # Independent neighboring observer time check in the perturbed model.
    htime=1e-5
    test=ray(.001,.8)
    plus=ray(.001,.8,t0=AGE+htime)
    minus=ray(.001,.8,t0=AGE-htime)
    check_error=abs((plus['at_r_z'](.8)-minus['at_r_z'](.8))/(2*htime)-test['at_r_Z'](.8))
    if check_error>1e-7:raise AssertionError('Neighboring-ray check failed')
    out={'Omega_m':OM,'Omega_Lambda':OL,'H0_age':AGE,
         'baseline_max_abs_error':baseline_error,
         'neighbor_ray_check_abs_error':float(check_error),'cases':{}}
    fig,ax=plt.subplots(figsize=(6.5,3.8))
    for w in (.8,1.2):
        D0=base['at_z_D'](ZG)
        secants=[]
        for h in (.002,.001,.0005):
            p=ray(h,w)
            secants.append((p['at_D_Z'](D0)-base['at_D_Z'](D0))/h)
        first=2*secants[1]-secants[0]
        second=2*secants[2]-secants[1]
        err=float(max(abs(first-second)))
        out['cases'][str(w)]={'z':ZG.tolist(),'equal_distance_response':second.tolist(),
                              'max_halfstep_change':err,
                              'minimum':float(min(second)),
                              'maximum':float(max(second))}
        ax.plot(ZG,second,marker='o',ms=2.5,label=f'w={w}')
    ax.axhline(0,color='black',lw=.7)
    ax.set(xlabel='Baseline redshift z (at equal $D_L$)',
           ylabel=r'$\partial_\eta(\dot z/H_0)|_{D_L}$',
           title=r'Flat $\Lambda$-LTB bang-time response')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/'lambda_response.pdf')
    fig.savefig(OUT/'lambda_response.png',dpi=180)
    (OUT/'lambda.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='cases'}|
                     {'cases':{k:{kk:vv for kk,vv in v.items() if isinstance(vv,float)}
                               for k,v in out['cases'].items()}},indent=2))


if __name__=='__main__':main()
