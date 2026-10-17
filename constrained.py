"""Local first-order exact-distance inversion in the EdS mass-coordinate gauge.

Solves for curvature e(r) so that delta D_L(z)=0 along the baseline light cone.
The evolution is singular at the EdS angular-diameter distance maximum r=1;
solutions are evaluated only for 0 < r < 1 and regularity there is inspected.
"""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analysis import shape, monotone_shape

OUT=Path(__file__).parent/'results'
OUT.mkdir(exist_ok=True)


def solve(width, end=0.985, r0=1e-5, profile=shape):
    def local(r,u):
        tp,zp,Zp=u
        q=1-r/3
        z=q**-2-1
        Z=(2/3)*((1+z)-(1+z)**1.5)
        f,fp=profile(r,width)
        gd=f+r*fp
        e=-(40/27)*q**-5*(tp-f)-(10/3)*(r-1)/r*zp
        Eu1=-40/(27*q**5)
        Eu2=-(10/3)*(r-1)/r
        Er=-(200/81)*q**-6*(tp-f)+(40/27)*q**-5*fp-(10/3)*zp/r**2
        At,B,Bt,Ct=2/(3*q),2/(3*q),-2/(9*q**4),8/(27*q**7)
        Ap0=9*e*q**4/20-2*gd/(3*q)-r*r*e*q*q/2
        Bp0=3*e*q/5+2*gd/(9*q**4)-r*r*e/(3*q)
        F1=-At*tp-Ap0
        F2=(1+z)*Bt*tp+B*zp+(1+z)*Bp0
        G1=-9*r*q**4/20
        G2=3*r/(5*q)
        denominator=1-Eu1*G1-Eu2*G2
        ep=(Er+Eu1*F1+Eu2*F2)/denominator
        gg=e+r*ep
        Bp=3*gg*q/5+2*gd/(9*q**4)-r*r*e/(3*q)
        Cp=gg/(5*q*q)-8*gd/(27*q**7)+r*r*e/(9*q**4)
        F3=(Bt*Z+Ct)*tp+B*Zp+Bp*Z+Cp
        return np.array([F1+G1*ep,F2+G2*ep,F3]),e,ep,denominator
    sol=solve_ivp(lambda r,u:local(r,u)[0],(r0,end),[0.,0.,0.],
                  rtol=1e-10,atol=1e-12,max_step=.001,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    grid=np.linspace(r0,end,500)
    u=sol.sol(grid).T
    e=np.array([local(r,v)[1] for r,v in zip(grid,u)])
    ep=np.array([local(r,v)[2] for r,v in zip(grid,u)])
    z=(1-grid/3)**-2-1
    # Independent algebraic check of the equal-distance constraint.
    q=1-grid/3
    f=np.array([profile(r,width)[0] for r in grid])
    Q=2*grid/(3*q**5)*(u[:,0]-f)+1.5*(grid-1)*u[:,1]+9*grid*e/20
    Z=(2/3)*((1+z)-(1+z)**1.5)
    Zr=2*Z/(3*q)-2/(9*q**4)
    zr=2/(3*q**3)
    drift=u[:,2]-Zr*u[:,1]/zr
    return {"r":grid,"z":z,"e":e,"ep":ep,"drift":drift,
            "max_constraint_error":float(np.max(np.abs(Q)))}


def main():
    results={}
    fig,ax=plt.subplots(1,2,figsize=(10,3.7))
    for w in (.8,1.2):
        x=solve(w,end=.98,profile=monotone_shape)
        check=solve(w,end=.9,r0=1e-6,profile=monotone_shape)
        e_at_09=float(CubicSpline(x['r'],x['e'])(.9))
        e_at_09_check=float(check['e'][-1])
        if abs(e_at_09-e_at_09_check)>2e-3:
            raise AssertionError('center-start convergence failed')
        results[str(w)]={k:(v.tolist() if isinstance(v,np.ndarray) else v)
                         for k,v in x.items()}
        results[str(w)]['e_at_r_0p9']=e_at_09
        results[str(w)]['e_at_r_0p98']=float(x['e'][-1])
        results[str(w)]['center_start_difference_at_r_0p9']=abs(e_at_09-e_at_09_check)
        ax[0].plot(x['z'],x['e'],label=f'w={w}')
        ax[1].plot(x['z'],x['drift'],label=f'w={w}')
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set(xlabel='Baseline redshift z',ylabel='Required curvature e(r)',
              title='Curvature needed to match distance')
    ax[1].set(xlabel='Baseline redshift z',ylabel=r'$\partial_\eta\dot z|_{D_L}$',
              title='Local drift after compensation')
    for a in ax:a.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/'constrained_response.pdf')
    fig.savefig(OUT/'constrained_response.png',dpi=180)
    (OUT/'constrained.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({k:{'zmax':v['z'][-1],'e_end':v['e'][-1],
                         'drift_end':v['drift'][-1],
                         'max_constraint_error':v['max_constraint_error']}
                      for k,v in results.items()},indent=2))

if __name__=='__main__':main()
