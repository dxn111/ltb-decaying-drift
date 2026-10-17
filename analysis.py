"""Reproduce the dimensionless parabolic LTB calculations in paper.tex.

Units: c=1, t_obs=1, 2M(r)=4r^3/9, E=Lambda=0.
Run: python analysis.py. Results are written to results/.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator, CubicSpline
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path(__file__).parent / "results"
OUT.mkdir(exist_ok=True)
ZGRID = np.linspace(0.1, 1.5, 15)


def shape(r, width):
    x = r / width
    return r*r*np.exp(-x*x), 2*r*np.exp(-x*x)*(1-x*x)


def geometry(t, r, eps, width):
    f, fp = shape(r, width)
    tau = t-eps*f
    if tau <= 0:
        raise ValueError("Ray reached the bang surface")
    rp = tau**(2/3) - 2*r*eps*fp/(3*tau**(1/3))
    drp = 2/(3*tau**(1/3)) + 2*r*eps*fp/(9*tau**(4/3))
    ddrp = -2/(9*tau**(4/3)) - 8*r*eps*fp/(27*tau**(7/3))
    R = r*tau**(2/3)
    return R, rp, drp, ddrp


def linear_response(width):
    """Integrate Eq. (7) at epsilon=0 and convert to equal-distance response."""
    def rhs(r, u):
        q = 1-r/3
        z = q**-2-1
        Z = (2/3)*((1+z)-(1+z)**1.5)
        f, fp = shape(r, width)
        g = f+r*fp
        A_t, B, B_t = 2/(3*q), 2/(3*q), -2/(9*q**4)
        C_t = 8/(27*q**7)
        A_p, B_p, C_p = -2*g/(3*q), 2*g/(9*q**4), -8*g/(27*q**7)
        return [-A_t*u[0]-A_p,
                (1+z)*B_t*u[0]+B*u[1]+(1+z)*B_p,
                (B_t*Z+C_t)*u[0]+B*u[2]+B_p*Z+C_p]
    rgrid = 3*(1-(1+ZGRID)**-0.5)
    sol = solve_ivp(rhs, (0, max(rgrid)), [0., 0., 0.], dense_output=True,
                    rtol=2e-12, atol=2e-14, max_step=0.005)
    if not sol.success:
        raise RuntimeError(sol.message)
    out = []
    for r, (tp, zp, Zp) in zip(rgrid, sol.sol(rgrid).T):
        q = 1-r/3
        z = q**-2-1
        Z = (2/3)*((1+z)-(1+z)**1.5)
        f, _ = shape(r, width)
        R = r*q*q
        Rt = 2*r/(3*q)
        Rr = q*q
        tr = -Rr
        zr = (1+z)*2/(3*q)
        Zr = 2*Z/(3*q)-2/(9*q**4)
        Dp = (1+z)**2*(Rt*tp-2*r*f/(3*q))+2*(1+z)*R*zp
        Dr = 2*(1+z)*zr*R+(1+z)**2*(Rt*tr+Rr)
        out.append(Zp-Zr*Dp/Dr)
    return np.array(out)


def ray(eps, width, rtol=2e-11, atol=2e-13, t0=1.):
    # y=(emission time, redshift, drift per unit observer proper time)
    def rhs(r, y):
        R, rp, drp, ddrp = geometry(y[0], r, eps, width)
        return [-rp, (1+y[1])*drp, drp*y[2]+ddrp]
    def stop(r, y):
        return y[1]-1.6
    stop.terminal = True
    stop.direction = 1
    sol = solve_ivp(rhs, (0, 1.8), [t0, 0., 0.], events=stop,
                    dense_output=True, rtol=rtol, atol=atol, max_step=0.005)
    if not sol.success or len(sol.t_events[0]) != 1:
        raise RuntimeError("Light cone failed to reach z=1.6")
    r = np.linspace(0, sol.t_events[0][0], 3001)
    t, z, drift = sol.sol(r)
    R = np.array([geometry(ti, ri, eps, width)[0] for ti, ri in zip(t, r)])
    geo = np.array([geometry(ti, ri, eps, width) for ti, ri in zip(t, r)])
    rp = geo[:, 1]
    if np.any(rp <= 0) or np.any(np.diff(z) <= 0):
        raise RuntimeError("Shell crossing or non-monotone redshift on the ray")
    D = (1+z)**2 * R
    # Endpoint formula is an FLRW expression, shown only to quantify its LTB error.
    endpoint = (1+z)*(2/(3*t0)) - geo[:, 2]/rp
    if np.any(np.diff(D) <= 0):
        raise RuntimeError("Distance is not monotone on the sampled ray")
    return {"r": r, "z": z, "D": D, "drift": drift,
            "at_r_z": CubicSpline(r, z), "at_r_drift": CubicSpline(r, drift),
            "at_z_D": PchipInterpolator(z, D),
            "at_z_drift": PchipInterpolator(z, drift),
            "at_z_endpoint": PchipInterpolator(z, endpoint),
            "at_D_drift": PchipInterpolator(D, drift)}


def run():
    base = ray(0, 0.8)
    eds = (2/3)*((1+ZGRID)-(1+ZGRID)**1.5)
    baseline_error = float(np.max(np.abs(base["at_z_drift"](ZGRID)-eds)))
    if baseline_error > 1e-8:
        raise AssertionError("EdS drift limit failed")
    # Independent two-observer-time check of the propagated drift at fixed r.
    check = ray(0.001, 0.8)
    rcheck = 0.8
    dt0 = 1e-5
    ptime = ray(0.001, 0.8, t0=1+dt0)
    mtime = ray(0.001, 0.8, t0=1-dt0)
    ztime = (ptime["at_r_z"](rcheck)-mtime["at_r_z"](rcheck))/(2*dt0)
    drift_check_error = float(abs(ztime-check["at_r_drift"](rcheck)))
    if drift_check_error > 1e-7:
        raise AssertionError("Neighboring-ray drift check failed")
    # This is an intentionally incorrect FLRW endpoint shortcut in LTB.
    cases = {}
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "savefig.dpi": 180})
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    fig2, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for width in (0.8, 1.2):
        rows = []
        for h in (2e-3, 1e-3, 5e-4):
            plus, minus = ray(h, width), ray(-h, width)
            D0 = base["at_z_D"](ZGRID)
            fixed_z = (plus["at_z_drift"](ZGRID)-minus["at_z_drift"](ZGRID))/(2*h)
            fixed_D = (plus["at_D_drift"](D0)-minus["at_D_drift"](D0))/(2*h)
            dist_resp = (plus["at_z_D"](ZGRID)-minus["at_z_D"](ZGRID))/(2*h)
            endpoint_resp = (plus["at_z_endpoint"](ZGRID)-minus["at_z_endpoint"](ZGRID))/(2*h)
            rows.append((h, fixed_z, fixed_D, dist_resp, endpoint_resp))
        h, fixed_z, fixed_D, dist_resp, endpoint_resp = rows[-1]
        convergence = float(np.max(np.abs(rows[-1][2]-rows[-2][2])))
        exact_linear = linear_response(width)
        variational_error = float(np.max(np.abs(fixed_D-exact_linear)))
        if variational_error > 1e-5:
            raise AssertionError("Variational response check failed")
        cases[str(width)] = {
            "z": ZGRID.tolist(), "fixed_z": fixed_z.tolist(),
            "fixed_D": fixed_D.tolist(), "distance_at_fixed_z": dist_resp.tolist(),
            "incorrect_endpoint_response_at_fixed_z": endpoint_resp.tolist(),
            "endpoint_error_at_z_1": float((endpoint_resp-fixed_z)[9]),
            "fixed_D_minus_fixed_z_at_z_1": float((fixed_D-fixed_z)[9]),
            "variational_fixed_D": exact_linear.tolist(),
            "max_variational_difference": variational_error,
            "min_fixed_D": float(np.min(fixed_D)),
            "max_fixed_D": float(np.max(fixed_D)),
            "max_halfstep_change": convergence,
            "max_rel_distance_change_for_eps_0p001": float(np.max(np.abs(0.001*dist_resp/base["at_z_D"](ZGRID))))
        }
        ax.plot(ZGRID, fixed_D, marker="o", markersize=2.5, label=fr"$w={width}$")
        color = "#1764ac" if width == 0.8 else "#d46b21"
        axes[0].plot(ZGRID, fixed_z, color=color, lw=1.6, label=fr"$w={width}$, fixed $z$")
        axes[0].plot(ZGRID, fixed_D, color=color, ls="--", lw=1.6,
                     label=fr"$w={width}$, equal $D_L$")
        axes[1].plot(ZGRID, endpoint_resp-fixed_z, color=color, lw=1.6,
                     label=fr"$w={width}$")
    ax.axhline(0, color="black", lw=0.7)
    ax.set(xlabel=r"Baseline redshift $z$ (at fixed baseline $D_L$)",
           ylabel=r"$\partial_\epsilon\dot z\vert_{D_L}$ (dimensionless)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT/"fixed_distance_response.pdf")
    fig.savefig(OUT/"fixed_distance_response.png", dpi=180)
    axes[0].axhline(0, color="black", lw=0.6)
    axes[0].set(xlabel="Baseline redshift $z$", ylabel=r"$\partial_\epsilon\dot z$",
                title="Choice of comparison observable")
    axes[0].legend(fontsize=8, frameon=False)
    axes[1].axhline(0, color="black", lw=0.6)
    axes[1].set(xlabel="Baseline redshift $z$",
                ylabel=r"endpoint shortcut $-$ exact response",
                title="Error from an FLRW endpoint formula")
    axes[1].legend(fontsize=8, frameon=False)
    fig2.tight_layout()
    fig2.savefig(OUT/"comparison_and_endpoint_error.pdf")
    fig2.savefig(OUT/"comparison_and_endpoint_error.png", dpi=180)
    result = {"baseline_max_abs_error": baseline_error,
              "neighbor_ray_check_abs_error": drift_check_error, "cases": cases,
              "note": "Derivative changes t_B at fixed M and E, then compares at equal D_L; the entire D_L(z) curve is not compensated."}
    (OUT/"numbers.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"baseline_max_abs_error": baseline_error,
                      "neighbor_ray_check_abs_error": drift_check_error,
                      "cases": {w:{k:v for k,v in c.items() if not isinstance(v,list)}
                                for w,c in cases.items()}}, indent=2))


if __name__ == "__main__":
    run()
