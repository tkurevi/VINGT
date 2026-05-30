"""
nodal.py
========

Nodal analysis for vertical CO2-brine geothermal production wells.

The classical nodal node is at the bottom of the well (perforations).
At the operating point:
   - the IPR (reservoir) and the VLP (wellbore + WHP boundary)
     give the same Pwf at the same flow rate q.

This module:
  * Generates IPR and VLP curves for plotting
  * Solves for the natural operating point (Brent's method on q)
  * Exports the operating-point depth profile to CSV in a format
    that can be diffed against engineering-software output

All flow rates are reported in m^3/h of brine at standard conditions.
"""

from __future__ import annotations
import csv
import numpy as np
from scipy.optimize import brentq

from pvt import (bar_to_Pa, Pa_to_bar, ppm_to_molality,
                 fmt_q, m3h_to_ls, m3h_to_m3d)
from vlp import (
    WellGeometry, ThermalContext, FluidStream,
    march_VLP,
)
from ipr import ReservoirModel


# =====================================================================
# 1. IPR curve generation
# =====================================================================
def ipr_curve(reservoir, n_points=40, q_max_m3h=None):
    """Generate (q, Pwf) arrays for the IPR.

    Parameters
    ----------
    reservoir   : ReservoirModel
    n_points    : number of points
    q_max_m3h   : optional upper bound on q (m^3/h); if None, use AOF

    Returns
    -------
    q_m3h, Pwf_bar : arrays
    """
    AOF_m3h = reservoir.AOF() * 3600.0
    q_high_m3h = q_max_m3h if q_max_m3h is not None else AOF_m3h * 1.05
    q_arr = np.linspace(0.0, q_high_m3h, n_points)
    Pwf_arr = np.array([reservoir.Pwf_at_q(q / 3600.0) for q in q_arr])
    return q_arr, Pwf_arr * 1.0e-5


# =====================================================================
# 2. VLP curve generation (at fixed WHP)
# =====================================================================
def vlp_curve(WHP_bar, well, GWR_std, m_NaCl, thermal,
              q_array_m3h, n_segments=50):
    """Generate (q, Pwf) arrays for the VLP at a fixed wellhead pressure.

    Parameters
    ----------
    WHP_bar     : wellhead pressure (bar)
    well        : WellGeometry
    GWR_std     : standard-conditions gas-water ratio
    m_NaCl      : brine NaCl molality
    thermal     : ThermalContext
    q_array_m3h : array of brine rates (m^3/h)
    n_segments  : VLP grid resolution

    Returns
    -------
    q_m3h, Pwf_bar : arrays
    """
    Pwf_arr = np.zeros_like(q_array_m3h, dtype=float)
    for i, q in enumerate(q_array_m3h):
        if q <= 0.0:
            # Limit: at q -> 0 there is no friction, only hydrostatic.
            # Use a tiny rate to avoid div-by-zero and recover the
            # limiting hydrostatic value naturally.
            q_use = max(q, 1.0e-3)
        else:
            q_use = q
        flu = FluidStream(q_use, GWR_std=GWR_std, m_NaCl=m_NaCl)
        r = march_VLP(bar_to_Pa(WHP_bar), well, flu, thermal,
                      n_segments=n_segments)
        Pwf_arr[i] = r['BHFP']
    return q_array_m3h, Pwf_arr * 1.0e-5


# =====================================================================
# 3. Operating-point solver
# =====================================================================
def solve_operating_point(WHP_bar, well, GWR_std, m_NaCl, thermal,
                          reservoir, n_segments=50,
                          q_min_m3h=1.0, q_max_m3h=None,
                          tol_q=0.5):
    """Find the natural operating point at the bottomhole node.

    f(q) = Pwf_VLP(q) - Pwf_IPR(q)
    Solve f(q*) = 0 by Brent's method.

    Returns
    -------
    dict with keys
      q_op_m3h        operating brine rate, m^3/h at standard conds
      Pwf_op_bar      operating bottomhole flowing pressure, bar
      WHP_bar         wellhead pressure (input)
      profile         the full march_VLP() result dict at the op point
      reservoir       the ReservoirModel passed in
      drawdown_bar    P_res - Pwf_op
      converged       True/False
      diagnostics     dict with extra info
    """
    AOF_m3h = reservoir.AOF() * 3600.0
    if q_max_m3h is None:
        q_max_m3h = AOF_m3h * 1.5

    def f(q_m3h):
        try:
            flu = FluidStream(q_m3h, GWR_std=GWR_std, m_NaCl=m_NaCl)
            r = march_VLP(bar_to_Pa(WHP_bar), well, flu, thermal,
                          n_segments=n_segments)
            Pwf_VLP = r['BHFP']
            if not np.isfinite(Pwf_VLP) or Pwf_VLP > 5.0e8:
                # Friction blew up; treat as VLP >> IPR
                return 1.0e10
            Pwf_IPR = reservoir.Pwf_at_q(q_m3h / 3600.0)
            return Pwf_VLP - Pwf_IPR
        except Exception:
            # Any PVT or numerical failure => treat as non-physical rate
            return 1.0e10

    # Bracket the root.  At q_min, VLP gives ~hydrostatic-only (low Pwf)
    # vs IPR which gives ~P_res (high) -> f<0.  At q_high, VLP gets very
    # high (friction) vs IPR which drops -> f>0.  Sweep upward to find
    # a sign change.
    q_lo = q_min_m3h
    f_lo = f(q_lo)
    if f_lo > 0.0:
        # Even at minimum rate the well will not flow naturally
        return dict(
            q_op_m3h=0.0, Pwf_op_bar=Pa_to_bar(reservoir.P_res),
            WHP_bar=WHP_bar, profile=None, reservoir=reservoir,
            drawdown_bar=0.0, converged=False,
            diagnostics=dict(reason="No flow: VLP > IPR even at q_min",
                             f_lo=f_lo, q_lo=q_lo))

    q_hi = max(q_max_m3h, q_lo * 4.0)
    # Cap iterations to avoid runaway
    for _ in range(20):
        f_hi = f(q_hi)
        if f_hi > 0.0:
            break
        q_hi *= 1.5
    else:
        return dict(
            q_op_m3h=q_hi, Pwf_op_bar=0.0, WHP_bar=WHP_bar,
            profile=None, reservoir=reservoir,
            drawdown_bar=Pa_to_bar(reservoir.P_res),
            converged=False,
            diagnostics=dict(reason="Could not bracket operating point",
                             q_hi=q_hi, f_hi=f_hi))

    # Brent's method
    q_op = brentq(f, q_lo, q_hi, xtol=tol_q, maxiter=60)

    # Final profile at the operating point
    flu = FluidStream(q_op, GWR_std=GWR_std, m_NaCl=m_NaCl)
    profile = march_VLP(bar_to_Pa(WHP_bar), well, flu, thermal,
                        n_segments=n_segments)
    Pwf_op = profile['BHFP']
    drawdown_bar = Pa_to_bar(reservoir.P_res - Pwf_op)

    return dict(
        q_op_m3h=float(q_op),
        Pwf_op_bar=float(Pwf_op * 1.0e-5),
        WHP_bar=float(WHP_bar),
        profile=profile,
        reservoir=reservoir,
        drawdown_bar=float(drawdown_bar),
        converged=True,
        diagnostics=dict(q_lo=q_lo, q_hi=q_hi),
    )


# =====================================================================
# 4. Pretty printer
# =====================================================================
def print_operating_point(op):
    """Print a clean operating-point summary."""
    print("-" * 66)
    print("OPERATING POINT")
    print("-" * 66)
    if not op['converged']:
        print(f"  Did NOT converge: {op['diagnostics'].get('reason','?')}")
        return
    p = op['profile']
    res = op['reservoir']
    print(f"  q_brine        = {fmt_q(op['q_op_m3h'])}  (standard)")
    print(f"  WHP            = {op['WHP_bar']:9.2f}  bar")
    print(f"  Pwf            = {op['Pwf_op_bar']:9.2f}  bar")
    print(f"  P_res          = {res.P_res*1e-5:9.2f}  bar")
    print(f"  Drawdown       = {op['drawdown_bar']:9.2f}  bar")
    print(f"  Total dP_well  = {p['dP_total']*1e-5:9.2f}  bar  "
          f"(grav {p['dP_grav']/p['dP_total']*100:4.1f} %, "
          f"fric {p['dP_fric']/p['dP_total']*100:4.1f} %)")
    print(f"  T (BH/WH)      = {p['T_BH_fluid']-273.15:6.1f} / "
          f"{p['T_WH_fluid']-273.15:6.1f}  degC")
    if p['bubble_point_depth'] is not None:
        print(f"  Bubble point   ~ {p['bubble_point_depth']:7.0f}  m TVD")
    else:
        print(f"  Bubble point   :  not crossed in well "
              "(under-saturated throughout)")
    print(f"  Skin (total)   = {res.s_total:9.2f}  "
          f"(d={res.s_d:.2f}, p={res.s_p:.2f}, gp={res.s_gp:.2f})")
    print("-" * 66)


# =====================================================================
# 5. CSV exporters - for engineering-software comparison
# =====================================================================
def export_profile_csv(op_or_profile, csv_path):
    """Export the depth profile to CSV.

    Columns: depth_m, P_bar, T_C, Te_C, H_l, rho_m_kgm3, v_m_ms,
             dpdL_grav_mbarpm, dpdL_fric_mbarpm, regime
    """
    if isinstance(op_or_profile, dict) and 'profile' in op_or_profile:
        p = op_or_profile['profile']
    else:
        p = op_or_profile

    rows = []
    n = len(p['z'])
    for i in range(n):
        rows.append({
            'depth_m':           f"{p['z'][i]:.2f}",
            'P_bar':             f"{p['P'][i]*1e-5:.4f}",
            'T_C':               f"{p['T'][i]-273.15:.3f}",
            'Te_C':              f"{p['T_e'][i]-273.15:.3f}",
            'H_l':               f"{p['H_l'][i]:.4f}",
            'rho_m_kgm3':        f"{p['rho_m'][i]:.2f}",
            'v_m_ms':            f"{p['v_m'][i]:.3f}",
            'dpdL_grav_mbarpm':  f"{p['dpdL_grav'][i]*1e-2:.3f}",
            'dpdL_fric_mbarpm':  f"{p['dpdL_fric'][i]*1e-2:.3f}",
            'regime':            p['regime'][i],
        })
    fields = list(rows[0].keys())
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def export_summary_csv(op, csv_path):
    """Export the operating-point summary as a single-row CSV."""
    if not op['converged']:
        with open(csv_path, 'w', newline='') as f:
            f.write("converged,reason\n")
            f.write(f"FALSE,{op['diagnostics'].get('reason','?')}\n")
        return
    p = op['profile']
    res = op['reservoir']
    summary = {
        'WHP_bar':           f"{op['WHP_bar']:.3f}",
        'q_brine_ls':        f"{m3h_to_ls(op['q_op_m3h']):.3f}",
        'q_brine_m3d':       f"{m3h_to_m3d(op['q_op_m3h']):.2f}",
        'q_brine_m3h':       f"{op['q_op_m3h']:.3f}",
        'Pwf_bar':           f"{op['Pwf_op_bar']:.3f}",
        'P_res_bar':         f"{res.P_res*1e-5:.3f}",
        'drawdown_bar':      f"{op['drawdown_bar']:.3f}",
        'dP_total_bar':      f"{p['dP_total']*1e-5:.3f}",
        'dP_grav_bar':       f"{p['dP_grav']*1e-5:.3f}",
        'dP_fric_bar':       f"{p['dP_fric']*1e-5:.3f}",
        'T_BH_fluid_C':      f"{p['T_BH_fluid']-273.15:.2f}",
        'T_WH_fluid_C':      f"{p['T_WH_fluid']-273.15:.2f}",
        'bubble_point_m':    (f"{p['bubble_point_depth']:.0f}"
                              if p['bubble_point_depth'] is not None
                              else "NA"),
        's_total':           f"{res.s_total:.3f}",
        's_damage':          f"{res.s_d:.3f}",
        's_perforation':     f"{res.s_p:.3f}",
        's_gravel_pack':     f"{res.s_gp:.3f}",
        'completion':        res.completion,
    }
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(summary.keys()))
        w.writeheader()
        w.writerow(summary)


def export_curves_csv(q_ipr, Pwf_ipr, q_vlp, Pwf_vlp, csv_path):
    """Export IPR + VLP curves to CSV (different lengths supported by
    side-by-side stacking with NaN padding)."""
    n = max(len(q_ipr), len(q_vlp))
    rows = []
    for i in range(n):
        rows.append({
            'q_ipr_ls':   (f"{m3h_to_ls(q_ipr[i]):.3f}"  if i < len(q_ipr)
                           else ""),
            'q_ipr_m3d':  (f"{m3h_to_m3d(q_ipr[i]):.2f}" if i < len(q_ipr)
                           else ""),
            'q_ipr_m3h':  f"{q_ipr[i]:.3f}"  if i < len(q_ipr)  else "",
            'Pwf_ipr_bar':f"{Pwf_ipr[i]:.3f}"if i < len(Pwf_ipr)else "",
            'q_vlp_ls':   (f"{m3h_to_ls(q_vlp[i]):.3f}"  if i < len(q_vlp)
                           else ""),
            'q_vlp_m3d':  (f"{m3h_to_m3d(q_vlp[i]):.2f}" if i < len(q_vlp)
                           else ""),
            'q_vlp_m3h':  f"{q_vlp[i]:.3f}"  if i < len(q_vlp)  else "",
            'Pwf_vlp_bar':f"{Pwf_vlp[i]:.3f}"if i < len(Pwf_vlp)else "",
        })
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


# =====================================================================
# 6. SELF-TEST / VALIDATION
# =====================================================================
def _validate():
    print("=" * 78)
    print("NODAL MODULE - SELF-CONSISTENCY TESTS")
    print("=" * 78)

    # --- Build a representative geothermal well ---
    well = WellGeometry(depth_TVD=2500.0, tubing_ID=0.114,
                        tubing_OD=0.127)
    therm = ThermalContext(T_surface=288.15, geo_gradient=0.080,
                           U_overall=20.0, time_seconds=60 * 86400.0)
    res = ReservoirModel(
        k_md=80.0, h_m=20.0,
        P_res_bar=320.0, T_res_C=215.0,
        r_w=0.108, r_e=300.0, m_NaCl=0.713,
        completion='perforated',
        perf_length=0.30, perf_diameter=0.012,
        perf_SPF_per_m=13.12, perf_phasing=90,
        damage_radius=0.30, damage_k_ratio=0.5,
        crushed_k_ratio=0.4, crushed_thickness=0.0127)

    GWR = 10.0
    m_NaCl = ppm_to_molality(40_000)

    print("\n--- Reservoir summary ---")
    print(res.describe())

    # --- Test 1: solve operating point at a single WHP ---
    print("\n[Test 1] Operating point at WHP = 8 bar")
    op = solve_operating_point(8.0, well, GWR, m_NaCl, therm, res)
    print_operating_point(op)

    # --- Test 2: WHP sweep showing operating point shifts ---
    print("\n[Test 2] WHP sweep")
    print(f"   {'WHP (bar)':>10}  {'q_op (m^3/h)':>13}  "
          f"{'Pwf (bar)':>10}  {'Drawdown':>10}")
    for whp in [3, 5, 8, 12, 18, 25]:
        op_i = solve_operating_point(float(whp), well, GWR, m_NaCl,
                                      therm, res, n_segments=40)
        if op_i['converged']:
            print(f"   {whp:10.1f}  {op_i['q_op_m3h']:13.1f}  "
                  f"{op_i['Pwf_op_bar']:10.1f}  "
                  f"{op_i['drawdown_bar']:10.1f}")
        else:
            print(f"   {whp:10.1f}  -- well does not flow --")

    # --- Test 3: IPR and VLP curves at same WHP, verify intersection ---
    print("\n[Test 3] Curve intersection sanity (WHP = 8 bar)")
    AOF = res.AOF() * 3600.0
    q_ipr, Pwf_ipr = ipr_curve(res, n_points=20)
    q_vlp = np.linspace(50.0, AOF * 1.05, 20)
    q_vlp_out, Pwf_vlp = vlp_curve(8.0, well, GWR, m_NaCl, therm,
                                    q_vlp, n_segments=40)
    # Find IPR at the operating q
    q_op = op['q_op_m3h']
    Pwf_ipr_at_op = res.Pwf_at_q(q_op / 3600.0) * 1.0e-5
    Pwf_vlp_at_op = op['Pwf_op_bar']
    print(f"   At q* = {q_op:.1f} m^3/h:")
    print(f"     Pwf_IPR = {Pwf_ipr_at_op:.2f} bar")
    print(f"     Pwf_VLP = {Pwf_vlp_at_op:.2f} bar")
    print(f"     |delta| = {abs(Pwf_ipr_at_op - Pwf_vlp_at_op):.4f} bar  "
          f"(should be near zero)")

    # --- Test 4: skin sensitivity ---
    print("\n[Test 4] Skin sensitivity (WHP = 8 bar)")
    print(f"   {'s_total':>8}  {'q_op (m^3/h)':>13}  {'Pwf':>8}  "
          f"{'Drawdown':>10}")
    for s in [-2.0, 0.0, 2.0, 5.0, 10.0]:
        res_s = ReservoirModel(
            k_md=80.0, h_m=20.0,
            P_res_bar=320.0, T_res_C=215.0,
            r_w=0.108, r_e=300.0, m_NaCl=0.713,
            completion='open_hole', skin_total=float(s))
        op_s = solve_operating_point(8.0, well, GWR, m_NaCl,
                                      therm, res_s, n_segments=40)
        if op_s['converged']:
            print(f"   {s:8.1f}  {op_s['q_op_m3h']:13.1f}  "
                  f"{op_s['Pwf_op_bar']:8.1f}  "
                  f"{op_s['drawdown_bar']:10.1f}")
        else:
            print(f"   {s:8.1f}  --  well does not flow")

    # --- Test 5: CSV exports ---
    print("\n[Test 5] Writing CSV exports for op @ WHP=8 bar")
    export_profile_csv(op, '/tmp/profile.csv')
    export_summary_csv(op, '/tmp/summary.csv')
    export_curves_csv(q_ipr, Pwf_ipr, q_vlp_out, Pwf_vlp,
                       '/tmp/curves.csv')
    import os
    for f in ('/tmp/profile.csv', '/tmp/summary.csv', '/tmp/curves.csv'):
        sz = os.path.getsize(f)
        print(f"   wrote {f}  ({sz} bytes)")
    # show first 3 lines of the profile CSV
    with open('/tmp/profile.csv') as fh:
        head = [next(fh) for _ in range(4)]
    print("   profile.csv header + first 3 rows:")
    for line in head:
        print("     " + line.rstrip())

    print("\n" + "=" * 78)


if __name__ == "__main__":
    _validate()
