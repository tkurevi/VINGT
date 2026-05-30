"""
nodal_multi.py
==============

Nodal analysis for commingled (multi-layer) production wells.

The reference depth for the nodal solution is the CommingledReservoir's
`z_ref` attribute (by default = deepest layer's mid-perforation).
The VLP profile must be marched from the surface DOWN to this same
depth -- not to the well bottom.  This module sets up the well
geometry accordingly.

The IPR/VLP intersection itself is found by Brent's method on
f(q) = Pwf_VLP(q) - Pwf_IPR_commingled(q) = 0

just as in single-layer nodal analysis.  The CommingledReservoir
exposes the same Pwf_at_q / q_at_Pwf / AOF interface, so this is a
drop-in replacement.
"""

from __future__ import annotations
import csv
import numpy as np
from scipy.optimize import brentq

from pvt import (bar_to_Pa, Pa_to_bar, ppm_to_molality, GRAVITY,
                 fmt_q, m3h_to_ls, m3h_to_m3d)
from vlp import (
    WellGeometry, ThermalContext, FluidStream,
    march_VLP,
)
from ipr_multilayer import CommingledReservoir


# =====================================================================
# 1. VLP at the reference depth
# =====================================================================
def vlp_curve_at_depth(WHP_bar, well_to_ref, GWR_std, m_NaCl,
                        thermal, q_array_m3h, n_segments=50, pump=None):
    """VLP Pwf curve at the reference depth (well_to_ref.depth_TVD).

    `well_to_ref` is a WellGeometry whose depth_TVD equals z_ref.
    Optional `pump` dict (z_intake_m, dP_Pa) adds an ESP boost.
    """
    Pwfs = np.zeros_like(q_array_m3h)
    for i, q_m3h in enumerate(q_array_m3h):
        flu = FluidStream(q_brine_std_m3h=float(q_m3h),
                          GWR_std=GWR_std, m_NaCl=m_NaCl)
        try:
            r = march_VLP(bar_to_Pa(WHP_bar), well_to_ref, flu, thermal,
                          n_segments=n_segments, pump=pump)
            Pwfs[i] = Pa_to_bar(r['BHFP'])
        except Exception:
            Pwfs[i] = np.nan
    return q_array_m3h, Pwfs


# =====================================================================
# 2. IPR curve for the commingled reservoir
# =====================================================================
def ipr_curve_commingled(comm, n_points=30, q_max_m3h=None):
    """Return (q_m3h_array, Pwf_bar_array) for the commingled IPR at
    the reference depth."""
    AOF_m3h = comm.AOF() * 3600.0
    if q_max_m3h is None:
        q_max_m3h = AOF_m3h * 1.05
    q_arr = np.linspace(0.0, q_max_m3h, n_points)
    Pwf_arr = np.array([
        comm.Pwf_at_q(q / 3600.0) * 1e-5 for q in q_arr
    ])
    return q_arr, Pwf_arr


# =====================================================================
# 3. Solve commingled operating point
# =====================================================================
def solve_operating_point_multi(WHP_bar, well_to_ref, GWR_std, m_NaCl,
                                 thermal, comm, n_segments=50,
                                 q_min_m3h=1.0, q_max_m3h=None,
                                 tol_q=0.5, pump=None):
    """Find the commingled-reservoir natural operating point.

    Optional ESP boost via `pump` dict (z_intake_m, dP_Pa) - the VLP
    is marched WITH the pump boost included, so the result is the
    *pumped* operating point.

    Returns
    -------
    dict with keys:
      converged       True/False
      WHP_bar
      q_op_m3h        total surface rate at op point
      Pwf_op_bar      Pwf at the reference depth
      drawdown_bar    P_res (deepest layer) - Pwf_op
      layer_rates     [(name, q_m3h), ...] at the op point
      profile         VLP profile from march_VLP (None if failed)
      reservoir       the CommingledReservoir
    """
    AOF_m3h = comm.AOF() * 3600.0
    if q_max_m3h is None:
        q_max_m3h = AOF_m3h * 1.5

    def f(q_m3h):
        try:
            flu = FluidStream(q_brine_std_m3h=float(q_m3h),
                              GWR_std=GWR_std, m_NaCl=m_NaCl)
            r = march_VLP(bar_to_Pa(WHP_bar), well_to_ref, flu, thermal,
                          n_segments=n_segments, pump=pump)
            Pwf_VLP = Pa_to_bar(r['BHFP'])
            Pwf_IPR = comm.Pwf_at_q(q_m3h / 3600.0) * 1e-5
            return Pwf_VLP - Pwf_IPR
        except Exception:
            return 1e10

    # Bracket: at q_min flow VLP should be ABOVE IPR (or close);
    # at q_max IPR should be below VLP.  Find a sign change.
    f_lo = f(q_min_m3h)
    f_hi = f(q_max_m3h)
    n_widen = 0
    while f_lo * f_hi > 0 and n_widen < 5:
        q_max_m3h *= 1.5
        if q_max_m3h > AOF_m3h * 3.0:
            break
        f_hi = f(q_max_m3h)
        n_widen += 1

    if f_lo * f_hi > 0:
        # No natural intersection
        return dict(
            converged=False,
            WHP_bar=WHP_bar,
            q_op_m3h=0.0, Pwf_op_bar=np.nan,
            drawdown_bar=np.nan,
            layer_rates=[],
            profile=None,
            reservoir=comm,
            reason=("No flow: VLP > IPR everywhere "
                    if f_lo > 0 else "No flow: IPR > VLP everywhere"),
        )

    try:
        q_op = brentq(f, q_min_m3h, q_max_m3h, xtol=tol_q)
    except Exception as e:
        return dict(
            converged=False,
            WHP_bar=WHP_bar,
            q_op_m3h=0.0, Pwf_op_bar=np.nan,
            drawdown_bar=np.nan,
            layer_rates=[],
            profile=None,
            reservoir=comm,
            reason=f"Brent failed: {e}",
        )

    # Run final VLP at op rate to get profile
    flu = FluidStream(q_brine_std_m3h=float(q_op),
                      GWR_std=GWR_std, m_NaCl=m_NaCl)
    profile = march_VLP(bar_to_Pa(WHP_bar), well_to_ref, flu, thermal,
                        n_segments=n_segments, pump=pump)
    Pwf_op_Pa = profile['BHFP']
    Pwf_op_bar = Pa_to_bar(Pwf_op_Pa)

    # Layer-by-layer breakdown
    layer_rates = []
    for name, q_si in comm.layer_rates_at_Pwf(Pwf_op_Pa):
        layer_rates.append((name, q_si * 3600.0))   # m^3/h

    return dict(
        converged=True,
        WHP_bar=WHP_bar,
        q_op_m3h=q_op,
        Pwf_op_bar=Pwf_op_bar,
        drawdown_bar=comm.P_res * 1e-5 - Pwf_op_bar,
        layer_rates=layer_rates,
        profile=profile,
        reservoir=comm,
    )


def print_operating_point_multi(op):
    """Pretty-print the commingled operating-point dict."""
    print("-" * 66)
    print("COMMINGLED OPERATING POINT")
    print("-" * 66)
    if not op['converged']:
        print(f"  Did NOT converge: {op.get('reason', 'unknown')}")
        print("-" * 66)
        return
    p = op['profile']
    print(f"  q_total        = {fmt_q(op['q_op_m3h'])}  (standard)")
    print(f"  WHP            = {op['WHP_bar']:9.2f}  bar")
    print(f"  Pwf @ z_ref    = {op['Pwf_op_bar']:9.2f}  bar")
    print(f"  Reference dep. = {op['reservoir'].z_ref:9.1f}  m TVD")
    print(f"  P_res (deepest)= {op['reservoir'].P_res*1e-5:9.2f}  bar")
    print(f"  Drawdown       = {op['drawdown_bar']:9.2f}  bar")
    print(f"  Total dP_well  = {p['dP_total']*1e-5:9.2f}  bar  "
          f"(grav {p['dP_grav']/p['dP_total']*100:.1f} %, "
          f"fric {p['dP_fric']/p['dP_total']*100:.1f} %)")
    print(f"  T (BH/WH)      = {p['T_BH_fluid']-273.15:5.1f} / "
          f"{p['T_WH_fluid']-273.15:5.1f}  degC")
    print()
    print("  Layer-by-layer breakdown at op point:")
    total_q_m3h = sum(q for _, q in op['layer_rates'])
    for name, q_m3h in op['layer_rates']:
        frac = q_m3h / total_q_m3h * 100 if total_q_m3h > 0 else 0
        print(f"    {name:<14s}: {fmt_q(q_m3h):<28s} ({frac:5.1f} %)")
    print("-" * 66)


def export_summary_csv_multi(op, csv_path):
    """Export commingled operating-point summary to CSV."""
    if not op['converged']:
        with open(csv_path, 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['status', 'reason'])
            w.writerow(['NO_FLOW', op.get('reason', 'unknown')])
        return
    p = op['profile']
    res = op['reservoir']
    rows = [
        ('WHP_bar',        f"{op['WHP_bar']:.3f}"),
        ('q_total_ls',     f"{m3h_to_ls(op['q_op_m3h']):.3f}"),
        ('q_total_m3d',    f"{m3h_to_m3d(op['q_op_m3h']):.2f}"),
        ('q_total_m3h',    f"{op['q_op_m3h']:.3f}"),
        ('Pwf_ref_bar',    f"{op['Pwf_op_bar']:.3f}"),
        ('z_ref_m',        f"{res.z_ref:.1f}"),
        ('P_res_deepest_bar', f"{res.P_res*1e-5:.3f}"),
        ('drawdown_bar',   f"{op['drawdown_bar']:.3f}"),
        ('dP_total_bar',   f"{p['dP_total']*1e-5:.3f}"),
        ('dP_grav_bar',    f"{p['dP_grav']*1e-5:.3f}"),
        ('dP_fric_bar',    f"{p['dP_fric']*1e-5:.3f}"),
        ('T_BH_fluid_C',   f"{p['T_BH_fluid']-273.15:.2f}"),
        ('T_WH_fluid_C',   f"{p['T_WH_fluid']-273.15:.2f}"),
        ('n_layers',       f"{len(op['layer_rates'])}"),
    ]
    for name, q_m3h in op['layer_rates']:
        rows.append((f'q_{name}_ls',  f"{m3h_to_ls(q_m3h):.3f}"))
        rows.append((f'q_{name}_m3d', f"{m3h_to_m3d(q_m3h):.2f}"))
        rows.append((f'q_{name}_m3h', f"{q_m3h:.3f}"))

    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['key', 'value'])
        for k, v in rows:
            w.writerow([k, v])


def export_layer_breakdown_csv(comm, csv_path, n_points=40):
    """Export a Pwf-q sweep showing each layer's contribution."""
    AOF_m3h = comm.AOF() * 3600.0
    q_arr = np.linspace(0.0, AOF_m3h * 1.05, n_points)
    rows = []
    layer_names = [L.name for L in comm.layers]
    for q_m3h in q_arr:
        Pwf_Pa = comm.Pwf_at_q(q_m3h / 3600.0)
        per_layer = dict(comm.layer_rates_at_Pwf(Pwf_Pa))
        row = {
            'q_total_ls':  f"{m3h_to_ls(q_m3h):.3f}",
            'q_total_m3d': f"{m3h_to_m3d(q_m3h):.2f}",
            'q_total_m3h': f"{q_m3h:.3f}",
            'Pwf_ref_bar': f"{Pwf_Pa*1e-5:.3f}",
        }
        for name in layer_names:
            q_layer_m3h = per_layer.get(name, 0.0) * 3600.0
            row[f'q_{name}_ls']  = f"{m3h_to_ls(q_layer_m3h):.3f}"
            row[f'q_{name}_m3d'] = f"{m3h_to_m3d(q_layer_m3h):.2f}"
            row[f'q_{name}_m3h'] = f"{q_layer_m3h:.3f}"
        rows.append(row)

    fields = list(rows[0].keys())
    with open(csv_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)
