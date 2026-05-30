"""
main_VinGT1.py
==============

Driver for the VinGT-1 (Vinkovci Geothermal-1) commingled production
nodal-analysis tool.  VinGT-1 is a vertical exploration geothermal well
in eastern Croatia with a single open-hole completion across
1504-2700 m TVD through breccia-conglomerate / tuffs / tuffitic
breccia ("podloga neogena").  The producing interval shows complex
crossflow that was confirmed by production logging (PLT) during the
84-hour second eruptive test in November 2025.

Five+ permeable zones contribute to flow, including two zones (L2 in
the 1820-1850 m cyan band and L4b at 2278-2295 m) that act as THIEF
zones during low-rate eruptive flow, taking injection from the
wellbore at the higher-pressure layers above and below.

Configurations:
  DEFAULT_CONFIG    -- 6-layer commingled, static-gradient P_res, no
                       Forchheimer, no ESP.  Use for shut-in / low-rate
                       behaviour and as the starting point for fitting.
  CALIBRATED_CONFIG -- Per-layer J fit to the PLT net contributions,
                       P_res shifts on thief layers, Forchheimer D on
                       L4a (the dominant deep producer) fit to the
                       65 Hz ESP stable point, two-segment wellbore
                       (4-1/2" tubing + 9-5/8" csg + 7" slotted liner)
                       with ESP at 737 m.

Run:
    python main_VinGT1.py             # default
    python main_VinGT1.py calibrated  # PLT + ESP-calibrated
"""

from __future__ import annotations
import os, sys, json, copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from pvt import (ppm_to_molality, bar_to_Pa, fmt_q, m3h_to_ls,
                 m3h_to_m3d, GRAVITY)
from vlp import WellGeometry, ThermalContext, FluidStream
from ipr_multilayer import Layer, CommingledReservoir
from nodal_multi import (
    solve_operating_point_multi, print_operating_point_multi,
    export_summary_csv_multi, export_layer_breakdown_csv,
    ipr_curve_commingled, vlp_curve_at_depth,
)


# ===================================================================
# Reference depth for the commingled IPR
# ===================================================================
# z_ref = 2080 m matches the dynamic-gauge depth used during the
# eruptive test ("ugradnja manometra @2080 m") and the directly
# measured static pressure P_static = 202.07 bar.  This is the
# cleanest pressure anchor in the well.  Two layers (L1-L3) sit
# above and three (L4a, L4b, L5) sit below this reference.
# ===================================================================
Z_REF_VINGT1 = 2080.0

# Wellbore brine density at the reference depth, used by
# CommingledReservoir to translate each layer's P_res to z_ref via
# hydrostatic.  At ~2080 m, ~109 degC, ~14000 ppm NaCl, ~200 bar:
# Batzle-Wang gives rho ~ 977 kg/m^3.  We use 970 kg/m^3 as a single
# value (slightly less to account for the warmer deeper section).
RHO_WB_VINGT1 = 970.0


# ===================================================================
# DEFAULT_CONFIG -- 6 layers, static-gradient P_res, no Forchheimer
# ===================================================================
DEFAULT_CONFIG = dict(
    case_name='VinGT-1 commingled (6-layer, static-gradient P_res, no Forchheimer)',

    # --------- Well geometry: eruptive (5" DP + 9-5/8" csg + 7" liner) ---------
    # During both eruptive tests an RTTS packer was set at 1438 m on
    # 5" NC50 X-95 19.5# drill pipe (ID 4.276" = 108.6 mm).  The 9-5/8"
    # 47# L-80 production csg runs 0-1504 m (ID 220.5 mm); below the
    # liner top at 1449 m the flow is inside the 7" 29# slotted
    # production liner (ID 154.8 mm) all the way to TD 2699 m.
    well=dict(
        segments=[
            (   0.0, 1438.0, 0.1086, 0.1270, 4.6e-5),  # 5"  DP
            (1438.0, 1499.0, 0.2205, 0.2445, 4.6e-5),  # 9-5/8" csg
            (1499.0, Z_REF_VINGT1, 0.1571, 0.1778, 4.6e-5),  # 7"  slotted liner
        ],
        pump_depth_m=None,   # natural flow / eruption
    ),

    # --------- Thermal context ---------
    # Surface T ~ 15 degC (annual mean).  Geothermal gradient calibrated
    # to the static-gradient survey: 21.45 degC @ 0 m -> 131.79 degC @
    # 2700 m  =>  0.041 K/m.
    thermal=dict(
        T_surface_C=15.0,
        geo_gradient_K_m=0.041,
        U_overall=20.0,
        t_prod_days=60.0,
    ),

    # --------- Fluid (commingled) ---------
    # Lab analysis on separator water: 14.1 g/L Cl-equivalent salinity;
    # density 1.009 g/cm^3 at 25 degC; pH 7.1.
    # GWR is operator-reported 0.40 sm^3/sm^3 from field separator metering
    # (lab PVT flash gives 0.124-0.143; reconciled by assuming the lab
    # bottle lost some gas in transit).  Gas composition is methane-
    # dominant (85 mol% CH4, 8 mol% CO2, 7 mol% N2) -- the existing
    # PVT module's CO2-brine flash will over-state solubility a bit
    # at depth, but at this low GWR the resulting VLP error is small.
    fluid=dict(
        GWR_std=0.40,
        NaCl_ppm=14000.0,
    ),

    # --------- Reservoir: 6 layers ---------
    # All in the open-hole "podloga neogena" interval (1504-2700 m).
    # Layer definitions follow the geologist-marked zones on the
    # composite log + the PLT-derived inflow profile.
    layers=[
        # ----- L1: Z01 area, 1660-1700 m, net producer -----
        # PLT net contribution at eruptive (Pwf~199 bar @ z_ref):
        # +27 m^3/d, dominated by the +104 producer at 1675-1688.
        dict(
            name='L1_Z01',
            top_depth_m=1660.0, bottom_depth_m=1700.0,
            h_net_m=20.0,           # ~50% of gross, log highlight
            k_md=8.0,               # initial guess, will fit on PLT
            P_res_bar=163.86,       # static gradient @ 1680 m
            T_res_C=95.5,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
        # ----- L2: Z02 area, 1820-1850 m, NET INJECTOR (thief) -----
        # Cyan-marked on the log (distinct lithology).  PLT shows
        # net -26 m^3/d at eruptive -> effective P_res is BELOW the
        # static-gradient reading; we start at the gradient value
        # and let the calibration shift it down.
        dict(
            name='L2_Z02_thief',
            top_depth_m=1820.0, bottom_depth_m=1850.0,
            h_net_m=15.0,
            k_md=5.0,
            P_res_bar=178.21,
            T_res_C=101.0,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
        # ----- L3: Z03 area, 1900-1985 m, net producer -----
        # Includes Zone 03 (1925-1945) + the strong +40 producer at
        # 1935-1943.  Net PLT contribution +44 m^3/d.
        dict(
            name='L3_Z03',
            top_depth_m=1900.0, bottom_depth_m=1985.0,
            h_net_m=30.0,
            k_md=10.0,
            P_res_bar=189.40,
            T_res_C=104.7,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
        # ----- L4a: deep-mid producer, 2230-2380 m -----
        # Lumps the two strong producing intervals (+60.6 at
        # 2230-2278 and +54.9 at 2291-2380) into one macro-layer.
        # Net PLT contribution from producers in 2200-2400: +115.5 m^3/d.
        # The dominant producer in the well and the main candidate for
        # Forchheimer D in the calibrated config.
        dict(
            name='L4a_deep_mid',
            top_depth_m=2230.0, bottom_depth_m=2380.0,
            h_net_m=35.0,
            k_md=15.0,
            P_res_bar=223.59,           # static gradient @ 2305 m
            T_res_C=117.7,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
        # ----- L4b: composite thief in 2210-2400 m -----
        # Lumps ALL injecting intervals in the deep-mid section:
        #   -22.3 m^3/d @ 2210-2220
        #   -26.5 m^3/d @ 2220-2230
        #   -68.4 m^3/d @ 2278-2292   <- dominant thief
        #    -6.6 m^3/d @ 2380-2400   (partial of -11 in 2380-2413)
        # Total -123.8 m^3/d.  Qp-weighted centroid depth = 2266 m.
        # The thief is geometrically interleaved with L4a (same broad
        # lithology, but a sub-set of the rock has lower effective
        # P_res due to prior local depletion).
        dict(
            name='L4b_thief',
            top_depth_m=2210.0, bottom_depth_m=2400.0,
            h_net_m=25.0,           # net of all thief intervals
            k_md=8.0,
            P_res_bar=221.04,       # static gradient @ z_mid=2266 m
            T_res_C=115.8,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
        # ----- L5: Z_st / deepest zone, 2500-2660 m -----
        # Includes the pink-marked Zone st (2625-2645) and the
        # strong producers at 2515-2540.  Net +33 m^3/d.
        dict(
            name='L5_Zst_deep',
            top_depth_m=2500.0, bottom_depth_m=2660.0,
            h_net_m=50.0,
            k_md=5.0,
            P_res_bar=249.99,
            T_res_C=125.6,
            r_w=0.108,  r_e=300.0,
            NaCl_ppm=14000.0,
            regime='pss', skin_total=0.0,
        ),
    ],

    # --------- Commingled options ---------
    commingled=dict(
        reference_depth_m=Z_REF_VINGT1,
        wellbore_density_kg_m3=RHO_WB_VINGT1,
    ),

    # --------- Operating conditions ---------
    # Eruptive operating point (the PLT condition): WHP ~ 0.22 bar,
    # surface rate ~ 0.77 l/s, Pwf @ z_ref ~ 199 bar.
    operating=dict(
        WHP_bar=0.22,
        n_segments=50,
        pump=None,
    ),
)


# ===================================================================
# Field calibration anchors (locked in from the report)
# ===================================================================
PLT_NET_QP_M3D = {
    'L1_Z01':       +27.0,
    'L2_Z02_thief': -26.0,
    'L3_Z03':       +44.0,
    'L4a_deep_mid': +115.5,    # sum of +60.6 and +54.9 producers
    'L4b_thief':    -123.8,    # sum of all 2200-2400 thieves
    'L5_Zst_deep':  +33.0,
}
PLT_PWF_AT_ZREF_BAR = 199.0      # central estimate (range 198.5-200.8)
PLT_SURFACE_Q_M3D   = 66.4       # net surface rate during PLT

# 4 ESP stable points at PSD = 737 m (Nov 5-6 2025)
ESP_STABLE_POINTS = [
    dict(Hz=50, q_lps=9.0,  WHP_bar=1.2, P_intake_bar=58.4, T_WH_C=95.0),
    dict(Hz=55, q_lps=10.6, WHP_bar=2.3, P_intake_bar=24.3, T_WH_C=97.0),
    dict(Hz=60, q_lps=11.4, WHP_bar=2.5, P_intake_bar=17.1, T_WH_C=99.0),
    dict(Hz=65, q_lps=15.0, WHP_bar=2.9, P_intake_bar=11.1, T_WH_C=100.0),
]

# PROSPER's linear-fit summary (for cross-check, NOT used as a hard anchor)
PROSPER_PI_M3D_BAR = 19.1        # = 0.221 l/s/bar
PROSPER_AOF_M3D    = 3837.0      # = 44.4 l/s
PROSPER_PRES_BAR   = 202.07      # @ 2080 m, matches static gauge


# ===================================================================
# Helpers
# ===================================================================
def build_commingled(config):
    layers = [Layer(**L) for L in config['layers']]
    return CommingledReservoir(
        layers,
        reference_depth_m=config['commingled'].get('reference_depth_m'),
        wellbore_density_kg_m3=config['commingled'].get('wellbore_density_kg_m3'),
    )


def build_well_to_ref(config, z_ref):
    w = config['well']
    if 'segments' in w and w['segments'] is not None:
        # Replace the deep boundary of the last segment with z_ref if needed
        segs = [list(s) for s in w['segments']]
        if abs(segs[-1][1] - z_ref) > 1e-6:
            segs[-1][1] = float(z_ref)
        return WellGeometry(
            depth_TVD=float(z_ref),
            segments=[tuple(s) for s in segs],
            pump_depth_m=w.get('pump_depth_m'),
        )
    return WellGeometry(
        depth_TVD=float(z_ref),
        tubing_ID=w['tubing_ID'],
        tubing_OD=w.get('tubing_OD'),
        wellbore_dia=w.get('wellbore_dia'),
        roughness=w.get('roughness', 46e-6),
        pump_depth_m=w.get('pump_depth_m'),
        casing_ID=w.get('casing_ID'),
        casing_OD=w.get('casing_OD'),
        casing_roughness=w.get('casing_roughness'),
    )


def build_thermal(config, z_ref):
    t = config['thermal']
    return ThermalContext(
        T_surface=t['T_surface_C'] + 273.15,
        geo_gradient=t['geo_gradient_K_m'],
        U_overall=t.get('U_overall', 20.0),
        time_seconds=t.get('t_prod_days', 60.0) * 86400.0,
    )


# ===================================================================
# CALIBRATED_CONFIG  -- PLT + ESP fitted, with Forchheimer D
# ===================================================================
# Calibration sequence (see plt_calibration.py and forchheimer_calibration.py):
#   1. Fit per-layer k_md for the 4 producing layers + per-layer P_res
#      shift for the 2 thief layers, holding Pwf=199 bar @ z_ref during
#      the PLT eruptive flow.  (plt_calibrated_config.json)
#   2. Fit a single uniform Forchheimer D = 7378 (m^3/s)^-1 to the
#      55/60/65 Hz ESP stable points.  RMS residual 1.3 l/s.
#      (calibrated_config.json)
#   3. ESP dP_pump fit by nodal-solver root-find on the 65 Hz point
#      (q=15 l/s, WHP=2.9 bar, T_WH=100 degC, P_intake=11.1 bar).
#
FORCHHEIMER_D = 7377.8

CALIBRATED_CONFIG = copy.deepcopy(DEFAULT_CONFIG)
CALIBRATED_CONFIG['case_name'] = ('VinGT-1 commingled (PLT+ESP calibrated, '
                                   f'uniform D={FORCHHEIMER_D:.0f})')

# Replace layer parameters with calibrated values (from plt_calibration.py
# + uniform Forchheimer D from forchheimer_calibration.py)
_CALIBRATED_OVERRIDES = {
    'L1_Z01':       dict(k_md=19.189, P_res_bar=163.86, D_nonDarcy=FORCHHEIMER_D),
    'L2_Z02_thief': dict(k_md= 5.000, P_res_bar=162.05, D_nonDarcy=FORCHHEIMER_D),
    'L3_Z03':       dict(k_md=15.978, P_res_bar=189.40, D_nonDarcy=FORCHHEIMER_D),
    'L4a_deep_mid': dict(k_md=35.119, P_res_bar=223.59, D_nonDarcy=FORCHHEIMER_D),
    'L4b_thief':    dict(k_md= 8.000, P_res_bar=199.13, D_nonDarcy=FORCHHEIMER_D),
    'L5_Zst_deep':  dict(k_md= 6.143, P_res_bar=249.99, D_nonDarcy=FORCHHEIMER_D),
}
for L in CALIBRATED_CONFIG['layers']:
    if L['name'] in _CALIBRATED_OVERRIDES:
        L.update(_CALIBRATED_OVERRIDES[L['name']])

# Switch wellbore to the ESP installation: 4-1/2" tubing (ID 100.5 mm)
# from surface down to ESP packer/intake at 737 m; 9-5/8" production
# csg below to 1499 m; 7" slotted liner from 1499 m to z_ref.  This
# is the geometry used during the November 2025 ESP test.
ESP_PSD_M = 737.0
CALIBRATED_CONFIG['well'] = dict(
    segments=[
        (   0.0,  ESP_PSD_M, 0.1005, 0.1143, 4.6e-5),   # 4-1/2" tubing
        (ESP_PSD_M, 1499.0,  0.2205, 0.2445, 4.6e-5),   # 9-5/8" csg
        (1499.0,   Z_REF_VINGT1, 0.1571, 0.1778, 4.6e-5),  # 7" slotted liner
    ],
    pump_depth_m=ESP_PSD_M,
)

# Set the calibrated 65 Hz operating point.  dP_pump in bar is the
# value that brings the nodal solver to q_op = 15 l/s at WHP=2.9 bar.
# Determined empirically by fit_pump_dP_for_qtarget() in run_calibrated.py;
# good starting guess is ~80-90 bar.
CALIBRATED_CONFIG['operating'] = dict(
    WHP_bar=2.9,
    n_segments=80,
    pump=dict(
        z_intake_m=ESP_PSD_M,
        dP_bar=85.0,           # initial guess; refit by run_calibrated.py
    ),
)
# Production thermal profile uses the same surface T but a longer
# production history during the ESP test.
CALIBRATED_CONFIG['thermal'] = dict(
    T_surface_C=15.0,
    geo_gradient_K_m=0.041,
    U_overall=20.0,
    t_prod_days=2.0,           # ESP test was ~2 days at 65 Hz stable
)


# ===================================================================
# Plotting helpers (per VGGT-1 template)
# ===================================================================
def plot_commingled_ipr(comm, save_path, n_points=40, title=None):
    """IPR plot showing total + per-layer contributions (production
    direction only -- standard nodal plot).  For the signed/crossflow
    view, see crossflow_analysis."""
    q_arr, Pwf_arr = ipr_curve_commingled(comm, n_points=n_points)
    q_ls = m3h_to_ls(q_arr)
    fig, ax = plt.subplots(figsize=(9, 7))
    ax.plot(q_ls, Pwf_arr, 'k-', linewidth=2.5, label='Total IPR (production)')
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C6', 'C7']
    layer_names = [L.name for L in comm.layers]
    per_layer_q = {n: [] for n in layer_names}
    for q_m3h in q_arr:
        Pwf_Pa = comm.Pwf_at_q(q_m3h / 3600.0)
        for name, q_si in comm.layer_rates_at_Pwf(Pwf_Pa):
            per_layer_q[name].append(q_si * 3600.0)
    for i, name in enumerate(layer_names):
        q_layer = np.array(per_layer_q[name])
        ax.plot(m3h_to_ls(q_layer), Pwf_arr, '--',
                color=colors[i % len(colors)], linewidth=1.5,
                label=f"{name} (k={comm.layers[i].k_md:.1f} mD)")
    # Mark P_res for each layer at z_ref
    for i, L in enumerate(comm.layers):
        dP_h = comm._dP_hydro_to_layer(L)
        Pres_at_ref = (L.P_res + dP_h) * 1e-5
        ax.axhline(Pres_at_ref, color=colors[i % len(colors)],
                    linestyle=':', alpha=0.55,
                    label=f"P_res({L.name}) @ z_ref = {Pres_at_ref:.1f} bar")
    AOF_ls = m3h_to_ls(comm.AOF() * 3600.0)
    ax.scatter([AOF_ls], [0.0], color='red', s=80, marker='v', zorder=5,
               label=f"AOF = {AOF_ls:.1f} l/s "
                     f"({m3h_to_m3d(comm.AOF()*3600.0):.0f} m^3/d)")
    ax.set_xlabel('q  (l/s)')
    ax.set_ylabel(f'Pwf  (bar)  at z_ref = {comm.z_ref:.0f} m TVD')
    ax.set_title(title or 'Commingled IPR with per-layer breakdown')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc='best', fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_commingled_nodal(comm, WHP_bar, well_to_ref, GWR, m_NaCl, thermal,
                          op, save_path, n_points=25, n_segments=50,
                          title=None, pump=None):
    """Nodal plot at z_ref: total IPR + VLP + operating point."""
    AOF_m3h = comm.AOF() * 3600.0
    q_max = max(AOF_m3h * 1.05, op.get('q_op_m3h', 0) * 2.5, 30.0)

    q_ipr, Pwf_ipr = ipr_curve_commingled(comm, n_points=30,
                                            q_max_m3h=q_max)
    q_vlp_in = np.linspace(max(1.0, q_max * 0.02), q_max, n_points)
    q_vlp, Pwf_vlp = vlp_curve_at_depth(
        WHP_bar, well_to_ref, GWR, m_NaCl, thermal,
        q_vlp_in, n_segments=n_segments, pump=pump)

    fig, ax = plt.subplots(figsize=(9, 7))
    ax.plot(m3h_to_ls(q_ipr), Pwf_ipr, 'C2-',
            label='Commingled IPR', linewidth=2)
    ax.plot(m3h_to_ls(q_vlp), Pwf_vlp, 'C0-',
            label=f'VLP @ WHP={WHP_bar:.1f} bar', linewidth=2)
    if op['converged']:
        q_ls  = m3h_to_ls(op['q_op_m3h'])
        q_m3d = m3h_to_m3d(op['q_op_m3h'])
        ax.scatter([q_ls], [op['Pwf_op_bar']],
                   color='red', s=80, zorder=5, marker='o',
                   label=f"Op pt: {q_ls:.2f} l/s ({q_m3d:.0f} m^3/d), "
                         f"Pwf={op['Pwf_op_bar']:.1f} bar")
    ax.set_xlabel('q  (l/s)')
    ax.set_ylabel(f'Pwf  (bar)  at z_ref = {comm.z_ref:.0f} m')
    ax.set_title(title or 'Nodal analysis (commingled)')
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


def plot_layer_split_vs_q(comm, save_path, n_points=40, title=None):
    """Stacked-area-style plot: for each total q, show layer split."""
    AOF_m3h = comm.AOF() * 3600.0
    q_arr = np.linspace(0.0, AOF_m3h * 0.99, n_points)
    layer_names = [L.name for L in comm.layers]
    per_layer_q = {n: [] for n in layer_names}
    for q_m3h in q_arr:
        Pwf_Pa = comm.Pwf_at_q(q_m3h / 3600.0)
        for name, q_si in comm.layer_rates_at_Pwf(Pwf_Pa):
            per_layer_q[name].append(q_si * 3600.0)
    q_ls = m3h_to_ls(q_arr)
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 9), sharex=True)
    for i, name in enumerate(layer_names):
        axes[0].plot(q_ls, m3h_to_ls(np.array(per_layer_q[name])),
                     'o-', label=name, markersize=4)
    axes[0].plot(q_ls, q_ls, 'k--', alpha=0.5, label='Total (identity)')
    axes[0].set_ylabel('q_layer  (l/s)')
    axes[0].set_title(title or 'Per-layer contribution vs total rate')
    axes[0].legend(loc='best', fontsize=9)
    axes[0].grid(True, alpha=0.3)
    totals = np.array([sum(per_layer_q[n][j] for n in layer_names)
                        for j in range(n_points)])
    totals[totals == 0] = np.nan
    for i, name in enumerate(layer_names):
        frac = np.array(per_layer_q[name]) / totals * 100
        axes[1].plot(q_ls, frac, 'o-', label=name, markersize=4)
    axes[1].set_xlabel('q_total  (l/s)')
    axes[1].set_ylabel('Fraction of total flow  (%)')
    axes[1].set_ylim(0, 100)
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc='best', fontsize=9)
    fig.tight_layout()
    fig.savefig(save_path)
    plt.close(fig)


# ===================================================================
# Driver
# ===================================================================
def run_full_analysis(config, out_dir, verbose=True):
    os.makedirs(out_dir, exist_ok=True)
    print("=" * 72)
    print(f"CASE: {config.get('case_name', '')}")
    print("=" * 72)

    comm = build_commingled(config)
    well = build_well_to_ref(config, comm.z_ref)
    thermal = build_thermal(config, comm.z_ref)
    WHP    = config['operating']['WHP_bar']
    n_seg  = config['operating'].get('n_segments', 50)
    GWR    = config['fluid']['GWR_std']
    m_NaCl = ppm_to_molality(config['fluid']['NaCl_ppm'])

    pump_cfg = config['operating'].get('pump')
    pump = None
    if pump_cfg is not None:
        pump = dict(
            z_intake_m=float(pump_cfg['z_intake_m']),
            dP_Pa=bar_to_Pa(float(pump_cfg['dP_bar'])),
        )

    if verbose:
        print(comm.describe())
        print(f"\nWell to reference: {well!r}")
        print(f"Thermal context:    {thermal!r}")
        print(f"GWR_std = {GWR:.2f}, NaCl_avg = "
              f"{config['fluid']['NaCl_ppm']:.0f} ppm, WHP = {WHP:.1f} bar\n")

    cfg_serialisable = copy.deepcopy(config)
    with open(os.path.join(out_dir, 'config.json'), 'w') as f:
        json.dump(cfg_serialisable, f, indent=2)

    op = solve_operating_point_multi(
        WHP, well, GWR, m_NaCl, thermal, comm,
        n_segments=n_seg, pump=pump)
    if verbose:
        print_operating_point_multi(op)

    export_summary_csv_multi(op, os.path.join(out_dir, 'summary.csv'))
    export_layer_breakdown_csv(comm,
        os.path.join(out_dir, 'layer_breakdown.csv'))

    plot_commingled_ipr(comm,
        save_path=os.path.join(out_dir, 'ipr_commingled.png'),
        title=f"Commingled IPR  -  {config.get('case_name','')}")
    plot_layer_split_vs_q(comm,
        save_path=os.path.join(out_dir, 'layer_split.png'),
        title=f"Layer-by-layer split  -  {config.get('case_name','')}")
    plot_commingled_nodal(comm, WHP, well, GWR, m_NaCl, thermal, op,
        save_path=os.path.join(out_dir, 'nodal_commingled.png'),
        title=f"Nodal plot  -  {config.get('case_name','')}",
        pump=pump)

    if op['converged']:
        from plotting import plot_pressure_profile, plot_temperature_profile
        plot_pressure_profile(op['profile'],
            title=f"P(z)  -  q={fmt_q(op['q_op_m3h'])}, "
                  f"WHP={WHP:.1f} bar",
            save_path=os.path.join(out_dir, 'profile_P.png'))
        plot_temperature_profile(op['profile'],
            title=f"T(z)",
            save_path=os.path.join(out_dir, 'profile_T.png'))

    files = sorted(os.listdir(out_dir))
    print(f"\n  Wrote {len(files)} files to {out_dir}/")
    for f in files:
        full = os.path.join(out_dir, f)
        print(f"    {f:<35s} ({os.path.getsize(full)//1024 + 1:4d} KB)")
    return op


if __name__ == "__main__":
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    DEFAULT_OUT = os.path.join(SCRIPT_DIR, "run_default")

    if len(sys.argv) > 1 and sys.argv[1] == 'calibrated':
        # Prefer the on-disk calibrated_config.json (fitted by
        # forchheimer_calibration.py and run_calibrated.py) over the
        # in-module CALIBRATED_CONFIG, which only holds initial seeds.
        cal_path = os.path.join(SCRIPT_DIR, 'calibrated_config.json')
        if os.path.exists(cal_path):
            with open(cal_path) as f:
                cfg = json.load(f)
        else:
            cfg = CALIBRATED_CONFIG
        out_dir = (sys.argv[2] if len(sys.argv) > 2
                   else os.path.join(SCRIPT_DIR, "run_calibrated"))
    elif len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            cfg = json.load(f)
        out_dir = sys.argv[2] if len(sys.argv) > 2 else DEFAULT_OUT
    else:
        cfg = DEFAULT_CONFIG
        out_dir = DEFAULT_OUT
    print(f"Output directory: {out_dir}")
    run_full_analysis(cfg, out_dir)
