"""
vlp.py
======

VLP (Vertical Lift Performance) module for vertical CO2-brine geothermal
production wells.

Implements:
  * Modified Hagedorn-Brown (1965) two-phase pressure-gradient correlation
    with the Griffith (1962) bubble-flow correction and the Brill-
    Mukherjee no-slip lower bound on liquid holdup.
  * Single-phase Darcy-Weisbach pressure gradient with Haaland (1983)
    friction factor for the under-saturated section of the wellbore.
  * Full Ramey (1962) wellbore heat-transfer model with the Hasan-Kabir
    (1991) transient time function f(t).
  * A top-down pressure-marching solver: given wellhead pressure and
    flow rate, returns the full P(z), T(z), H_L(z), regime(z) profile
    and the bottomhole flowing pressure (BHFP).

All inputs and outputs are SI; the module uses the PVT functions in
pvt.py for fluid properties.

References
----------
Hagedorn, A.R. and Brown, K.E. (1965). Experimental study of pressure
    gradients occurring during continuous two-phase flow in small-
    diameter vertical conduits. JPT 17(4), 475-484.
Griffith, P. and Wallis, G.B. (1961). Two-phase slug flow.
    J. Heat Transfer 83(3), 307-318.
Brill, J.P. and Mukherjee, H. (1999). Multiphase Flow in Wells.
    SPE Monograph Vol. 17, Henry L. Doherty Series.
Haaland, S.E. (1983). Simple and explicit formulas for the friction
    factor in turbulent pipe flow. J. Fluids Eng. 105, 89-90.
Ramey, H.J. (1962). Wellbore heat transmission. JPT 14(4), 427-435.
Hasan, A.R. and Kabir, C.S. (1991). Heat transfer during two-phase
    flow in wellbores: Part I - formation temperature. SPE 22866.
"""

from __future__ import annotations
import numpy as np

from pvt import (
    R_GAS, M_CO2, M_NACL, T_STD, P_STD,
    C_to_K, K_to_C, bar_to_Pa, Pa_to_bar, ppm_to_molality,
    molality_to_mass_fraction,
    density_water_BW, density_brine_BW, viscosity_brine,
    density_CO2, viscosity_CO2,
    flash_CO2_brine, GWR_std_to_mass_fraction,
)

GRAVITY = 9.80665      # m/s^2
P_ATM   = 101325.0     # Pa


# =====================================================================
# 1. Configuration containers (lightweight - just hold inputs)
# =====================================================================
class WellGeometry:
    """Vertical well geometry.

    Parameters
    ----------
    depth_TVD     : true vertical depth of perforations (m)
    tubing_ID     : flow conduit inner diameter (m).  In a two-conduit
                    well (with `pump_depth_m` set), this is the ID of
                    the tubing ABOVE the pump.
    tubing_OD     : tubing outer diameter (m); defaults to 1.14 * ID
                    (typical L-80 schedule).  Used as the reference
                    surface for the Ramey heat-transfer coefficient U.
    wellbore_dia  : drilled hole diameter (m); defaults to 2 * tubing_OD
                    (used in the Hasan-Kabir transient time function).
    roughness     : absolute pipe roughness (m); default 4.6e-5 (new
                    commercial steel).
    pump_depth_m  : optional ESP intake depth (m TVD).  When set, the
                    flow conduit switches at this depth from the
                    `tubing_*` parameters (above) to the `casing_*`
                    parameters (below).
    casing_ID     : inner diameter of the flow conduit BELOW the pump
                    (m).  Default = tubing_ID (single-conduit fallback).
    casing_OD     : outer diameter of the casing string (m).  Default
                    = 1.14 * casing_ID.  Used for Ramey U evaluation
                    in the section below the pump.
    casing_roughness : pipe roughness below pump (m).  Default = roughness.
    segments     : optional list of (z_top, z_bot, ID, OD, eps) tuples
                   describing N arbitrary flow conduits stacked along
                   the TVD axis.  When provided, this OVERRIDES the
                   two-conduit (tubing/casing) logic and is the
                   preferred way to model a multi-conduit completion
                   (e.g. tubing + casing + slotted liner).  Segments
                   must cover [0, depth_TVD] contiguously, ordered
                   from shallow to deep, with no gaps or overlaps.
                   The tubing_ID argument is then used only as a
                   fallback for any z outside the covered range.
                   Example for VinGT-1 eruptive (5" DP + 9-5/8" csg
                   + 7" slotted liner):
                       segments=[
                           (0,    1438,   0.1086, 0.1270, 4.6e-5),
                           (1438, 1499,   0.2205, 0.2445, 4.6e-5),
                           (1499, z_ref,  0.1571, 0.1778, 4.6e-5),
                       ]
                   pump_depth_m can still be set independently to
                   inject the ESP boost at the right depth.
    """

    def __init__(self, depth_TVD, tubing_ID=None, tubing_OD=None,
                 wellbore_dia=None, roughness=4.6e-5,
                 pump_depth_m=None,
                 casing_ID=None, casing_OD=None,
                 casing_roughness=None,
                 segments=None):
        self.L         = float(depth_TVD)
        # ---- Optional N-segment description (preferred for multi-conduit) ----
        # Each entry is (z_top, z_bot, ID, OD, eps), with ID/OD/eps in metres.
        self.segments = None
        if segments is not None:
            segs = []
            for entry in segments:
                if len(entry) == 4:
                    z_top, z_bot, ID, OD = entry
                    eps = roughness
                elif len(entry) == 5:
                    z_top, z_bot, ID, OD, eps = entry
                else:
                    raise ValueError(
                        f"WellGeometry.segments: each entry must be "
                        f"(z_top, z_bot, ID, OD[, eps]); got {entry!r}")
                segs.append((float(z_top), float(z_bot),
                              float(ID), float(OD), float(eps)))
            # Sort shallow -> deep and sanity-check coverage
            segs.sort(key=lambda s: s[0])
            for i in range(1, len(segs)):
                if abs(segs[i][0] - segs[i-1][1]) > 1e-6:
                    raise ValueError(
                        f"WellGeometry.segments: gap or overlap between "
                        f"segment {i-1} ({segs[i-1][:2]}) and "
                        f"segment {i} ({segs[i][:2]})")
            self.segments = segs
            # Pick a fallback (used by geom_at when z is outside coverage,
            # e.g. just under the wellhead before the deepest node is set):
            # use the SHALLOWEST segment's geometry.
            top_seg = segs[0]
            if tubing_ID is None:
                tubing_ID = top_seg[2]
            if tubing_OD is None:
                tubing_OD = top_seg[3]
            if roughness is None:
                roughness = top_seg[4]
        elif tubing_ID is None:
            raise ValueError(
                "WellGeometry: either tubing_ID or segments must be provided.")

        # ---- above pump (or single conduit) ----
        self.d_ti      = float(tubing_ID)
        self.d_to      = float(tubing_OD) if tubing_OD is not None \
                         else 1.14 * self.d_ti
        self.d_w       = float(wellbore_dia) if wellbore_dia is not None \
                         else 2.0 * self.d_to
        self.eps       = float(roughness)
        self.r_to      = 0.5 * self.d_to
        self.r_w       = 0.5 * self.d_w
        self.A_p       = 0.25 * np.pi * self.d_ti ** 2

        # ---- below pump (two-conduit) ----
        self.pump_depth_m = (float(pump_depth_m)
                              if pump_depth_m is not None else None)
        if casing_ID is not None:
            self.d_ti_below = float(casing_ID)
        else:
            self.d_ti_below = self.d_ti
        if casing_OD is not None:
            self.d_to_below = float(casing_OD)
        else:
            self.d_to_below = 1.14 * self.d_ti_below
        self.d_w_below  = max(self.d_w, 2.0 * self.d_to_below)
        self.eps_below  = (float(casing_roughness)
                           if casing_roughness is not None else self.eps)
        self.r_to_below = 0.5 * self.d_to_below
        self.r_w_below  = 0.5 * self.d_w_below
        self.A_p_below  = 0.25 * np.pi * self.d_ti_below ** 2

    def is_two_conduit(self):
        return (self.pump_depth_m is not None
                 and abs(self.d_ti_below - self.d_ti) > 1e-9)

    def geom_at(self, z_m):
        """Return (A_p, d_ti, eps) for the flow conduit at depth z.

        Priority order:
            1. If `segments` was given, look up the segment containing z.
            2. Else if `pump_depth_m` is set, switch tubing/casing at it.
            3. Else: uniform single-conduit (tubing_* everywhere).
        """
        # --- (1) N-segment lookup ---
        if self.segments is not None:
            for (z_top, z_bot, ID, OD, eps) in self.segments:
                # Use z >= z_top and z < z_bot, with the deepest segment
                # also covering its bottom node exactly.
                if z_top - 1e-9 <= z_m < z_bot + 1e-9:
                    A_p = 0.25 * np.pi * ID * ID
                    return A_p, ID, eps
            # z outside coverage -> use shallowest segment as fallback
            z_top, z_bot, ID, OD, eps = self.segments[0]
            A_p = 0.25 * np.pi * ID * ID
            return A_p, ID, eps
        # --- (2) two-conduit (pump) fallback ---
        if self.pump_depth_m is None or z_m < self.pump_depth_m:
            return self.A_p, self.d_ti, self.eps
        return self.A_p_below, self.d_ti_below, self.eps_below

    def __repr__(self):
        if self.segments is not None:
            parts = []
            for (z_top, z_bot, ID, OD, eps) in self.segments:
                parts.append(f"{z_top:.0f}-{z_bot:.0f}m: ID={ID*1000:.1f}mm")
            s = (f"WellGeometry(L={self.L:.1f} m, "
                 f"segments=[{', '.join(parts)}]")
            if self.pump_depth_m is not None:
                s += f", pump_depth={self.pump_depth_m:.0f}m"
            s += ")"
            return s
        s = (f"WellGeometry(L={self.L:.1f} m, d_ti={self.d_ti*1000:.1f} mm, "
             f"d_to={self.d_to*1000:.1f} mm, d_w={self.d_w*1000:.1f} mm, "
             f"eps={self.eps*1e6:.0f} um")
        if self.is_two_conduit():
            s += (f"; below {self.pump_depth_m:.0f} m: "
                  f"d_ti={self.d_ti_below*1000:.1f} mm)")
        else:
            s += ")"
        return s


class ThermalContext:
    """Thermal / formation properties for the Ramey heat-transfer model.

    Parameters
    ----------
    T_surface    : surface ambient temperature (K)
    geo_gradient : geothermal gradient (K/m, positive downward)
    T_BH         : bottom-hole flowing-fluid temperature at the
                   perforations (K).  Default = T_surface + g * L (i.e.
                   undisturbed formation T at the perfs - assumes no
                   significant Joule-Thomson cooling on inflow).
    k_formation  : formation thermal conductivity (W/m/K), typ. 1.5-3.0
                   for sediments, 2.5-4.0 for granitic basement.
    alpha_formation : formation thermal diffusivity (m^2/s), typ. ~1e-6.
    U_overall    : overall heat-transfer coefficient referenced to the
                   tubing OD (W/m^2/K).  Typical ranges:
                      - bare tubing in cemented casing : 15-30
                      - vacuum-insulated tubing       :  1-3
                      - fully insulated, gas-filled annulus : 5-10
    time_seconds : production time used in the Hasan-Kabir transient
                   function f(t).  Use a representative steady-state
                   duration (e.g. 30-90 days for design).
    """

    def __init__(self, T_surface, geo_gradient,
                 T_BH=None,
                 k_formation=2.5, alpha_formation=1.0e-6,
                 U_overall=20.0,
                 time_seconds=30.0 * 86400.0):
        self.T_surface = float(T_surface)
        self.geo_grad  = float(geo_gradient)
        self.T_BH      = float(T_BH) if T_BH is not None else None
        self.k_e       = float(k_formation)
        self.alpha_e   = float(alpha_formation)
        self.U         = float(U_overall)
        self.t         = float(time_seconds)


class FluidStream:
    """Composition and mass flow rate of the produced stream.

    The stream is described by:
       - brine volumetric rate at standard conditions, m^3/h
       - GWR at standard conditions (m^3 CO2 / m^3 brine)
       - NaCl molality of the brine

    Provides the total mass flow (W_total) and total CO2 mass fraction
    (w_CO2) used by the marching solver.
    """

    def __init__(self, q_brine_std_m3h, GWR_std=0.0, m_NaCl=0.0):
        self.q_brine_std = float(q_brine_std_m3h) / 3600.0   # m^3/s
        self.GWR_std     = float(GWR_std)
        self.m_NaCl      = float(m_NaCl)

        rho_brine_std = density_brine_BW(T_STD, P_STD, m_NaCl)
        rho_CO2_std   = density_CO2(T_STD, P_STD, phase='gas')

        self.W_brine = self.q_brine_std * rho_brine_std         # kg/s
        self.W_CO2   = self.GWR_std * self.q_brine_std * rho_CO2_std
        self.W_total = self.W_brine + self.W_CO2
        self.w_CO2   = self.W_CO2 / self.W_total if self.W_total > 0 else 0.0

    def __repr__(self):
        q_m3h = self.q_brine_std * 3600.0
        q_ls  = q_m3h / 3.6
        q_m3d = q_m3h * 24.0
        return (f"FluidStream(q_brine={q_ls:.2f} l/s "
                f"({q_m3d:.0f} m^3/d), "
                f"GWR={self.GWR_std:.2f}, m_NaCl={self.m_NaCl:.3f} mol/kg, "
                f"W_total={self.W_total:.2f} kg/s, "
                f"w_CO2={self.w_CO2*100:.3f} wt%)")


# =====================================================================
# 2. Specific heat of brine (used by Ramey)
# =====================================================================
def cp_brine(T_K, m_NaCl):
    """Brine specific heat (J/kg/K).  Engineering correlation.

    Pure-water cp from a polynomial fit to NIST data (within 1 % over
    273-573 K).  NaCl correction following Phillips et al. (1981) -
    monotonic decrease with salinity.
    """
    T_C = T_K - 273.15
    # Pure water cp - shallow minimum near 35 degC, rises above 100 degC
    cp_w = 4180.0 + 0.013 * (T_C - 35.0) ** 2
    # Salinity correction: each kg NaCl contributes ~850 J/K instead of
    # ~4180; linearise as a small reduction in mixture cp.
    S = molality_to_mass_fraction(m_NaCl)
    cp = cp_w * (1.0 - 0.45 * S)
    return cp


# =====================================================================
# 3. Ramey heat transfer
# =====================================================================
def hasan_kabir_time_function(t_seconds, alpha_e, r_w):
    """Hasan-Kabir (1991) transient time function f(t_D).

    t_D = alpha_e * t / r_w^2  (dimensionless time)

    Smoothly covers both early-time (t_D < 1.5) and late-time regimes.
    """
    tD = alpha_e * t_seconds / r_w ** 2
    if tD <= 0:
        return 0.0
    if tD <= 1.5:
        return 1.1281 * np.sqrt(tD) * (1.0 - 0.3 * np.sqrt(tD))
    return (0.4063 + 0.5 * np.log(tD)) * (1.0 + 0.6 / tD)


def ramey_relaxation_distance(W_dot, cp, U, r_to, k_e, f_t):
    """Ramey (1962) relaxation distance A (m).

    A = (W * cp) / (2 pi) * [ 1/(U r_to) + f(t)/k_e ]

    The bracketed term is the total thermal resistance per unit length
    from flowing fluid to undisturbed formation, in m K / W.
    """
    R_per_L = 1.0 / (U * r_to) + f_t / k_e          # m K / W per unit length factor
    return W_dot * cp * R_per_L / (2.0 * np.pi)


def ramey_temperature_profile(z_array, L, T_surf, geo_grad, T_BH, A):
    """Steady-state Ramey production-well fluid temperature.

    Parameters
    ----------
    z_array  : array of depths from surface, m (z=0 surface, z=L BH)
    L        : total TVD, m
    T_surf   : surface ambient temperature, K
    geo_grad : geothermal gradient, K/m (positive downward)
    T_BH     : bottom-hole flowing-fluid temperature, K
    A        : Ramey relaxation distance, m

    Returns
    -------
    T_f : array of fluid temperatures at z_array nodes (K)
    """
    z = np.asarray(z_array, dtype=float)
    T_e_BH = T_surf + geo_grad * L
    s = L - z                                            # distance from BH (up)
    T_e = T_surf + geo_grad * z                          # local formation T
    T_f = T_e + geo_grad * A * (1.0 - np.exp(-s / A))
    # Optional anomaly: fluid leaves reservoir at T_BH != T_e_BH
    if abs(T_BH - T_e_BH) > 0.01:
        T_f = T_f + (T_BH - T_e_BH) * np.exp(-s / A)
    return T_f


# =====================================================================
# 4. Friction factor (Haaland 1983)
# =====================================================================
def darcy_friction_factor(Re, eps_over_d):
    """Darcy friction factor.  Laminar for Re<2300, Haaland otherwise."""
    if Re < 1.0:
        return 64.0
    if Re < 2300.0:
        return 64.0 / Re
    inv_sqrt_f = -1.8 * np.log10((eps_over_d / 3.7) ** 1.11 + 6.9 / Re)
    return 1.0 / inv_sqrt_f ** 2


# =====================================================================
# 5. Modified Hagedorn-Brown two-phase holdup
# =====================================================================
def _CN_L(N_L):
    """Hagedorn-Brown viscosity-number correction CN_L.
    Brill-Mukherjee polynomial fit to the original H&B chart."""
    X = np.log10(max(N_L, 1.0e-12)) + 3.0
    Y = (-2.69851
         + 0.15841 * X
         - 0.55100 * X ** 2
         + 0.54785 * X ** 3
         - 0.12195 * X ** 4)
    return 10.0 ** Y


def _HL_over_psi(phi):
    """Hagedorn-Brown holdup correlating function H_L / psi.
    Brill-Mukherjee polynomial fit to the H&B holdup chart."""
    X = np.log10(max(phi, 1.0e-12)) + 6.0
    Y = (-0.10307
         + 0.61777 * X
         - 0.63295 * X ** 2
         + 0.29598 * X ** 3
         - 0.04010 * X ** 4)
    # The fit overshoots beyond its calibration range - clamp.
    return float(np.clip(Y, 0.0, 1.0))


def _psi_correction(chi):
    """Hagedorn-Brown psi (secondary) correction factor."""
    if chi <= 0.012:
        return 1.0
    psi = (0.91163
           - 4.82176 * chi
           + 1232.25 * chi ** 2
           - 22253.6 * chi ** 3
           + 116174.3 * chi ** 4)
    return max(1.0, psi)


def hagedorn_brown_holdup(v_sl, v_sg, rho_l, rho_g, mu_l, mu_g, sigma,
                          P_Pa, d_m, g=GRAVITY):
    """Modified Hagedorn-Brown liquid holdup with Griffith bubble check.

    Returns
    -------
    H_l    : liquid volumetric holdup (-)
    regime : 'bubble' or 'slug' (informational)
    """
    v_m = v_sl + v_sg
    if v_m <= 0.0:
        return 1.0, 'no_flow'
    lambda_l = v_sl / v_m

    # Griffith bubble-flow boundary (transition to slug)
    L_B = max(1.071 - 0.2218 * v_m ** 2 / d_m, 0.13)
    if v_sg / v_m < L_B:
        # Bubble flow (Griffith-Wallis)
        v_s = 0.244          # m/s, bubble slip velocity
        rad = (1.0 + v_m / v_s) ** 2 - 4.0 * v_sg / v_s
        if rad < 0.0:
            rad = 0.0
        H_l = 1.0 - 0.5 * (1.0 + v_m / v_s - np.sqrt(rad))
        H_l = max(H_l, lambda_l)
        H_l = min(H_l, 1.0)
        return H_l, 'bubble'

    # Slug / intermittent: standard Hagedorn-Brown
    # Dimensionless H&B numbers (consistent SI form)
    sg_term = (rho_l / (g * sigma)) ** 0.25
    N_lv = v_sl * sg_term
    N_gv = v_sg * sg_term
    N_d  = d_m * np.sqrt(rho_l * g / sigma)
    N_l  = mu_l * (g / (rho_l * sigma ** 3)) ** 0.25

    CN_l = _CN_L(N_l)

    phi = ((N_lv / N_gv ** 0.575)
           * (P_Pa / P_ATM) ** 0.1
           * CN_l / N_d)

    HL_psi = _HL_over_psi(phi)

    chi = N_gv * N_l ** 0.380 / N_d ** 2.14
    psi = _psi_correction(chi)

    H_l = HL_psi * psi
    # Brill-Mukherjee modification: never go below the no-slip holdup
    H_l = max(H_l, lambda_l)
    H_l = min(H_l, 1.0)
    return H_l, 'slug'


# =====================================================================
# 6. Total pressure gradient (single-phase or H&B)
# =====================================================================
def pressure_gradient(P_Pa, T_K, m_NaCl, w_CO2_total,
                      W_dot, A_p, d_m, eps,
                      g=GRAVITY):
    """Local steady-state pressure gradient (Pa/m, positive going DOWN).

    Branches automatically between single-phase liquid (Darcy-Weisbach)
    and two-phase (Modified Hagedorn-Brown) based on the local CO2
    flash result.

    Returns
    -------
    dPdL_total : Pa/m (positive means pressure rises with depth)
    flash      : the dict returned by `flash_CO2_brine`
    info       : dict with hydrostatic and friction breakdown plus any
                 two-phase quantities
    """
    flash = flash_CO2_brine(T_K, P_Pa, m_NaCl, w_CO2_total)

    if not flash['is_two_phase']:
        # ---------- Single-phase liquid ----------
        rho = flash['rho_liq']
        mu  = flash['mu_liq']
        v   = W_dot / (rho * A_p)
        Re  = rho * v * d_m / mu
        f   = darcy_friction_factor(Re, eps / d_m)
        dp_g = rho * g
        dp_f = f * rho * v ** 2 / (2.0 * d_m)
        return (dp_g + dp_f, flash,
                dict(regime='single_phase_L', H_l=1.0,
                     rho_m=rho, rho_ns=rho,
                     v_sl=v, v_sg=0.0, v_m=v,
                     mu_m=mu, Re=Re, f=f,
                     dp_grav=dp_g, dp_fric=dp_f))

    # ---------- Two-phase (Modified Hagedorn-Brown) ----------
    rho_l = flash['rho_liq']
    rho_g = flash['rho_gas']
    mu_l  = flash['mu_liq']
    mu_g  = flash['mu_gas']
    sigma = flash['sigma_gl']

    W_l = W_dot * flash['x_liq']
    W_g = W_dot * flash['x_gas']
    q_l = W_l / rho_l
    q_g = W_g / rho_g
    v_sl = q_l / A_p
    v_sg = q_g / A_p
    v_m  = v_sl + v_sg

    H_l, regime = hagedorn_brown_holdup(v_sl, v_sg, rho_l, rho_g,
                                        mu_l, mu_g, sigma, P_Pa, d_m, g)
    H_g = 1.0 - H_l
    lambda_l = v_sl / v_m
    lambda_g = 1.0 - lambda_l

    rho_m  = rho_l * H_l + rho_g * H_g                   # slip mixture density
    rho_ns = rho_l * lambda_l + rho_g * lambda_g          # no-slip
    mu_ns  = mu_l * lambda_l + mu_g * lambda_g            # no-slip viscosity

    # H&B Reynolds and friction:
    Re = rho_ns * v_m * d_m / mu_ns
    f  = darcy_friction_factor(Re, eps / d_m)

    dp_g = rho_m * g
    # H&B friction form: f * G^2 / (2 d rho_m), with G = rho_ns * v_m
    dp_f = f * rho_ns ** 2 * v_m ** 2 / (2.0 * d_m * rho_m)

    return (dp_g + dp_f, flash,
            dict(regime=regime, H_l=H_l,
                 rho_m=rho_m, rho_ns=rho_ns,
                 v_sl=v_sl, v_sg=v_sg, v_m=v_m,
                 mu_m=mu_ns, Re=Re, f=f,
                 dp_grav=dp_g, dp_fric=dp_f))


# =====================================================================
# 7. Marching solver: WHP -> BHFP (top-down)
# =====================================================================
def march_VLP(WHP_Pa, well, fluid, thermal,
              n_segments=50, tol=50.0, max_iter=20,
              pump=None):
    """March pressure from wellhead down to perforations.

    Parameters
    ----------
    WHP_Pa     : wellhead pressure (Pa)
    well       : WellGeometry
    fluid      : FluidStream
    thermal    : ThermalContext
    n_segments : number of axial segments along TVD
    tol        : per-segment pressure-convergence tolerance (Pa)
    max_iter   : max iterations per segment
    pump       : optional dict describing an ESP boost, with keys
                   z_intake_m : float   - intake depth (m TVD).  Below
                                          this depth, the marched
                                          pressure DROPS by dP_Pa to
                                          model the suction side of
                                          the pump.
                   dP_Pa      : float   - pump head expressed as a
                                          pressure rise (Pa) between
                                          intake and discharge.  At
                                          the segment containing
                                          z_intake_m, the marched
                                          pressure is decreased by
                                          dP_Pa (since we march top-
                                          down: above pump we have the
                                          pump-boosted column with
                                          relatively *low* P at WH;
                                          below the pump we have the
                                          un-boosted reservoir column,
                                          which is HIGHER than what
                                          continued marching from
                                          above would give -- hence we
                                          ADD dP_Pa back to keep
                                          Pwf physically consistent).
                 If None, no pump is modelled.

    Returns
    -------
    dict with arrays:  z, P, T, T_e, H_l, rho_m, v_m, regime,
       dpdL, dpdL_grav, dpdL_fric  (length n+1)
    plus scalars  BHFP, A_relax, T_e_BH, x_two_phase_top,
       and -- if pump is not None --
       z_pump_m, dP_pump_Pa, P_intake_Pa, P_discharge_Pa.
    """
    L = well.L
    z_nodes = np.linspace(0.0, L, n_segments + 1)
    dz = L / n_segments

    # ---- Ramey: precompute T(z) profile ----
    T_e_BH = thermal.T_surface + thermal.geo_grad * L
    T_BH = thermal.T_BH if thermal.T_BH is not None else T_e_BH
    T_avg_est = 0.5 * (T_BH + thermal.T_surface)
    cp_avg = cp_brine(T_avg_est, fluid.m_NaCl)

    f_t = hasan_kabir_time_function(thermal.t, thermal.alpha_e, well.r_w)
    A_relax = ramey_relaxation_distance(
        fluid.W_total, cp_avg, thermal.U, well.r_to, thermal.k_e, f_t)

    T_array  = ramey_temperature_profile(z_nodes, L, thermal.T_surface,
                                         thermal.geo_grad, T_BH, A_relax)
    T_e_arr  = thermal.T_surface + thermal.geo_grad * z_nodes

    # ---- Pressure marching (top down, segment by segment) ----
    P_array       = np.zeros(n_segments + 1)
    H_l_array     = np.zeros(n_segments + 1)
    rho_m_array   = np.zeros(n_segments + 1)
    v_m_array     = np.zeros(n_segments + 1)
    dpdL_array    = np.zeros(n_segments + 1)
    dpdL_grav_arr = np.zeros(n_segments + 1)
    dpdL_fric_arr = np.zeros(n_segments + 1)
    regime_array  = [''] * (n_segments + 1)

    # Top boundary
    P_array[0] = WHP_Pa
    Ap_top, di_top, eps_top = well.geom_at(z_nodes[0])
    grad_top, _, info_top = pressure_gradient(
        P_array[0], T_array[0], fluid.m_NaCl, fluid.w_CO2,
        fluid.W_total, Ap_top, di_top, eps_top)
    H_l_array[0]   = info_top['H_l']
    rho_m_array[0] = info_top['rho_m']
    v_m_array[0]   = info_top['v_m']
    dpdL_array[0]  = grad_top
    dpdL_grav_arr[0] = info_top['dp_grav']
    dpdL_fric_arr[0] = info_top['dp_fric']
    regime_array[0]  = info_top['regime']

    bubble_point_depth = None    # deepest depth at which 2-phase exists

    # Pump-boost bookkeeping
    P_intake_Pa = None
    P_discharge_Pa = None
    z_pump = float(pump['z_intake_m']) if pump is not None else None
    dP_pump = float(pump['dP_Pa']) if pump is not None else 0.0
    pump_applied = False

    for i in range(n_segments):
        P_in = P_array[i]
        T_in = T_array[i]
        T_out = T_array[i + 1]
        T_avg = 0.5 * (T_in + T_out)
        z_avg = 0.5 * (z_nodes[i] + z_nodes[i + 1])
        Ap_loc, di_loc, eps_loc = well.geom_at(z_avg)

        # Initial guess for P_out: liquid hydrostatic
        P_out = P_in + 1000.0 * GRAVITY * dz
        for it in range(max_iter):
            P_avg = 0.5 * (P_in + P_out)
            # Sanity: if P has run away to non-physical values, bail.
            # This protects PR-EOS / flash from NaN inputs and lets the
            # outer nodal Brent solver catch the failure.
            if (not np.isfinite(P_avg)) or P_avg <= 0.0 or P_avg > 5.0e8:
                raise ValueError(
                    f"march_VLP: P diverged at z={z_nodes[i+1]:.0f} m "
                    f"(P_avg={P_avg:.3e} Pa); rate likely too high "
                    f"for this conduit.")
            dpdL, _flash, _info = pressure_gradient(
                P_avg, T_avg, fluid.m_NaCl, fluid.w_CO2,
                fluid.W_total, Ap_loc, di_loc, eps_loc)
            P_out_new = P_in + dpdL * dz
            if abs(P_out_new - P_out) < tol:
                P_out = P_out_new
                break
            # mild damping for stability when crossing bubble point
            P_out = 0.6 * P_out_new + 0.4 * P_out

        # ---- Apply ESP pump boost at the segment containing z_pump ----
        # We march TOP-DOWN.  Just above the pump (in the tubing) the
        # fluid is at the DISCHARGE pressure -- relatively low, because
        # the discharge side only has to support the tubing column up
        # to the wellhead.  Just below the pump (in the casing/annulus,
        # being fed by the reservoir below) the fluid is at the INTAKE
        # pressure, which is dP_pump LOWER than the discharge.  So
        # marching top-down across the pump we SUBTRACT dP_pump.
        #
        # The reservoir then has to lift the fluid from Pwf_ref up
        # through the casing to the intake; the marched pressure
        # naturally builds back up with hydrostatic + friction in
        # 9-5/8" casing below z_pump.
        if (pump is not None and (not pump_applied)
                and z_nodes[i + 1] >= z_pump):
            P_discharge_Pa = float(P_out)
            P_out = P_out - dP_pump
            P_intake_Pa = float(P_out)
            pump_applied = True

        P_array[i + 1] = P_out

        # Recompute properties at the new (output) node for storage
        Ap_out, di_out, eps_out = well.geom_at(z_nodes[i + 1])
        grad_out, _flash_out, info_out = pressure_gradient(
            P_array[i + 1], T_array[i + 1], fluid.m_NaCl, fluid.w_CO2,
            fluid.W_total, Ap_out, di_out, eps_out)
        H_l_array[i + 1]   = info_out['H_l']
        rho_m_array[i + 1] = info_out['rho_m']
        v_m_array[i + 1]   = info_out['v_m']
        dpdL_array[i + 1]  = grad_out
        dpdL_grav_arr[i + 1] = info_out['dp_grav']
        dpdL_fric_arr[i + 1] = info_out['dp_fric']
        regime_array[i + 1]  = info_out['regime']

        # Track the deepest depth at which two-phase still exists.
        # Marching top-down, the last 2-phase node we see is just above
        # the bubble-point crossing.
        if info_out['regime'] != 'single_phase_L':
            bubble_point_depth = z_nodes[i + 1]

    # Integrated pressure-drop breakdown (trapezoidal over the profile).
    # np.trapezoid is the NumPy 2.x name; fall back to np.trapz on older.
    _trapz = getattr(np, 'trapezoid', getattr(np, 'trapz', None))
    dP_grav_total = float(_trapz(dpdL_grav_arr, z_nodes))
    dP_fric_total = float(_trapz(dpdL_fric_arr, z_nodes))
    dP_total      = float(P_array[-1] - P_array[0])

    out = dict(
        z=z_nodes, P=P_array, T=T_array, T_e=T_e_arr,
        H_l=H_l_array, rho_m=rho_m_array, v_m=v_m_array,
        regime=regime_array,
        dpdL=dpdL_array,
        dpdL_grav=dpdL_grav_arr, dpdL_fric=dpdL_fric_arr,
        WHP=float(WHP_Pa),
        BHFP=float(P_array[-1]),
        Pwf=float(P_array[-1]),                      # alias used in IPR
        dP_total=dP_total,
        dP_grav=dP_grav_total,
        dP_fric=dP_fric_total,
        A_relax=float(A_relax),
        T_e_BH=float(T_e_BH),
        T_BH_fluid=float(T_array[-1]),
        T_WH_fluid=float(T_array[0]),
        bubble_point_depth=bubble_point_depth,
    )
    if pump is not None:
        out['z_pump_m']        = z_pump
        out['dP_pump_Pa']      = dP_pump
        out['P_discharge_Pa']  = P_discharge_Pa
        out['P_intake_Pa']     = P_intake_Pa
    return out


def Pwf_from_WHP(WHP_Pa, well, fluid, thermal, n_segments=50, pump=None):
    """Convenience wrapper: just return BHFP for a given WHP."""
    return march_VLP(WHP_Pa, well, fluid, thermal,
                     n_segments=n_segments, pump=pump)['BHFP']


# =====================================================================
# 8. Pretty-printer for a marching result
# =====================================================================
def print_profile(result, every=5):
    """Print every Nth row of a march_VLP() result dict."""
    z = result['z']
    P = result['P']
    T = result['T']
    Te = result['T_e']
    Hl = result['H_l']
    rm = result['rho_m']
    vm = result['v_m']
    rg = result['regime']
    dpg = result['dpdL_grav']
    dpf = result['dpdL_fric']

    print(f"  {'z (m)':>7} {'P (bar)':>8} {'T (degC)':>9} "
          f"{'Te (degC)':>9} {'H_l':>5} {'rho_m':>7} "
          f"{'v_m':>6} {'dPg':>6} {'dPf':>6}  {'regime':<14}")
    print(f"  {'-'*7} {'-'*8} {'-'*9} {'-'*9} {'-'*5} {'-'*7} "
          f"{'-'*6} {'-'*6} {'-'*6}  {'-'*14}")
    n = len(z)
    idx = list(range(0, n, every))
    if (n - 1) not in idx:
        idx.append(n - 1)
    for i in idx:
        print(f"  {z[i]:7.1f} {P[i]*1e-5:8.2f} {T[i]-273.15:9.1f} "
              f"{Te[i]-273.15:9.1f} {Hl[i]:5.3f} {rm[i]:7.1f} "
              f"{vm[i]:6.2f} {dpg[i]*1e-2:6.1f} {dpf[i]*1e-2:6.1f}  "
              f"{rg[i]:<14}")
    print(f"  (dPg, dPf in mbar/m)")


# =====================================================================
# 9. SELF-TEST / VALIDATION BLOCK
# =====================================================================
def _validate():
    print("=" * 78)
    print("VLP MODULE - SELF-CONSISTENCY TESTS")
    print("=" * 78)

    # --------------------------------------------------------------
    # TEST 1: Pure single-phase brine -> must equal rho*g*L + friction
    # --------------------------------------------------------------
    print("\n[Test 1] Pure brine, no CO2 - hydrostatic + friction check")
    well1 = WellGeometry(depth_TVD=2000.0, tubing_ID=0.1)
    fluid1 = FluidStream(q_brine_std_m3h=200.0, GWR_std=0.0,
                         m_NaCl=ppm_to_molality(50_000))
    therm1 = ThermalContext(T_surface=C_to_K(20.0),
                            geo_gradient=0.090,         # 180 K over 2000 m
                            U_overall=20.0)
    res1 = march_VLP(bar_to_Pa(5.0), well1, fluid1, therm1, n_segments=80)
    print_profile(res1, every=10)
    print(f"  BHFP = {res1['BHFP']*1e-5:.2f} bar  "
          f"(WHP = 5.00 bar; depth = {well1.L:.0f} m)")
    print(f"  T at surface = {res1['T_WH_fluid']-273.15:.1f} degC  "
          f"(formation T at surface = {therm1.T_surface-273.15:.1f} degC)")
    print(f"  T at BH      = {res1['T_BH_fluid']-273.15:.1f} degC  "
          f"(formation T at BH = {res1['T_e_BH']-273.15:.1f} degC)")
    print(f"  Ramey relaxation distance A = {res1['A_relax']:.0f} m")
    print(f"  --- dP breakdown ---")
    print(f"     dP_total      = {res1['dP_total']*1e-5:7.2f} bar")
    print(f"     dP_grav (rho g L) = {res1['dP_grav']*1e-5:7.2f} bar  "
          f"({res1['dP_grav']/res1['dP_total']*100:.1f} %)")
    print(f"     dP_fric           = {res1['dP_fric']*1e-5:7.2f} bar  "
          f"({res1['dP_fric']/res1['dP_total']*100:.1f} %)")

    # --------------------------------------------------------------
    # TEST 2: Low-GWR geothermal case (the actual use case)
    # --------------------------------------------------------------
    print("\n[Test 2] Low-GWR geothermal: GWR=10 Sm^3/Sm^3")
    well2 = WellGeometry(depth_TVD=2500.0, tubing_ID=0.114)
    fluid2 = FluidStream(q_brine_std_m3h=250.0, GWR_std=10.0,
                         m_NaCl=ppm_to_molality(40_000))
    therm2 = ThermalContext(T_surface=C_to_K(15.0),
                            geo_gradient=0.080,         # 200 K over 2500 m
                            U_overall=20.0,
                            time_seconds=60.0 * 86400.0)
    res2 = march_VLP(bar_to_Pa(8.0), well2, fluid2, therm2, n_segments=100)
    print(fluid2)
    print_profile(res2, every=10)
    print(f"  BHFP = {res2['BHFP']*1e-5:.2f} bar")
    if res2['bubble_point_depth'] is not None:
        print(f"  Bubble-point crossing depth ~ "
              f"{res2['bubble_point_depth']:.0f} m TVD "
              f"(deepest 2-phase node)")
    print(f"  --- dP breakdown ---")
    print(f"     dP_total = {res2['dP_total']*1e-5:7.2f} bar  "
          f"(grav {res2['dP_grav']/res2['dP_total']*100:.1f} %, "
          f"fric {res2['dP_fric']/res2['dP_total']*100:.1f} %)")

    # --------------------------------------------------------------
    # TEST 3: GWR sensitivity - at high rates friction dominates,
    # so increasing GWR raises BHFP despite reduced hydrostatic head.
    # At low rates the opposite can happen (hydrostatic dominates).
    # --------------------------------------------------------------
    print("\n[Test 3] GWR sensitivity at fixed rate (250 m^3/h, WHP=8 bar)")
    fluid3 = FluidStream(q_brine_std_m3h=250.0, GWR_std=50.0,
                         m_NaCl=ppm_to_molality(40_000))
    res3 = march_VLP(bar_to_Pa(8.0), well2, fluid3, therm2, n_segments=100)
    print(fluid3)
    print(f"  BHFP(GWR=10) = {res2['BHFP']*1e-5:7.2f} bar   "
          f"grav {res2['dP_grav']*1e-5:6.1f} bar   "
          f"fric {res2['dP_fric']*1e-5:6.1f} bar")
    print(f"  BHFP(GWR=50) = {res3['BHFP']*1e-5:7.2f} bar   "
          f"grav {res3['dP_grav']*1e-5:6.1f} bar   "
          f"fric {res3['dP_fric']*1e-5:6.1f} bar")
    if res3['BHFP'] > res2['BHFP']:
        print("  => at this high flow rate in a 4.5in tubing, friction")
        print("     dominates: the extra gas raises v_m and dP_fric more")
        print("     than it cuts dP_grav, so BHFP rises with GWR.")
    else:
        print("  => at this rate, hydrostatic relief from the extra gas")
        print("     outweighs the added friction: BHFP drops with GWR.")

    # --------------------------------------------------------------
    # TEST 4: Rate sensitivity (just to verify monotonicity)
    # --------------------------------------------------------------
    print("\n[Test 4] Rate sensitivity (fixed WHP=8 bar, GWR=10)")
    print(f"  {'q_brine (m^3/h)':>16}  {'BHFP (bar)':>11}  "
          f"{'T at WH (degC)':>15}")
    for q in [50, 100, 200, 400, 800]:
        flu = FluidStream(q_brine_std_m3h=q, GWR_std=10.0,
                          m_NaCl=ppm_to_molality(40_000))
        r = march_VLP(bar_to_Pa(8.0), well2, flu, therm2, n_segments=60)
        print(f"  {q:16.0f}  {r['BHFP']*1e-5:11.2f}  "
              f"{r['T'][0]-273.15:15.1f}")
    print("  As rate goes up:")
    print("    - Friction loss increases (BHFP rises)")
    print("    - Less heat lost per kg fluid (WHP T closer to BH T)")

    # --------------------------------------------------------------
    # TEST 5: Heat-transfer coefficient sensitivity
    # --------------------------------------------------------------
    print("\n[Test 5] U sensitivity (fixed q=200 m^3/h, GWR=10)")
    print(f"  {'U (W/m2/K)':>11}  {'A_relax (m)':>12}  "
          f"{'T at WH (degC)':>15}  {'BHFP (bar)':>11}")
    fluid_test = FluidStream(q_brine_std_m3h=200.0, GWR_std=10.0,
                             m_NaCl=ppm_to_molality(40_000))
    for U in [2.0, 10.0, 30.0, 100.0]:
        therm = ThermalContext(T_surface=C_to_K(15.0), geo_gradient=0.080,
                               U_overall=U, time_seconds=60.0 * 86400.0)
        r = march_VLP(bar_to_Pa(8.0), well2, fluid_test, therm,
                      n_segments=60)
        print(f"  {U:11.1f}  {r['A_relax']:12.0f}  "
              f"{r['T'][0]-273.15:15.1f}  {r['BHFP']*1e-5:11.2f}")
    print("  Lower U (better insulation) -> larger A -> hotter WH fluid")

    print("\n" + "=" * 78)
    print("MODEL NOTES:")
    print("  * Hagedorn-Brown is calibrated for vertical bubble/slug flow.")
    print("    Outside this regime (annular/mist) errors grow.  For low-")
    print("    GWR geothermal we are firmly in the bubble/slug range.")
    print("  * Joule-Thomson cooling on inflow/expansion is neglected.")
    print("    For low GWR the JT effect is small (a few K).")
    print("  * Ramey assumes constant cp and U.  For wide T excursions")
    print("    the relaxation distance can be re-evaluated iteratively.")
    print("  * Acceleration term in dP/dL is omitted (typically <1 % for")
    print("    the conditions of interest).")
    print("=" * 78)


if __name__ == "__main__":
    _validate()
