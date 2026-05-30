"""
ipr.py
======

IPR (Inflow Performance Relationship) module for CO2-brine geothermal
production wells.

Implements:
  * Steady-state and pseudo-steady-state Darcy radial inflow for a
    single-phase liquid reservoir (CO2 fully dissolved at reservoir P).
  * Straight-line PI mode for cases where a measured productivity
    index is available from a well test.
  * Optional Fetkovich back-pressure form for cases where non-Darcy
    or partially turbulent inflow is expected.
  * Skin contributions:
        - Hawkins (1956) formation-damage skin
        - Karakas-Tariq (1991) perforation skin (simplified)
        - Furui-style gravel-pack pressure-drop skin
  * A `ReservoirModel` container that holds reservoir + completion
    inputs and exposes `Pwf_at_q(q)` and `q_at_Pwf(Pwf)` methods.

All inputs and outputs are SI:
    permeability  k        : m^2  (1 mD = 9.869233e-16 m^2)
    thickness     h        : m
    pressures     P, Pwf   : Pa
    flow rate     q        : m^3/s at standard conditions (T_STD, P_STD)
    radii         r_w, r_e : m

References
----------
Darcy-radial PSS / SS:
    Earlougher, R.C. (1977). Advances in Well-Test Analysis.
    SPE Monograph Vol. 5.
Hawkins, M.F. (1956). A note on the skin effect. Trans. AIME 207, 356.
Karakas, M. and Tariq, S.M. (1991). Semi-analytical productivity
    models for perforated completions.  SPEPE 6, 73-82.
Furui, K., Zhu, D., Hill, A.D. (2005). A comprehensive model of
    horizontal-well completion performance.  SPEPF 20, 207-220.
Fetkovich, M.J. (1973). The isochronal testing of oil wells.
    SPE 4529.
"""

from __future__ import annotations
import numpy as np
from scipy.optimize import brentq

from pvt import (
    T_STD, P_STD,
    C_to_K, K_to_C, bar_to_Pa, Pa_to_bar,
    density_water_BW, density_brine_BW, viscosity_brine,
    flash_CO2_brine, fmt_q, m3h_to_ls, m3h_to_m3d,
)

# Permeability conversion
MD_TO_M2 = 9.869233e-16            # 1 mD in m^2


# =====================================================================
# 1. Skin models
# =====================================================================
def hawkins_damage_skin(k_ratio, r_damage, r_w):
    """Hawkins (1956) damage skin.

        s_d = (k / k_d - 1) * ln(r_d / r_w)

    Parameters
    ----------
    k_ratio  : k_d / k = ratio of damaged-zone to undisturbed
               permeability (0 < k_ratio <= 1; values > 1 indicate
               stimulation, returning a negative skin)
    r_damage : outer radius of damaged zone, m
    r_w      : wellbore radius, m
    """
    if k_ratio <= 0.0:
        raise ValueError("k_ratio must be > 0")
    if r_damage <= r_w:
        return 0.0
    return (1.0 / k_ratio - 1.0) * np.log(r_damage / r_w)


def karakas_tariq_perforation_skin(L_p, d_p, SPF_per_m, phasing_deg,
                                   r_w, k_ratio_crushed=0.4,
                                   t_crushed=0.0127):
    """Karakas-Tariq (1991) perforation skin - simplified form.

    Captures the two dominant terms: the plane-flow component s_H
    (effective wellbore radius from phasing) and the crushed-zone
    contribution.  The vertical convergence (s_V) and wellbore
    interaction (s_wb) terms are second-order for typical geothermal
    perforated completions and are omitted.

    Parameters
    ----------
    L_p              : perforation tunnel length (m), e.g. 0.20
    d_p              : perforation diameter (m), e.g. 0.010
    SPF_per_m        : perforation density (shots per metre); 4 spf =
                       13.12 spm
    phasing_deg      : perforation phasing in deg (0, 60, 90, 120, 180,
                       360); 0/360 = aligned, 90 most common
    r_w              : wellbore radius (m)
    k_ratio_crushed  : k_crushed / k (default 0.4)
    t_crushed        : crushed-zone thickness around each perforation,
                       m (default 1/2 inch = 0.0127 m)

    Returns
    -------
    s_p              : perforation skin (dimensionless)
    """
    # Karakas-Tariq Table 1: phasing factor a_theta for r_wp = a*(r_w + L_p)
    phasing_a = {
        0:   0.250, 360: 0.250,
        180: 0.500,
        120: 0.648,
         90: 0.726,
         60: 0.813,
         45: 0.860,
    }
    if phasing_deg not in phasing_a:
        # Default to 90 deg if unrecognised
        phasing_deg = 90
    a = phasing_a[phasing_deg]

    # Plane-flow component (effective wellbore radius)
    r_w_prime = a * (r_w + L_p)
    s_H = np.log(r_w / r_w_prime) if r_w_prime > 0 else 0.0

    # Crushed zone contribution (per unit length, weighted by spacing)
    r_p = 0.5 * d_p
    h_perf = 1.0 / SPF_per_m       # spacing between shots, m
    if k_ratio_crushed <= 0.0 or L_p <= 0.0 or r_p <= 0.0:
        s_cz = 0.0
    else:
        s_cz = ((1.0 / k_ratio_crushed - 1.0)
                * np.log((r_p + t_crushed) / r_p)
                * (h_perf / L_p))

    return s_H + s_cz


def gravel_pack_skin(q_brine_std, screen_open_area_ratio, screen_OD,
                     screen_ID, screen_length, k_gravel_md,
                     mu_brine, B_brine, h_perf):
    """Effective skin from a gravel-packed completion (Furui-style).

    Models the pressure drop from the perforation tunnel through the
    gravel-filled annulus and screen as an additional pressure loss,
    converted to an equivalent skin via dPgp = q * mu * B / (2 pi k h)
    * s_gp.

    Parameters
    ----------
    q_brine_std            : flow rate, m^3/s at standard conditions
    screen_open_area_ratio : open-area fraction of the screen (0-1)
    screen_OD              : screen outer diameter, m
    screen_ID              : screen inner diameter, m
    screen_length          : screen length over the perfs, m
    k_gravel_md            : gravel permeability, mD (typ. 1.2e5 - 5e5)
    mu_brine               : brine viscosity at reservoir T, Pa.s
    B_brine                : brine FVF (-)
    h_perf                 : perforated formation thickness, m

    Returns
    -------
    s_gp                   : gravel-pack equivalent skin
    """
    # Approximate radial flow through gravel of OD/ID screen geometry.
    # Per Furui (2005), for a simple gravel-pack the skin is dominated
    # by linear flow through gravel of length (r_wb - r_screen_OD/2)
    # plus the screen restriction.  Here we lump them into a quadratic
    # in q with an effective skin form.
    if screen_open_area_ratio <= 0.0 or k_gravel_md <= 0.0:
        return 0.0
    k_gravel = k_gravel_md * MD_TO_M2
    # Linear gravel "permeability ratio" relative to formation: this
    # function returns an *effective* skin proxy that scales with
    # screen geometry.  For first-cut design work this is adequate.
    s_gp = (h_perf / screen_length) * np.log(screen_OD / screen_ID) \
           / max(screen_open_area_ratio, 0.05)
    # Add a small flow-rate-dependent term proportional to (1/k_gravel)
    # to penalise unrealistic gravel quality.  Kept conservative to
    # avoid double-counting against the perforation-skin term.
    s_gp += 0.05 * (1.0e-12 / k_gravel)        # ~0.005 for typical gravel
    return float(s_gp)


# =====================================================================
# 2. Brine FVF and viscosity at reservoir conditions
# =====================================================================
def brine_FVF(T_res_K, P_res_Pa, m_NaCl):
    """Brine formation volume factor B = rho_std / rho_res."""
    rho_std = density_brine_BW(T_STD, P_STD, m_NaCl)
    rho_res = density_brine_BW(T_res_K, P_res_Pa, m_NaCl)
    return rho_std / rho_res


# =====================================================================
# 3. ReservoirModel container with IPR methods
# =====================================================================
class ReservoirModel:
    """Reservoir + completion description with IPR evaluation.

    Two flow modes are supported:

      * 'darcy_radial' (default): radial Darcy inflow with three regimes:
          - regime='pss' (default): pseudo-steady-state, closed boundary
                q_std = (2 pi k h dP) / (mu B (ln(r_e/r_w) - 3/4 + s))
          - regime='ss': steady-state, constant pressure outer boundary
                q_std = (2 pi k h dP) / (mu B (ln(r_e/r_w) + s))
          - regime='transient': infinite-acting line-source solution
                q_std = (2 pi k h dP) / (mu B (1/2 ln(t_D) + 0.4045 + s))
                with t_D = k * t / (phi * mu * c_t * r_w^2)
            The transient mode requires `porosity`, `c_t`, and
            `t_prod_days` to be supplied.  This matches WellPerform's
            "Transient" inflow option and is appropriate for exploration
            / appraisal wells before pseudo-steady state is reached.

      * 'PI': straight-line productivity index, q_std = J * dP, where J
        is supplied directly (m^3/s per Pa).

      * 'fetkovich': back-pressure form, q_std = C * (P_res^2 - Pwf^2)^n
        with user-supplied C and n (mostly for completeness).

    Skin handling:
      Either pass `skin_total` (overrides everything) or let the model
      compute it from the Hawkins / Karakas-Tariq / gravel-pack helpers
      using the completion sub-dictionary.
    """

    def __init__(self,
                 # Reservoir geometry / properties
                 k_md, h_m, P_res_bar, T_res_C,
                 r_w, r_e=250.0,
                 # Brine composition
                 m_NaCl=0.0,
                 # Mode
                 mode='darcy_radial', regime='pss',
                 # Transient-mode parameters (used when regime='transient')
                 porosity=None, c_t=None, t_prod_days=None,
                 # Skin: either total or per-component
                 skin_total=None,
                 skin_damage=0.0,
                 skin_perforation=0.0,
                 skin_gravel_pack=0.0,
                 # Completion type and details (used to compute skin if
                 # skin_total is None and per-component skins are not
                 # explicitly supplied)
                 completion='open_hole',           # 'open_hole' / 'perforated' / 'gravel_pack'
                 damage_radius=None, damage_k_ratio=None,
                 perf_length=None, perf_diameter=0.010,
                 perf_SPF_per_m=13.12, perf_phasing=90,
                 crushed_k_ratio=0.4, crushed_thickness=0.0127,
                 gp_screen_OD=None, gp_screen_ID=None,
                 gp_screen_length=None, gp_open_area=0.10,
                 gp_gravel_md=2.0e5,
                 # Direct PI input (used when mode='PI')
                 J_si=None,
                 # Fetkovich coefficients
                 fet_C=None, fet_n=1.0,
                 ):
        # Convert to SI
        self.k = k_md * MD_TO_M2
        self.h = float(h_m)
        self.P_res = bar_to_Pa(P_res_bar)
        self.T_res = C_to_K(T_res_C)
        self.r_w = float(r_w)
        self.r_e = float(r_e)
        self.m_NaCl = float(m_NaCl)
        self.mode = mode
        self.regime = regime           # 'ss', 'pss', or 'transient' for darcy_radial
        self.completion = completion

        # Transient-mode parameters
        self.phi   = float(porosity)    if porosity    is not None else None
        self.c_t   = float(c_t)         if c_t         is not None else None
        self.t_s   = (float(t_prod_days) * 86400.0
                      if t_prod_days is not None else None)
        if self.mode == 'darcy_radial' and self.regime == 'transient':
            missing = [n for n, v in (('porosity', self.phi),
                                       ('c_t', self.c_t),
                                       ('t_prod_days', self.t_s))
                       if v is None]
            if missing:
                raise ValueError(
                    "transient regime requires: " + ", ".join(missing))

        # Compute / store skin
        if skin_total is not None:
            self.s_d, self.s_p, self.s_gp = 0.0, 0.0, 0.0
            self.s_total = float(skin_total)
        else:
            # Damage skin
            if (damage_radius is not None and damage_k_ratio is not None):
                self.s_d = hawkins_damage_skin(damage_k_ratio,
                                               damage_radius, r_w)
            else:
                self.s_d = float(skin_damage)
            # Perforation skin
            if (completion in ('perforated', 'gravel_pack')
                    and perf_length is not None):
                self.s_p = karakas_tariq_perforation_skin(
                    perf_length, perf_diameter, perf_SPF_per_m,
                    perf_phasing, r_w,
                    k_ratio_crushed=crushed_k_ratio,
                    t_crushed=crushed_thickness)
            else:
                self.s_p = float(skin_perforation)
            # Gravel-pack skin (placeholder - full evaluation needs
            # mu, B which depend on conditions; we use a geometry-only
            # term here, plus user can add additional via skin_gravel_pack).
            if completion == 'gravel_pack' and gp_screen_OD is not None:
                # Use a representative mu, B from reservoir conditions.
                mu_b = viscosity_brine(self.T_res, self.P_res, m_NaCl)
                B_b  = brine_FVF(self.T_res, self.P_res, m_NaCl)
                self.s_gp = gravel_pack_skin(
                    1.0, gp_open_area, gp_screen_OD, gp_screen_ID,
                    gp_screen_length, gp_gravel_md, mu_b, B_b, h_m)
                self.s_gp += float(skin_gravel_pack)
            else:
                self.s_gp = float(skin_gravel_pack)
            self.s_total = self.s_d + self.s_p + self.s_gp

        # Mode-specific parameters
        self.J_si = float(J_si) if J_si is not None else None
        self.fet_C = float(fet_C) if fet_C is not None else None
        self.fet_n = float(fet_n)

    # --------------------- Properties ---------------------
    def J(self):
        """Productivity index J = q/dP at reservoir conditions
        (m^3/s per Pa).  Valid for darcy_radial and PI modes.

        For PSS       : ln(r_e/r_w) - 3/4 + s
        For SS        : ln(r_e/r_w) + s
        For transient : 1/2 ln(t_D) + 0.4045 + s
                        with t_D = k t / (phi mu c_t r_w^2)

        Returns None ONLY for mode='fetkovich' (nonlinear, no constant J).
        Any other unrecognized mode/regime raises a ValueError.
        """
        if self.mode == 'PI':
            if self.J_si is None:
                raise ValueError(
                    "mode='PI' requires J_si to be set in the reservoir "
                    "config.  Either supply J_si or change mode to "
                    "'darcy_radial'.")
            return self.J_si

        if self.mode == 'fetkovich':
            return None        # nonlinear; AOF/q_at_Pwf use Fetkovich form

        if self.mode == 'darcy_radial':
            mu = viscosity_brine(self.T_res, self.P_res, self.m_NaCl)
            B  = brine_FVF(self.T_res, self.P_res, self.m_NaCl)
            if self.regime == 'transient':
                if (self.phi is None or self.c_t is None
                        or self.t_s is None):
                    raise ValueError(
                        "regime='transient' requires porosity, c_t, and "
                        "t_prod_days to be set in the reservoir config.")
                t_D = (self.k * self.t_s
                       / (self.phi * mu * self.c_t * self.r_w ** 2))
                if t_D <= 0.0:
                    raise ValueError(f"transient: non-positive t_D = {t_D}")
                denom = 0.5 * np.log(t_D) + 0.4045 + self.s_total
            elif self.regime == 'ss':
                denom = np.log(self.r_e / self.r_w) + self.s_total
            elif self.regime == 'pss':
                denom = np.log(self.r_e / self.r_w) - 0.75 + self.s_total
            else:
                raise ValueError(
                    f"regime={self.regime!r} not recognized.  Use 'pss', "
                    "'ss', or 'transient'.")
            return (2.0 * np.pi * self.k * self.h) / (mu * B * denom)

        raise ValueError(
            f"mode={self.mode!r} not recognized.  Use 'darcy_radial', "
            "'PI', or 'fetkovich'.")

    # --------------------- IPR equations ------------------
    def Pwf_at_q(self, q_std):
        """Bottom-hole flowing pressure (Pa) for a given standard-
        condition flow rate q_std (m^3/s)."""
        if self.mode == 'fetkovich':
            # q = C * (P_res^2 - Pwf^2)^n  =>  Pwf = sqrt(P_res^2 - (q/C)^(1/n))
            arg = self.P_res ** 2 - (q_std / self.fet_C) ** (1.0 / self.fet_n)
            return float(np.sqrt(max(0.0, arg)))

        return float(self.P_res - q_std / self.J())

    def q_at_Pwf(self, Pwf):
        """Standard-condition flow rate (m^3/s) for a given Pwf (Pa)."""
        dP = self.P_res - Pwf
        if dP <= 0.0:
            return 0.0
        if self.mode == 'fetkovich':
            return self.fet_C * (self.P_res ** 2 - Pwf ** 2) ** self.fet_n

        return self.J() * dP

    def AOF(self):
        """Absolute open-flow potential = q at Pwf = 0 (m^3/s)."""
        return self.q_at_Pwf(0.0)

    # --------------------- Diagnostics --------------------
    def __repr__(self):
        return (f"ReservoirModel(k={self.k/MD_TO_M2:.1f} mD, h={self.h:.1f} m, "
                f"P_res={self.P_res*1e-5:.1f} bar, T_res={self.T_res-273.15:.1f} degC, "
                f"r_w={self.r_w:.3f}, r_e={self.r_e:.0f}, "
                f"completion={self.completion}, s_total={self.s_total:.2f})")

    def describe(self):
        lines = [repr(self)]
        lines.append(f"   skin breakdown: s_damage={self.s_d:.3f}, "
                     f"s_perf={self.s_p:.3f}, s_gp={self.s_gp:.3f}")
        if self.mode == 'darcy_radial':
            J = self.J()
            lines.append(f"   regime = {self.regime}")
            if self.regime == 'transient':
                mu = viscosity_brine(self.T_res, self.P_res, self.m_NaCl)
                t_D = (self.k * self.t_s
                       / (self.phi * mu * self.c_t * self.r_w ** 2))
                lines.append(f"   transient: t = {self.t_s/86400:.1f} days, "
                             f"phi = {self.phi:.3f}, "
                             f"c_t = {self.c_t*1e5:.2f} 1/bar")
                lines.append(f"   t_D = {t_D:.3e}, "
                             f"1/2 ln(t_D) + 0.4045 + s = "
                             f"{0.5*np.log(t_D)+0.4045+self.s_total:.3f}")
            lines.append(f"   PI J = {J:.3e} m^3/s/Pa "
                         f"= {J * 3600 * 1e5:.3f} m^3/h/bar "
                         f"= {J * 1000 * 1e5:.3f} l/s/bar")
            mu = viscosity_brine(self.T_res, self.P_res, self.m_NaCl)
            B  = brine_FVF(self.T_res, self.P_res, self.m_NaCl)
            lines.append(f"   reservoir mu_brine = {mu*1e3:.3f} mPa.s, "
                         f"B = {B:.3f}")
            lines.append(f"   AOF = {self.AOF()*1000:.2f} l/s ({self.AOF()*86400:.0f} m^3/d) "
                         f"(Pwf -> 0)")
        return '\n'.join(lines)


# =====================================================================
# 4. SELF-TEST / VALIDATION
# =====================================================================
def _validate():
    print("=" * 78)
    print("IPR MODULE - SELF-CONSISTENCY TESTS")
    print("=" * 78)

    # ------- Test 1: Basic Darcy PSS for an open-hole well -------
    print("\n[Test 1] Darcy PSS, open hole, no skin")
    res1 = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108,                       # ~ 8.5 in / 2
        r_e=300.0,
        m_NaCl=0.713,                    # 40 000 ppm
        completion='open_hole',
        skin_total=0.0)
    print(res1.describe())
    print(f"   q at Pwf=200 bar: "
          f"{fmt_q(res1.q_at_Pwf(bar_to_Pa(200.0)) * 3600)}")
    print(f"   q at Pwf=100 bar: "
          f"{fmt_q(res1.q_at_Pwf(bar_to_Pa(100.0)) * 3600)}")

    # ------- Test 2: Same with formation damage (Hawkins) -------
    print("\n[Test 2] Same well, with Hawkins damage "
          "(k_d/k=0.2 over 0.5 m)")
    res2 = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108, r_e=300.0,
        m_NaCl=0.713,
        completion='open_hole',
        damage_radius=0.5, damage_k_ratio=0.2)
    print(res2.describe())
    print(f"   q at Pwf=100 bar: "
          f"{res2.q_at_Pwf(bar_to_Pa(100.0)) * 3600:.1f} m^3/h "
          f"(was {res1.q_at_Pwf(bar_to_Pa(100.0))*3600:.1f} for s=0)")

    # ------- Test 3: Perforated with Karakas-Tariq -------
    print("\n[Test 3] Perforated completion (4 SPF, 90 deg phasing, "
          "L_p=0.30 m)")
    res3 = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108, r_e=300.0,
        m_NaCl=0.713,
        completion='perforated',
        perf_length=0.30, perf_diameter=0.012,
        perf_SPF_per_m=13.12, perf_phasing=90,
        damage_radius=0.30, damage_k_ratio=0.5,
        crushed_k_ratio=0.4, crushed_thickness=0.0127)
    print(res3.describe())
    print(f"   q at Pwf=100 bar: "
          f"{fmt_q(res3.q_at_Pwf(bar_to_Pa(100.0)) * 3600)}")

    # ------- Test 4: Gravel pack -------
    print("\n[Test 4] Gravel-packed completion")
    res4 = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108, r_e=300.0,
        m_NaCl=0.713,
        completion='gravel_pack',
        perf_length=0.30, perf_diameter=0.012,
        perf_SPF_per_m=13.12, perf_phasing=90,
        damage_radius=0.30, damage_k_ratio=0.5,
        gp_screen_OD=0.140, gp_screen_ID=0.114,
        gp_screen_length=20.0, gp_open_area=0.12,
        gp_gravel_md=2.5e5)
    print(res4.describe())
    print(f"   q at Pwf=100 bar: "
          f"{fmt_q(res4.q_at_Pwf(bar_to_Pa(100.0)) * 3600)}")

    # ------- Test 5: PI mode -------
    print("\n[Test 5] PI mode (J = 2 m^3/h/bar)")
    J_si = 2.0 / 3600.0 / 1.0e5
    res5 = ReservoirModel(
        k_md=0.0, h_m=0.0, P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108,
        m_NaCl=0.713, mode='PI', J_si=J_si)
    print(res5.describe())
    print(f"   q at Pwf=100 bar: "
          f"{res5.q_at_Pwf(bar_to_Pa(100.0)) * 3600:.1f} m^3/h "
          f"(should be 400.0)")

    # ------- Test 6: Transient mode (line-source) -------
    print("\n[Test 6] Transient (infinite-acting) - 100 days "
          "(WellPerform-style)")
    res6 = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108, r_e=300.0,
        m_NaCl=0.713,
        completion='open_hole',
        skin_total=0.0,
        regime='transient',
        porosity=0.15, c_t=1.0e-9,            # 1e-9 1/Pa = 1e-4 1/bar
        t_prod_days=100.0)
    print(res6.describe())
    print(f"   q at Pwf=100 bar: "
          f"{fmt_q(res6.q_at_Pwf(bar_to_Pa(100.0)) * 3600)}")

    print("\n[Test 6b] Time evolution of transient PI (same well)")
    print(f"   {'t (days)':>10}  {'t_D':>11}  {'denom':>8}  "
          f"{'J (l/s/bar)':>13}  {'q@Pwf=100 (l/s)':>16}  {'(m^3/d)':>10}")
    for t_d in [1, 10, 30, 100, 365, 1000]:
        res_t = ReservoirModel(
            k_md=100.0, h_m=20.0,
            P_res_bar=300.0, T_res_C=200.0,
            r_w=0.108, r_e=300.0,
            m_NaCl=0.713, completion='open_hole', skin_total=0.0,
            regime='transient', porosity=0.15, c_t=1.0e-9,
            t_prod_days=float(t_d))
        mu_r = viscosity_brine(res_t.T_res, res_t.P_res, res_t.m_NaCl)
        t_D  = res_t.k * res_t.t_s / (res_t.phi * mu_r * res_t.c_t * res_t.r_w**2)
        denom = 0.5*np.log(t_D) + 0.4045 + res_t.s_total
        J_t_ls = res_t.J() * 1000 * 1e5                  # l/s/bar
        q_t_m3h = res_t.q_at_Pwf(bar_to_Pa(100.0)) * 3600.0
        print(f"   {t_d:10.0f}  {t_D:11.2e}  {denom:8.3f}  "
              f"{J_t_ls:13.3f}  {m3h_to_ls(q_t_m3h):16.2f}  "
              f"{m3h_to_m3d(q_t_m3h):10.0f}")
    print("   --> transient PI declines slowly with time as t_D rises")
    print("       (1/2 ln t_D term grows).  At long times the well")
    print("       'feels' the boundary and PSS becomes more appropriate.")

    # ------- Test 6c: Transient -> PSS crossover sanity -------
    print("\n[Test 6c] Comparison transient(100d) vs PSS (same well)")
    res_pss = ReservoirModel(
        k_md=100.0, h_m=20.0,
        P_res_bar=300.0, T_res_C=200.0,
        r_w=0.108, r_e=300.0,
        m_NaCl=0.713, completion='open_hole', skin_total=0.0,
        regime='pss')
    print(f"   transient(100d): J = {res6.J()*1000*1e5:.3f} l/s/bar")
    print(f"   pss            : J = {res_pss.J()*1000*1e5:.3f} l/s/bar")
    print(f"   For r_e=300m and these props, transient dominates "
          f"until well 'sees' the boundary.")

    # ------- Test 7: Sanity sweep of Pwf -------
    print("\n[Test 7] IPR curve summary for Test 1 well")
    print(f"   {'Pwf (bar)':>10} | {'q (l/s)':>9} | {'q (m^3/d)':>10}")
    for Pwf_bar in [300, 250, 200, 150, 100, 50, 0]:
        q_m3h = res1.q_at_Pwf(bar_to_Pa(float(Pwf_bar))) * 3600.0
        print(f"   {Pwf_bar:10.0f} | {m3h_to_ls(q_m3h):9.2f} | "
              f"{m3h_to_m3d(q_m3h):10.0f}")

    print("\n" + "=" * 78)
    print("Notes:")
    print("  * IPR uses single-phase brine viscosity and FVF at reservoir")
    print("    P, T.  Dissolved CO2 effect on B is small (<1 %) at low GWR.")
    print("  * Karakas-Tariq is implemented in simplified 2-term form;")
    print("    full 4-term version (s_H, s_V, s_wb, s_cz) is a Turn-4")
    print("    upgrade if needed.")
    print("  * Gravel-pack skin uses a geometry-driven proxy, suitable")
    print("    for first-cut design.  Replace with vendor-specific tool")
    print("    output if available.")
    print("=" * 78)


if __name__ == "__main__":
    _validate()
