"""
ipr_multilayer.py
=================

Commingled-production IPR for two or more layers producing into the same
wellbore.  Each layer has its own permeability, thickness, mid-perforation
depth, static pressure, temperature, and (optionally) salinity and skin.

The total surface rate at any common-Pwf condition (measured at the
*reference depth*, by default the deepest mid-perforation) is:

    q_total(Pwf_ref) = sum_i  q_i( Pwf_ref - dP_hydro_i )

where dP_hydro_i = rho_b_i * g * (z_ref - z_i) is the wellbore-hydrostatic
adjustment from the reference depth up to layer i's mid-perforation.
Layers physically above the reference depth (z_i < z_ref) get dP_hydro_i > 0,
so each individual layer "sees" a *higher* Pwf than the reference depth
(less drawdown), but the *deepest* layer is the reference and sees no
correction.

Internally, each layer holds its own `_LayerIPR` object which is just a
thin wrapper around the single-layer Darcy/PI/Fetkovich physics from
`ipr.py`.  The CommingledReservoir class exposes the same interface as
single-layer `ReservoirModel`, so the existing nodal solver and the
plotting routines can be used unchanged.

References:
  Brown, K.E. (1984). The Technology of Artificial Lift Methods, Vol 4.
  Beggs, H.D. (1991). Production Optimization, OGCI publications.
  Bahrami, H. et al. (2014). Commingled reservoirs: A review.
"""

from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Callable, List

from pvt import (
    T_STD, P_STD, GRAVITY,
    C_to_K, K_to_C, bar_to_Pa, Pa_to_bar,
    density_water_BW, density_brine_BW, viscosity_brine,
    flash_CO2_brine,
    fmt_q, m3h_to_ls, m3h_to_m3d,
    ppm_to_molality,
)

# Permeability conversion (1 mD in m^2)
MD_TO_M2 = 9.869233e-16


# =====================================================================
# 1. Single-layer Darcy / PI helper (lightweight version of ipr.py)
# =====================================================================
@dataclass
class Layer:
    """One reservoir layer in a commingled completion.

    Parameters
    ----------
    name : str
        Label, e.g. "sandstone" or "limestone".
    top_depth_m, bottom_depth_m : float
        Top and base of the perforated/completed interval, TVD in metres.
        Mid-perforation depth z_mid = 0.5*(top + bottom) is used as the
        layer's "pressure reference" for both its IPR (P_res is reported
        at z_mid) and for the wellbore-hydrostatic correction.
    h_net_m : float
        Net pay thickness (m).  May be less than (bottom - top) if some
        of the gross interval is non-reservoir.  For VGGT-1 sandstone:
        h_gross = 30.5 m but h_efektivno = 53 m  (effective vertical pay
        from log analysis).  Use h_net = h_ef.
    k_md : float
        Reservoir permeability (mD).
    P_res_bar : float
        Static reservoir pressure at z_mid (bar).
    T_res_C : float
        Reservoir temperature at z_mid (deg C).
    r_w : float
        Wellbore radius at the layer (m).  Usually the same for all layers.
    r_e : float
        Drainage radius (m).  Default 300 m for geothermal wells.
    NaCl_ppm : float
        Salinity in ppm.  Each layer may have its own brine composition.
    mode : {'darcy_radial', 'PI', 'fetkovich'}
        Inflow form.  Default 'darcy_radial'.
    regime : {'pss', 'ss', 'transient'}
        Boundary regime for 'darcy_radial'.
    skin_total : float
        Total skin factor (-).  Default 0.
    porosity, c_t, t_prod_days : Optional[float]
        Required only when regime='transient'.
    J_si : Optional[float]
        Productivity index (m^3/s/Pa).  Required only when mode='PI'.
    fet_C, fet_n : Optional[float]
        Required only when mode='fetkovich'.
    """
    name: str
    top_depth_m: float
    bottom_depth_m: float
    h_net_m: float
    k_md: float
    P_res_bar: float
    T_res_C: float
    r_w: float = 0.108
    r_e: float = 300.0
    NaCl_ppm: float = 0.0
    mode: str = 'darcy_radial'
    regime: str = 'pss'
    skin_total: float = 0.0
    porosity: Optional[float] = None
    c_t: Optional[float] = None
    t_prod_days: Optional[float] = None
    J_si: Optional[float] = None
    fet_C: Optional[float] = None
    fet_n: Optional[float] = None
    # Forchheimer / non-Darcy rate-dependent skin coefficient (1/(m^3/s)).
    # If non-zero, the effective skin used in the radial-Darcy inflow is
    #   s_eff = skin_total + D_nonDarcy * q_layer,
    # making q_layer the positive root of a quadratic in q (see q_at_Pwf).
    # Only used when mode == 'darcy_radial'.  Default 0.0 = linear PI.
    D_nonDarcy: float = 0.0

    def __post_init__(self):
        self.z_mid_m = 0.5 * (self.top_depth_m + self.bottom_depth_m)
        # Convert to SI internally
        self.k     = self.k_md * MD_TO_M2          # m^2
        self.P_res = bar_to_Pa(self.P_res_bar)     # Pa
        self.T_res = C_to_K(self.T_res_C)          # K
        self.m_NaCl = ppm_to_molality(self.NaCl_ppm)
        if self.t_prod_days is not None:
            self.t_s = self.t_prod_days * 86400.0
        else:
            self.t_s = None
        if self.mode == 'darcy_radial' and self.regime == 'transient':
            missing = [n for n, v in (('porosity', self.porosity),
                                       ('c_t', self.c_t),
                                       ('t_prod_days', self.t_s))
                       if v is None]
            if missing:
                raise ValueError(
                    f"Layer '{self.name}': transient regime requires: "
                    + ", ".join(missing))

    # ----- Fluid properties at reservoir conditions -----
    def mu_brine(self):
        """Brine viscosity at the layer's P_res, T_res (Pa.s)."""
        return viscosity_brine(self.T_res, self.P_res, self.m_NaCl)

    def B_brine(self):
        """Formation-volume factor (-) = rho_std / rho_res."""
        rho_std = density_brine_BW(T_STD, P_STD, self.m_NaCl)
        rho_res = density_brine_BW(self.T_res, self.P_res, self.m_NaCl)
        return rho_std / rho_res

    def rho_brine_res(self):
        """Brine density at reservoir conditions (kg/m^3)."""
        return density_brine_BW(self.T_res, self.P_res, self.m_NaCl)

    # ----- Productivity index -----
    def J(self):
        """Productivity index J = dq/dP at reservoir conditions
        (m^3/s/Pa).  Returns None for Fetkovich (nonlinear)."""
        if self.mode == 'PI':
            if self.J_si is None:
                raise ValueError(
                    f"Layer '{self.name}': mode='PI' requires J_si.")
            return self.J_si

        if self.mode == 'fetkovich':
            return None

        if self.mode == 'darcy_radial':
            mu = self.mu_brine()
            B  = self.B_brine()
            if self.regime == 'transient':
                t_D = (self.k * self.t_s
                       / (self.porosity * mu * self.c_t * self.r_w ** 2))
                if t_D <= 0.0:
                    raise ValueError(
                        f"Layer '{self.name}': non-positive t_D = {t_D}")
                denom = 0.5 * np.log(t_D) + 0.4045 + self.skin_total
            elif self.regime == 'ss':
                denom = np.log(self.r_e / self.r_w) + self.skin_total
            elif self.regime == 'pss':
                denom = (np.log(self.r_e / self.r_w) - 0.75
                          + self.skin_total)
            else:
                raise ValueError(
                    f"Layer '{self.name}': unknown regime "
                    f"{self.regime!r}.  Use 'pss', 'ss', 'transient'.")
            return (2.0 * np.pi * self.k * self.h_net_m) / (mu * B * denom)

        raise ValueError(
            f"Layer '{self.name}': unknown mode {self.mode!r}.")

    def q_at_Pwf(self, Pwf_layer):
        """Standard-condition flow rate (m^3/s) at the given Pwf at the
        layer's mid-perf depth (Pa).  Returns 0 if Pwf >= P_res.

        If D_nonDarcy != 0 and mode == 'darcy_radial', the layer rate is
        the positive root of the quadratic

            (D / denom0) * q^2  +  q  -  J0 * (P_res - Pwf)  =  0,

        where denom0 = ln(re/rw) - 0.75 + s_total  (pss regime), and
        J0 is the linear-PI productivity index at zero rate.  This is
        the standard Forchheimer rate-dependent skin formulation
        (Jones et al., 1976; Golan & Whitson, 1991).
        """
        dP = self.P_res - Pwf_layer
        if dP <= 0.0:
            return 0.0
        if self.mode == 'fetkovich':
            return self.fet_C * (self.P_res ** 2 - Pwf_layer ** 2) ** self.fet_n

        if (self.mode == 'darcy_radial'
                and getattr(self, 'D_nonDarcy', 0.0) != 0.0):
            # denom0 for the steady/pss regime; for transient use the
            # log-time-based denom from J() consistently.
            if self.regime == 'pss':
                denom0 = np.log(self.r_e / self.r_w) - 0.75 + self.skin_total
            elif self.regime == 'ss':
                denom0 = np.log(self.r_e / self.r_w) + self.skin_total
            else:
                mu = self.mu_brine()
                t_D = (self.k * self.t_s
                       / (self.porosity * mu * self.c_t * self.r_w ** 2))
                denom0 = 0.5 * np.log(t_D) + 0.4045 + self.skin_total
            J0 = self.J()
            a = self.D_nonDarcy / max(denom0, 1e-6)
            b = 1.0
            c = -J0 * dP
            disc = b * b - 4.0 * a * c
            if disc < 0.0:
                return 0.0
            if a <= 0.0:
                return max(0.0, -c / b)
            q = (-b + np.sqrt(disc)) / (2.0 * a)
            return max(0.0, q)

        return self.J() * dP

    def q_at_Pwf_signed(self, Pwf_layer):
        """Signed layer rate (m^3/s), allowing INJECTION (negative q) when
        Pwf_layer > P_res.

        This is the "crossflow" variant of q_at_Pwf -- it does NOT clip
        at zero.  Use it for shut-in / crossflow analysis and for
        building the per-layer IPR plot that shows the breakpoints at
        each layer's P_res.

        Sign convention: q > 0 = production (formation -> well),
                         q < 0 = injection / crossflow (well -> formation).

        Forchheimer is APPLIED in the production direction only.  When
        the layer is taking injection (Pwf > P_res), it is assumed to
        do so linearly (rate-dependent skin in the negative direction
        is poorly defined and typically negligible at these low
        crossflow rates).
        """
        dP = self.P_res - Pwf_layer
        if self.mode == 'fetkovich':
            # Symmetric Fetkovich (not standard; treat same as production)
            if dP > 0:
                return self.fet_C * (self.P_res ** 2 - Pwf_layer ** 2) ** self.fet_n
            return -self.fet_C * (Pwf_layer ** 2 - self.P_res ** 2) ** self.fet_n
        if dP >= 0:
            return self.q_at_Pwf(Pwf_layer)
        # Pwf > P_res -> linear injection at the same J
        return self.J() * dP    # negative

    def Pwf_at_q(self, q_std):
        """Pwf (Pa) at layer mid-perf for a given standard rate (m^3/s)."""
        if self.mode == 'fetkovich':
            inside = (self.P_res ** 2
                      - (q_std / self.fet_C) ** (1.0 / self.fet_n))
            if inside <= 0:
                return 0.0
            return float(np.sqrt(inside))

        if (self.mode == 'darcy_radial'
                and getattr(self, 'D_nonDarcy', 0.0) != 0.0):
            # Re-arrange quadratic for dP given q:
            #   J0*dP = q + (D/denom0)*q^2
            #   dP   = (q + (D/denom0)*q^2) / J0
            if self.regime == 'pss':
                denom0 = np.log(self.r_e / self.r_w) - 0.75 + self.skin_total
            elif self.regime == 'ss':
                denom0 = np.log(self.r_e / self.r_w) + self.skin_total
            else:
                mu = self.mu_brine()
                t_D = (self.k * self.t_s
                       / (self.porosity * mu * self.c_t * self.r_w ** 2))
                denom0 = 0.5 * np.log(t_D) + 0.4045 + self.skin_total
            J0 = self.J()
            dP = (q_std + (self.D_nonDarcy / max(denom0, 1e-6)) * q_std**2) / J0
            return float(self.P_res - dP)

        return float(self.P_res - q_std / self.J())

    def AOF(self):
        """Absolute open flow potential at the layer's mid-perf
        (q at Pwf_layer = 0)  (m^3/s)."""
        return self.q_at_Pwf(0.0)

    def describe(self):
        """Pretty-printed multi-line summary of this layer."""
        J_si = self.J() if self.mode != 'fetkovich' else None
        lines = [
            f"  Layer '{self.name}':",
            f"    interval     : {self.top_depth_m:.1f} - "
                f"{self.bottom_depth_m:.1f} m  (mid = {self.z_mid_m:.1f} m)",
            f"    h_net        : {self.h_net_m:.1f} m",
            f"    k            : {self.k_md:.2f} mD  "
                f"(k*h = {self.k_md * self.h_net_m:.0f} mD.m)",
            f"    P_res        : {self.P_res_bar:.2f} bar @ z_mid",
            f"    T_res        : {self.T_res_C:.2f} degC",
            f"    NaCl         : {self.NaCl_ppm:.0f} ppm  "
                f"(m = {self.m_NaCl:.3f} mol/kg)",
            f"    mode/regime  : {self.mode}/{self.regime}",
            f"    skin_total   : {self.skin_total:.3f}",
        ]
        if getattr(self, 'D_nonDarcy', 0.0) != 0.0:
            lines.append(
                f"    D_nonDarcy   : {self.D_nonDarcy:.3e} (m^3/s)^-1"
                f"  ({self.D_nonDarcy/86400.0:.3e} (m^3/d)^-1)")
        if J_si is not None:
            lines.append(
                f"    J            : {J_si:.3e} m^3/s/Pa  "
                f"= {J_si * 1000 * 1e5:.3f} l/s/bar  "
                f"= {J_si * 3600 * 1e5:.3f} m^3/h/bar")
            AOF_m3h = self.AOF() * 3600.0
            lines.append(
                f"    AOF          : {fmt_q(AOF_m3h)}")
        mu = self.mu_brine()
        B  = self.B_brine()
        rho = self.rho_brine_res()
        lines.append(
            f"    mu_b,res     : {mu*1e3:.3f} mPa.s ({mu:.3e} Pa.s)")
        lines.append(
            f"    B            : {B:.4f}")
        lines.append(
            f"    rho_b,res    : {rho:.1f} kg/m^3")
        return '\n'.join(lines)


# =====================================================================
# 2. Commingled reservoir
# =====================================================================
class CommingledReservoir:
    """Two or more layers producing through a common wellbore.

    The class exposes the *same interface* as single-layer
    `ReservoirModel`, so the existing nodal solver, plotter and
    sensitivity sweeps work unchanged on commingled reservoirs.

    Parameters
    ----------
    layers : list of Layer
        Layers, in any order.
    reference_depth_m : Optional[float]
        TVD at which the commingled (q, Pwf) IPR is evaluated.  Defaults
        to the deepest layer's mid-perforation.  This is the node depth
        for nodal analysis; the VLP profile must be marched to this
        same depth.
    wellbore_density_kg_m3 : Optional[float]
        Wellbore brine density used to compute inter-layer hydrostatic
        offsets.  Defaults to the density of the deepest layer's brine
        at its reservoir conditions (a good approximation at modest
        rates when CO2 stays dissolved).  Pass a different value if
        you have a better estimate (e.g. measured static gradient or
        average of layers).
    """
    def __init__(self, layers, reference_depth_m=None,
                 wellbore_density_kg_m3=None):
        if not layers:
            raise ValueError("CommingledReservoir needs at least 1 layer.")
        self.layers = list(layers)
        # Reference depth = deepest layer's mid-perf by default
        if reference_depth_m is None:
            reference_depth_m = max(L.z_mid_m for L in self.layers)
        self.z_ref = float(reference_depth_m)
        # Wellbore density for inter-layer hydrostatic
        if wellbore_density_kg_m3 is None:
            # Use deepest layer's brine density at reservoir conditions
            deepest = max(self.layers, key=lambda L: L.z_mid_m)
            self.rho_wb = deepest.rho_brine_res()
        else:
            self.rho_wb = float(wellbore_density_kg_m3)

        # ------------------------------------------------------------------
        # Exposed attributes matching the single-layer ReservoirModel API
        # ------------------------------------------------------------------
        # P_res "of the commingled system" -- shut-in Pwf at the reference
        # depth (i.e. q=0).  At shut-in each layer is at its own static
        # P_res; the wellbore communicates via hydrostatic.  The static
        # Pwf at z_ref is well-defined only if all layers are in
        # hydrostatic equilibrium with each other -- this is an idealised
        # value.  We use the deepest layer's P_res as the canonical
        # "P_res" (since it lives at z_ref by default).
        deepest = max(self.layers, key=lambda L: L.z_mid_m)
        self.P_res = deepest.P_res
        self.T_res = deepest.T_res
        self.r_w   = deepest.r_w
        self.r_e   = deepest.r_e
        self.m_NaCl = deepest.m_NaCl
        self.mode  = 'commingled'
        self.regime = 'commingled'
        self.completion = 'commingled'
        # Aggregate "effective" skin = kh-weighted average
        kh_sum = sum(L.k_md * L.h_net_m for L in self.layers)
        if kh_sum > 0:
            self.s_total = (
                sum(L.skin_total * L.k_md * L.h_net_m for L in self.layers)
                / kh_sum)
        else:
            self.s_total = 0.0
        # Per-component skin breakdown not tracked at commingled level
        self.s_d = 0.0
        self.s_p = 0.0
        self.s_gp = 0.0
        # PI-mode field for compatibility (None means we use commingled
        # piecewise-linear evaluation, not a single J)
        self.J_si = None
        self.fet_C = None
        self.fet_n = None
        # Transient fields kept None at commingled level
        self.phi = None
        self.c_t = None
        self.t_s = None

    # ------------------------------------------------------------------
    # Inter-layer hydrostatic offset
    # ------------------------------------------------------------------
    def _dP_hydro_to_layer(self, layer):
        """Wellbore-hydrostatic pressure drop (Pa) from the reference
        depth UP to the given layer's mid-perforation.

        Convention: positive value if z_ref > z_layer (layer is *above*
        reference depth), in which case the Pwf at the shallower layer
        is *less* than at the reference depth by this amount.

        Sign convention is consistent with:
            Pwf_layer = Pwf_ref - dP_hydro_to_layer
        """
        dz = self.z_ref - layer.z_mid_m         # positive if layer above ref
        return self.rho_wb * GRAVITY * dz

    # ------------------------------------------------------------------
    # Inflow / outflow at the reference depth
    # ------------------------------------------------------------------
    def q_at_Pwf(self, Pwf_ref):
        """Total standard-condition rate (m^3/s) for a given Pwf
        (Pa) at the reference depth.  Sums each layer's individual
        contribution accounting for inter-layer hydrostatic offset."""
        q_total = 0.0
        for L in self.layers:
            Pwf_L = Pwf_ref - self._dP_hydro_to_layer(L)
            q_total += L.q_at_Pwf(Pwf_L)
        return q_total

    def Pwf_at_q(self, q_total):
        """Invert q_at_Pwf to find the Pwf (Pa) at the reference depth
        that gives the requested total rate.  Solved by bisection
        between Pwf=0 and Pwf_max = max(P_res_i + dP_hydro_i)."""
        from scipy.optimize import brentq
        # Upper bracket: max Pwf at which any layer is still injecting
        # (i.e. q would be zero or negative); add a small margin.
        Pwf_hi = max(L.P_res + self._dP_hydro_to_layer(L)
                     for L in self.layers)
        f_lo = self.q_at_Pwf(0.0) - q_total
        f_hi = self.q_at_Pwf(Pwf_hi) - q_total
        if f_lo * f_hi > 0:
            # q_total outside the achievable range
            if q_total > self.AOF():
                return 0.0   # capped at AOF
            return Pwf_hi
        return float(brentq(lambda P: self.q_at_Pwf(P) - q_total,
                             0.0, Pwf_hi))

    def AOF(self):
        """Total absolute open-flow potential at Pwf_ref = 0."""
        return self.q_at_Pwf(0.0)

    def J(self):
        """Effective productivity index = sum of layer Js (m^3/s/Pa).
        For a strictly linear (Darcy) commingled IPR, this is the
        slope of q vs Pwf_ref *near small drawdowns* (i.e. when no
        layer has shut off).  Returns None if any layer is Fetkovich."""
        total = 0.0
        for L in self.layers:
            J = L.J()
            if J is None:
                return None
            total += J
        return total

    # ------------------------------------------------------------------
    # Per-layer breakdown at a given Pwf
    # ------------------------------------------------------------------
    def layer_rates_at_Pwf(self, Pwf_ref):
        """Returns list of (layer_name, q_std_m3_s) for each layer at
        the given Pwf_ref."""
        out = []
        for L in self.layers:
            Pwf_L = Pwf_ref - self._dP_hydro_to_layer(L)
            out.append((L.name, L.q_at_Pwf(Pwf_L)))
        return out

    def layer_rates_at_Pwf_signed(self, Pwf_ref):
        """Like layer_rates_at_Pwf, but returns SIGNED rates: positive
        for layers producing, negative for layers injecting (crossflow).

        Use this for the crossflow analysis and for the IPR plot that
        walks each layer's P_res breakpoint."""
        out = []
        for L in self.layers:
            Pwf_L = Pwf_ref - self._dP_hydro_to_layer(L)
            out.append((L.name, L.q_at_Pwf_signed(Pwf_L)))
        return out

    def q_total_signed(self, Pwf_ref):
        """Net rate at Pwf_ref summing SIGNED layer contributions
        (production minus crossflow injection).  Equals the surface
        rate the well actually delivers under steady-state crossflow.
        """
        return sum(q for _, q in self.layer_rates_at_Pwf_signed(Pwf_ref))

    def describe(self):
        lines = [
            f"CommingledReservoir: {len(self.layers)} layers",
            f"  reference depth     : {self.z_ref:.1f} m TVD",
            f"  wellbore brine rho  : {self.rho_wb:.1f} kg/m^3  "
                "(used for inter-layer hydrostatic)",
            f"  total AOF           : {fmt_q(self.AOF() * 3600.0)}",
        ]
        J = self.J()
        if J is not None:
            lines.append(
                f"  effective J         : "
                f"{J * 1000 * 1e5:.3f} l/s/bar  "
                f"({J * 3600 * 1e5:.3f} m^3/h/bar)  "
                f"(slope at low draw-down)")
        lines.append("")
        for L in self.layers:
            lines.append(L.describe())
            dP_h = self._dP_hydro_to_layer(L)
            lines.append(
                f"    dP_hydro to ref: "
                f"{dP_h*1e-5:+.2f} bar  "
                f"(layer is {abs(self.z_ref - L.z_mid_m):.1f} m "
                f"{'above' if L.z_mid_m < self.z_ref else 'below'} "
                "reference)")
            lines.append("")
        return '\n'.join(lines)


# =====================================================================
# 3. Self-test on simple cases
# =====================================================================
def _validate():
    """Sanity tests on commingled IPR."""
    print("=" * 78)
    print("IPR MULTILAYER MODULE - SELF-CONSISTENCY TESTS")
    print("=" * 78)

    # ------- Test 1: degenerate single-layer case must match single-layer -------
    print("\n[Test 1] Single layer in CommingledReservoir matches stand-alone")
    L1 = Layer(
        name='only', top_depth_m=1890.0, bottom_depth_m=1920.0,
        h_net_m=30.0, k_md=70.0, P_res_bar=194.0, T_res_C=100.0,
        r_w=0.108, r_e=300.0, NaCl_ppm=5000.0, regime='pss',
        skin_total=0.0)
    com = CommingledReservoir([L1])
    J_layer = L1.J()
    J_com   = com.J()
    AOF_layer = L1.AOF()
    AOF_com   = com.AOF()
    print(f"  J  (Layer)     = {J_layer*1000*1e5:.4f} l/s/bar")
    print(f"  J  (Commingled)= {J_com*1000*1e5:.4f} l/s/bar")
    print(f"  AOF (Layer)    = {fmt_q(AOF_layer*3600)}")
    print(f"  AOF (Commingled)= {fmt_q(AOF_com*3600)}")
    assert abs(J_com - J_layer) < 1e-15
    assert abs(AOF_com - AOF_layer) < 1e-10
    print("  PASS - identical results.")

    # ------- Test 2: VGGT-1 two-layer case -------
    print("\n[Test 2] VGGT-1 commingled (sandstone + limestone)")

    # Sandstone -- exact PBU-derived values from the report:
    #   k = 71.29 mD, h_ef = 53 m, P_res = 194.05 bar @ 1888 m,
    #   T_res = 99.73 C, phi = 0.237, skin = -0.135.
    # Salinity: final salinity from the 3-stage test = 4.95 g/L; convert
    # to ppm by mass ~ 4.95 / 1.00 = 4950 ppm  (low salinity).
    sand = Layer(
        name='sandstone',
        top_depth_m=1889.5, bottom_depth_m=1920.0,
        h_net_m=53.0,        # h_efektivno from PBU report
        k_md=71.29,
        P_res_bar=194.05,
        T_res_C=99.73,
        r_w=0.108, r_e=300.0,
        NaCl_ppm=4950.0,
        regime='pss',
        skin_total=-0.135)

    # Limestone -- k is uncertain (mini-frak indicated "weak permeability",
    # acid jobs improved it).  Use a screening k from kh observed after
    # 2nd acid job: ~2.7 l/s from the natural-flow test, modelled later.
    # For now, plug k_md = 5 mD (representative of the matrix after acid)
    # and let the user override.  P_res = 201.5 bar @ 2280 m.  We need
    # P_res at z_mid = 2211.5 m: extrapolate up by ~ rho*g*68.5 m.  Using
    # rho ~ 980 kg/m^3 and g = 9.81: dP ~ 6.6 bar, so P_res ~ 194.9 bar
    # at z_mid.  T_res at z_mid: 113.9 C at 2280 m, gradient ~0.04 C/m,
    # so ~111.2 C at 2211.5 m.
    lime = Layer(
        name='limestone',
        top_depth_m=2144.0, bottom_depth_m=2279.0,
        h_net_m=135.0,
        k_md=5.0,            # representative post-acid; user can update
        P_res_bar=194.9,     # extrapolated to z_mid = 2211.5 m
        T_res_C=111.2,
        r_w=0.108, r_e=300.0,
        NaCl_ppm=28000.0,    # salinity after acid washout
        regime='pss',
        skin_total=0.0)

    com = CommingledReservoir([sand, lime])
    print(com.describe())

    # Layer-by-layer at AOF (Pwf_ref = 0)
    print(f"\n  Layer-by-layer breakdown at Pwf_ref = 0 (AOF condition):")
    for name, q in com.layer_rates_at_Pwf(0.0):
        print(f"    {name:<12s}: {fmt_q(q*3600.0)}")

    # And at a 30-bar drawdown (realistic operating point)
    Pwf_test = com.P_res - bar_to_Pa(30.0)
    print(f"\n  Layer-by-layer breakdown at Pwf_ref = "
          f"{Pwf_test*1e-5:.1f} bar (30 bar drawdown):")
    for name, q in com.layer_rates_at_Pwf(Pwf_test):
        print(f"    {name:<12s}: {fmt_q(q*3600.0)}")
    print(f"    total       : {fmt_q(com.q_at_Pwf(Pwf_test)*3600.0)}")

    # ------- Test 3: IPR curve sweep -------
    print("\n[Test 3] Commingled IPR curve at the reference depth")
    print(f"  {'Pwf_ref (bar)':>14} | {'q_total (l/s)':>14} | "
          f"{'q_sand (l/s)':>13} | {'q_lime (l/s)':>13}")
    for Pwf_bar in [200, 190, 180, 170, 150, 100, 50, 0]:
        Pwf_Pa = bar_to_Pa(float(Pwf_bar))
        q_tot = com.q_at_Pwf(Pwf_Pa) * 3600.0
        rates = dict(com.layer_rates_at_Pwf(Pwf_Pa))
        q_s = rates['sandstone'] * 3600.0
        q_l = rates['limestone'] * 3600.0
        print(f"  {Pwf_bar:14.0f} | {m3h_to_ls(q_tot):14.2f} | "
              f"{m3h_to_ls(q_s):13.2f} | {m3h_to_ls(q_l):13.2f}")

    # ------- Test 4: inversion consistency -------
    print("\n[Test 4] q -> Pwf -> q round-trip")
    q_target_ls = 23.0   # ESP test average
    q_target_si = q_target_ls / 1000.0
    Pwf = com.Pwf_at_q(q_target_si)
    q_back = com.q_at_Pwf(Pwf) * 1000.0    # back to l/s
    print(f"  Target q     : {q_target_ls:.2f} l/s")
    print(f"  -> Pwf_ref   : {Pwf*1e-5:.2f} bar  "
          f"(drawdown = {(com.P_res - Pwf)*1e-5:.2f} bar)")
    print(f"  -> q (back)  : {q_back:.4f} l/s")
    print(f"  Round-trip OK: {abs(q_back - q_target_ls) < 0.01}")


if __name__ == "__main__":
    _validate()
