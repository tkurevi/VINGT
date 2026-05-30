"""
pvt.py
======

PVT (Pressure-Volume-Temperature) module for the CO2-brine system in
geothermal wells.

Provides thermodynamic and transport properties for:
  * Pure water and NaCl brine (density, viscosity, surface tension)
  * Pure CO2 (density via Peng-Robinson EOS, viscosity via Fenghour et al.)
  * CO2 solubility in NaCl brine (Duan & Sun 2003)
  * "Live" brine - brine + dissolved CO2 - density via Garcia (2001)
    apparent molar volume of CO2
  * A two-phase equilibrium flash that, given total CO2 loading and local
    (P, T, salinity), returns phase split and per-phase properties
  * Standard-condition gas-water ratio (GWR) conversions

Units (SI used internally everywhere):
  Pressure          P : Pa
  Temperature       T : K
  Density           rho : kg/m^3
  Viscosity         mu : Pa.s
  NaCl salinity     m_NaCl : mol NaCl per kg of H2O (molality)
  Mass fraction     w : kg solute per kg of mixture
  Surface tension   sigma : N/m
  Molar volume      V : m^3/mol

Convenience converters between common engineering units (bar, MPa, degC,
ppm) are provided at the top of the module.

Validation references for the test block at the end of the file are taken
from the NIST Webbook (CO2) and IAPWS Industrial Formulation (water).

References
----------
Batzle, M. and Wang, Z. (1992). Seismic properties of pore fluids.
    Geophysics 57(11), 1396-1408.
Mao, S. and Duan, Z. (2009). The viscosity of aqueous alkali-chloride
    solutions up to 623 K, 1000 bar, and high ionic strength.
    Int. J. Thermophys. 30, 1510-1523.
Peng, D.Y. and Robinson, D.B. (1976). A new two-constant equation of
    state. Ind. Eng. Chem. Fundam. 15, 59-64.
Fenghour, A., Wakeham, W.A., Vesovic, V. (1998). The viscosity of carbon
    dioxide. J. Phys. Chem. Ref. Data 27, 31-44.
Duan, Z. and Sun, R. (2003). An improved model calculating CO2 solubility
    in pure water and aqueous NaCl solutions from 273 to 533 K and from
    0 to 2000 bar. Chem. Geol. 193, 257-271.
Garcia, J.E. (2001). Density of aqueous solutions of CO2.
    Lawrence Berkeley National Lab. report LBNL-49023.
Wagner, W. and Pruss, A. (1993). International equations for the
    saturation properties of ordinary water substance. J. Phys. Chem.
    Ref. Data 22, 783.
IAPWS R1-76(2014). Revised release on surface tension of ordinary
    water substance.
"""

from __future__ import annotations
import numpy as np

# =====================================================================
# Physical constants and molar masses
# =====================================================================
R_GAS  = 8.314472         # J / (mol.K)  - universal gas constant
M_CO2  = 44.0098e-3       # kg / mol     - CO2 molar mass
M_H2O  = 18.01528e-3      # kg / mol     - water molar mass
M_NACL = 58.4428e-3       # kg / mol     - NaCl molar mass

# CO2 critical properties (NIST)
TC_CO2    = 304.1282      # K
PC_CO2    = 7.3773e6      # Pa
OMEGA_CO2 = 0.22394       # acentric factor

# Water critical properties (IAPWS)
TC_H2O    = 647.096       # K
PC_H2O    = 22.064e6      # Pa

# Standard conditions used for GWR conversions (15 degC, 1 atm).
# Common in petroleum / geothermal volumetric reporting.
T_STD = 288.15            # K  (15 degC)
P_STD = 101325.0          # Pa (1 atm)
GRAVITY = 9.81            # m/s^2  (gravitational acceleration)


# =====================================================================
# Unit converters
# =====================================================================
def C_to_K(T_C):
    return T_C + 273.15


def K_to_C(T_K):
    return T_K - 273.15


def bar_to_Pa(P_bar):
    return P_bar * 1.0e5


def Pa_to_bar(P_Pa):
    return P_Pa / 1.0e5


def MPa_to_Pa(P_MPa):
    return P_MPa * 1.0e6


def Pa_to_MPa(P_Pa):
    return P_Pa / 1.0e6


# ---------------------------------------------------------------------
# Flow-rate unit conversions (geothermal/water-well convention)
# Internal data is in m^3/h; user-facing display is l/s with m^3/d in
# brackets, e.g. "1.0 l/s (86.4 m^3/d)".
# ---------------------------------------------------------------------
def m3h_to_ls(q_m3h):
    """m^3/h -> litres per second.  1 m^3/h = 0.27778 l/s."""
    return q_m3h / 3.6


def m3h_to_m3d(q_m3h):
    """m^3/h -> m^3/day.  1 m^3/h = 24 m^3/d."""
    return q_m3h * 24.0


def ls_to_m3h(q_ls):
    """litres/second -> m^3/h."""
    return q_ls * 3.6


def m3d_to_m3h(q_m3d):
    """m^3/day -> m^3/h."""
    return q_m3d / 24.0


def fmt_q(q_m3h, ls_dec=2, m3d_dec=1):
    """Format a flow rate in m^3/h as 'X.XX l/s (Y.Y m^3/d)'."""
    q_ls  = m3h_to_ls(q_m3h)
    q_m3d = m3h_to_m3d(q_m3h)
    return f"{q_ls:.{ls_dec}f} l/s ({q_m3d:.{m3d_dec}f} m^3/d)"


def ppm_to_molality(ppm_NaCl):
    """NaCl mass concentration in ppm (mg/kg of solution) -> molality
    (mol NaCl per kg of water)."""
    w = ppm_NaCl * 1.0e-6
    if w >= 1.0:
        raise ValueError("Salinity mass fraction must be < 1.")
    return w / (M_NACL * (1.0 - w))


def molality_to_mass_fraction(m_NaCl):
    """NaCl molality -> NaCl mass fraction (kg salt / kg solution)."""
    return m_NaCl * M_NACL / (1.0 + m_NaCl * M_NACL)


def mass_fraction_to_molality(w_NaCl):
    """NaCl mass fraction (kg salt / kg solution) -> molality."""
    if w_NaCl >= 1.0:
        raise ValueError("Mass fraction must be < 1.")
    return w_NaCl / (M_NACL * (1.0 - w_NaCl))


# =====================================================================
# 1. WATER SATURATION PRESSURE (Wagner & Pruss 1993, IAPWS-IF97 form)
# =====================================================================
def saturation_pressure_water(T_K):
    """Pure water saturation pressure (Pa) - Wagner & Pruss (1993).
    Valid 273.16 K - 647.096 K, accuracy <0.05 %."""
    if T_K >= TC_H2O:
        return PC_H2O
    if T_K < 250.0:
        T_K = 250.0  # very cold - clamp to avoid numerical issues
    tau = 1.0 - T_K / TC_H2O
    a = (-7.85951783, 1.84408259, -11.7866497,
         22.6807411, -15.9618719, 1.80122502)
    p = (1.0, 1.5, 3.0, 3.5, 4.0, 7.5)
    s = sum(a[i] * tau ** p[i] for i in range(6))
    return PC_H2O * np.exp((TC_H2O / T_K) * s)


# =====================================================================
# 2. PURE WATER & BRINE DENSITY (Batzle & Wang 1992)
# =====================================================================
def density_water_BW(T_K, P_Pa):
    """Pure water density (kg/m^3) via Batzle & Wang (1992).
    Valid ~5 - 350 degC, 0.1 - 100 MPa; accuracy ~0.5 %."""
    T = T_K - 273.15            # degC
    P = P_Pa * 1.0e-6           # MPa
    rho_w = 1.0 + 1.0e-6 * (
        -80.0 * T - 3.3 * T ** 2 + 0.00175 * T ** 3
        + 489.0 * P - 2.0 * T * P + 0.016 * T ** 2 * P
        - 1.3e-5 * T ** 3 * P - 0.333 * P ** 2 - 0.002 * T * P ** 2
    )
    return rho_w * 1000.0       # g/cm^3 -> kg/m^3


def density_brine_BW(T_K, P_Pa, m_NaCl):
    """Pure NaCl brine density (kg/m^3) via Batzle & Wang (1992).
    Salinity m_NaCl in mol/kg of H2O.  Does NOT include the effect of
    dissolved CO2 - use `density_live_brine` for that.
    Valid ~5 - 350 degC, 0.1 - 100 MPa, 0 - 320 000 ppm NaCl."""
    T = T_K - 273.15
    P = P_Pa * 1.0e-6
    S = molality_to_mass_fraction(m_NaCl)               # mass fraction
    rho_w = density_water_BW(T_K, P_Pa) / 1000.0        # back to g/cm^3
    rho_b = rho_w + S * (
        0.668 + 0.44 * S
        + 1.0e-6 * (
            300.0 * P - 2400.0 * P * S
            + T * (80.0 + 3.0 * T - 3300.0 * S - 13.0 * P + 47.0 * P * S)
        )
    )
    return rho_b * 1000.0


# =====================================================================
# 3. WATER & BRINE VISCOSITY  (Vogel + Mao-Duan 2009 ratio)
# =====================================================================
def viscosity_water(T_K):
    """Pure water dynamic viscosity (Pa.s).
    Vogel-type fit; accuracy ~1 % from 273 - 573 K.
    Pressure dependence is weak below ~50 MPa and is neglected."""
    return 2.414e-5 * 10.0 ** (247.8 / (T_K - 140.0))


def viscosity_brine(T_K, P_Pa, m_NaCl):
    """NaCl brine dynamic viscosity (Pa.s) via Mao & Duan (2009).
    Pressure dependence below ~100 MPa is weak (<5 %) and is omitted
    here for simplicity, matching the original Mao-Duan ratio form."""
    T = T_K
    a = -0.21319213 + 0.13651589e-2 * T - 0.12191756e-5 * T ** 2
    b =  0.69161945e-1 - 0.27292263e-3 * T + 0.20852448e-6 * T ** 2
    c = -0.25988855e-2 + 0.77989227e-5 * T
    m = m_NaCl
    ln_ratio = a * m + b * m ** 2 + c * m ** 3
    return viscosity_water(T_K) * np.exp(ln_ratio)


# =====================================================================
# 4. PURE CO2 - PENG-ROBINSON EOS  (density, fugacity coefficient)
# =====================================================================
def _PR_AB_CO2(T_K, P_Pa):
    """Return (a, b, A, B) for pure CO2 in PR-EOS at (T, P)."""
    Tr = T_K / TC_CO2
    kappa = 0.37464 + 1.54226 * OMEGA_CO2 - 0.26992 * OMEGA_CO2 ** 2
    alpha = (1.0 + kappa * (1.0 - np.sqrt(Tr))) ** 2
    a = 0.45724 * R_GAS ** 2 * TC_CO2 ** 2 / PC_CO2 * alpha
    b = 0.07780 * R_GAS * TC_CO2 / PC_CO2
    A = a * P_Pa / (R_GAS * T_K) ** 2
    B = b * P_Pa / (R_GAS * T_K)
    return a, b, A, B


def Z_factor_CO2(T_K, P_Pa, phase='auto'):
    """Compressibility factor Z of pure CO2 from PR-EOS.

    phase :
        'gas'    - largest real root above B
        'liquid' - smallest real root above B
        'auto'   - root with lower Gibbs energy when 2 real roots exist
    """
    if not (np.isfinite(T_K) and np.isfinite(P_Pa)) or P_Pa <= 0.0 or T_K <= 0.0:
        raise ValueError(
            f"PR-EOS: invalid inputs T={T_K} K, P={P_Pa} Pa")
    _a, _b, A, B = _PR_AB_CO2(T_K, P_Pa)
    if not (np.isfinite(A) and np.isfinite(B)):
        raise ValueError(
            f"PR-EOS: non-finite A,B at T={T_K} K, P={P_Pa} Pa")
    # Z^3 - (1-B) Z^2 + (A - 3 B^2 - 2 B) Z - (A B - B^2 - B^3) = 0
    coeffs = [1.0,
              -(1.0 - B),
              (A - 3.0 * B ** 2 - 2.0 * B),
              -(A * B - B ** 2 - B ** 3)]
    if not all(np.isfinite(c) for c in coeffs):
        raise ValueError(
            f"PR-EOS: non-finite cubic coeffs at T={T_K} K, P={P_Pa} Pa")
    roots = np.roots(coeffs)
    real_roots = sorted(r.real for r in roots
                        if abs(r.imag) < 1.0e-9 and r.real > B)
    if not real_roots:
        # Fall back: take largest real root regardless of sign check
        real_roots = sorted(r.real for r in roots if abs(r.imag) < 1.0e-9)
        if not real_roots:
            raise RuntimeError(
                f"PR-EOS: no real root for CO2 at T={T_K} K, P={P_Pa} Pa")
        return real_roots[-1]
    if len(real_roots) == 1:
        return real_roots[0]
    if phase == 'gas':
        return real_roots[-1]
    if phase == 'liquid':
        return real_roots[0]
    # auto - pick min Gibbs energy
    Z_l, Z_v = real_roots[0], real_roots[-1]
    sq2 = np.sqrt(2.0)
    # Difference in dimensionless Gibbs energy (Michelsen)
    dG = ((Z_v - Z_l) - np.log((Z_v - B) / (Z_l - B))
          - A / (2.0 * sq2 * B) * np.log(
              ((Z_v + (1 + sq2) * B) * (Z_l + (1 - sq2) * B))
              / ((Z_v + (1 - sq2) * B) * (Z_l + (1 + sq2) * B))
          ))
    return Z_v if dG < 0 else Z_l


def density_CO2(T_K, P_Pa, phase='auto'):
    """CO2 mass density (kg/m^3) via Peng-Robinson EOS.

    Note
    ----
    PR-EOS is known to under-predict dense supercritical CO2 density by
    typically 5 - 12 %.  For higher fidelity use Span-Wagner (1996),
    not implemented here.
    """
    Z = Z_factor_CO2(T_K, P_Pa, phase=phase)
    Vm = Z * R_GAS * T_K / P_Pa            # m^3 / mol
    return M_CO2 / Vm                       # kg / m^3


def fugacity_coefficient_CO2(T_K, P_Pa, phase='auto'):
    """Fugacity coefficient phi of pure CO2 from PR-EOS."""
    _a, _b, A, B = _PR_AB_CO2(T_K, P_Pa)
    Z = Z_factor_CO2(T_K, P_Pa, phase=phase)
    sq2 = np.sqrt(2.0)
    ln_phi = ((Z - 1.0) - np.log(Z - B)
              - A / (2.0 * sq2 * B) * np.log(
                  (Z + (1 + sq2) * B) / (Z + (1 - sq2) * B)
              ))
    return np.exp(ln_phi)


# =====================================================================
# 5. CO2 VISCOSITY  (Fenghour, Wakeham, Vesovic 1998)
# =====================================================================
def viscosity_CO2(T_K, rho_kgm3):
    """CO2 dynamic viscosity (Pa.s) via Fenghour et al. (1998).
    Critical enhancement term is omitted (small except very near
    the critical point)."""
    # Zero-density viscosity (uPa.s)
    a_coef = (0.235156, -0.491266, 5.211155e-2,
              5.347906e-2, -1.537102e-2)
    Tstar = T_K / 251.196
    lnTs = np.log(Tstar)
    lnGstar = sum(a_coef[i] * lnTs ** i for i in range(5))
    Gstar = np.exp(lnGstar)
    eta0 = 1.00697 * np.sqrt(T_K) / Gstar             # uPa.s

    # Excess viscosity (uPa.s) - Fenghour et al. d-coefficients
    d11, d21 = 0.4071119e-2, 0.7198037e-4
    d64       = 0.2411697e-16
    d81, d82 = 0.2971072e-22, -0.1627888e-22
    rho = rho_kgm3
    eta_ex = (d11 * rho
              + d21 * rho ** 2
              + d64 * rho ** 6 / Tstar ** 3
              + d81 * rho ** 8
              + d82 * rho ** 8 / Tstar)

    return (eta0 + eta_ex) * 1.0e-6                   # uPa.s -> Pa.s


# =====================================================================
# 6. CO2 SOLUBILITY IN BRINE  (Duan & Sun 2003)
# =====================================================================
def _DS_mu_CO2_RT(T_K, P_bar):
    """Duan-Sun reference chemical potential mu_CO2^l(0)/RT  (T in K,
    P in bar). Coefficients from Duan & Sun (2003) Table 2."""
    T, P = T_K, P_bar
    return (28.9447706
            - 0.0354581768 * T
            - 4770.67077 / T
            + 1.02782768e-5 * T ** 2
            + 33.8126098 / (630.0 - T)
            + 9.04037140e-3 * P
            - 1.14934031e-3 * P * np.log(T)
            - 0.307405726 * P / T
            - 0.0907301486 * P / (630.0 - T)
            + 9.32713393e-4 * P ** 2 / (630.0 - T) ** 2)


def _DS_lambda_CO2_Na(T_K, P_bar):
    """Duan-Sun lambda_{CO2-Na} interaction parameter."""
    T, P = T_K, P_bar
    return (-0.411370585
            + 6.07632013e-4 * T
            + 97.5347708 / T
            - 0.0237622469 * P / T
            + 0.0170656236 * P / (630.0 - T)
            + 1.41335834e-5 * T * np.log(P))


def _DS_zeta_CO2_NaCl(T_K, P_bar):
    """Duan-Sun zeta_{CO2-Na-Cl} ternary interaction parameter."""
    T, P = T_K, P_bar
    return (3.36389723e-4
            - 1.98298980e-5 * T
            + 2.12220830e-3 * P / T
            - 5.24873303e-3 * P / (630.0 - T))


def CO2_solubility_in_brine(T_K, P_Pa, m_NaCl):
    """Maximum CO2 molality in NaCl brine at saturation, via Duan & Sun
    (2003).  Returns mol CO2 per kg of H2O.

    Valid range: 273 - 533 K, 0 - 2000 bar, 0 - 4.3 mol/kg NaCl.
    """
    P_bar = P_Pa * 1.0e-5
    if P_bar < 1.0e-3:        # avoid log(0) in lambda term
        P_bar = 1.0e-3

    # Vapor-phase mole fraction of CO2 (water vapor in CO2-rich phase
    # treated by Raoult-like correction; valid for P >> P_sat,H2O)
    P_sat_w = saturation_pressure_water(T_K)
    y_H2O = min(0.999, P_sat_w / P_Pa)
    y_CO2 = 1.0 - y_H2O

    phi = fugacity_coefficient_CO2(T_K, P_Pa, phase='auto')
    mu0 = _DS_mu_CO2_RT(T_K, P_bar)
    lam = _DS_lambda_CO2_Na(T_K, P_bar)
    zet = _DS_zeta_CO2_NaCl(T_K, P_bar)

    ln_m = (np.log(y_CO2 * phi * P_bar)
            - mu0
            - 2.0 * lam * m_NaCl
            - zet * m_NaCl ** 2)
    return np.exp(ln_m)


# =====================================================================
# 7. LIVE BRINE DENSITY  (Garcia 2001)
# =====================================================================
def density_live_brine(T_K, P_Pa, m_NaCl, m_CO2_dissolved):
    """Density of brine + dissolved CO2 (kg/m^3).

    Uses the Garcia (2001) apparent molar volume of CO2 in brine to
    correct the density of the CO2-free brine.
    """
    rho_b0 = density_brine_BW(T_K, P_Pa, m_NaCl)         # CO2-free brine

    # Garcia (2001) apparent molar volume of CO2 in brine, cm^3/mol
    T = T_K - 273.15
    V_phi_cm3 = 37.51 - 9.585e-2 * T + 8.74e-4 * T ** 2 - 5.044e-7 * T ** 3
    V_phi = V_phi_cm3 * 1.0e-6                           # m^3/mol

    # Per 1 kg of H2O:
    #   mass of CO2-free brine  = 1 + m_NaCl * M_NaCl   kg
    #   mass of dissolved CO2   = m_CO2_dissolved * M_CO2 kg
    m_brine0 = 1.0 + m_NaCl * M_NACL
    m_CO2 = m_CO2_dissolved * M_CO2
    m_total = m_brine0 + m_CO2

    V_brine0 = m_brine0 / rho_b0                         # m^3
    V_CO2 = m_CO2_dissolved * V_phi                       # m^3
    V_total = V_brine0 + V_CO2

    return m_total / V_total


# =====================================================================
# 8. INTERFACIAL TENSION (engineering correlation)
# =====================================================================
def surface_tension_water(T_K):
    """Air-water surface tension (N/m) - IAPWS R1-76(2014)."""
    if T_K >= TC_H2O:
        return 0.0
    tau = 1.0 - T_K / TC_H2O
    return 235.8e-3 * tau ** 1.256 * (1.0 - 0.625 * tau)


def interfacial_tension_brine_CO2(T_K, P_Pa, m_NaCl):
    """Brine-CO2 interfacial tension (N/m).

    Engineering correlation: IAPWS pure-water surface tension reduced
    by an empirical pressure trend fit to Bachu & Bennion (2009) brine-
    CO2 data, with a small linear salinity correction in the spirit
    of Chalbaud et al. (2009).  Approximate accuracy ~15 % - adequate
    for two-phase flow correlations where surface tension enters with
    a weak (typically 0.25) exponent.
    """
    sigma_w = surface_tension_water(T_K) * 1000.0        # mN/m
    P_bar = max(0.0, P_Pa * 1.0e-5)                       # clamp non-negative
    # Asymptote toward ~25 mN/m at high P (CO2 saturation)
    sigma = sigma_w - 47.0 * (1.0 - np.exp(-P_bar / 80.0))
    sigma = max(20.0, sigma)
    sigma *= (1.0 + 0.02 * m_NaCl)                        # salinity bump
    return sigma * 1.0e-3                                  # mN/m -> N/m


# =====================================================================
# 9. EQUILIBRIUM FLASH OF CO2 + BRINE  (the workhorse for the VLP)
# =====================================================================
def flash_CO2_brine(T_K, P_Pa, m_NaCl, w_CO2_total):
    """Equilibrium flash of CO2 + NaCl brine at local (T, P).

    Parameters
    ----------
    T_K           : temperature (K)
    P_Pa          : pressure (Pa)
    m_NaCl        : NaCl molality of the brine (mol / kg H2O)
    w_CO2_total   : total CO2 mass fraction in the produced mixture
                    (kg CO2 per kg of total mixture)

    Returns
    -------
    dict with keys:
        x_gas           - mass fraction of free-gas (CO2) phase
        x_liq           - mass fraction of liquid (live brine) phase
        rho_gas         - free CO2 density (kg/m^3); None if no gas
        rho_liq         - live brine density (kg/m^3)
        mu_gas          - CO2 viscosity (Pa.s); None if no gas
        mu_liq          - live brine viscosity (Pa.s)
        sigma_gl        - gas-liquid interfacial tension (N/m)
        m_CO2_dis       - dissolved CO2 molality (mol/kg H2O)
        m_CO2_total_aq  - molality if all CO2 were dissolved (mol/kg H2O)
        m_CO2_max       - solubility limit at (T, P) (mol/kg H2O)
        is_two_phase    - True iff free-gas phase is present
    """
    if not (0.0 <= w_CO2_total < 1.0):
        raise ValueError("w_CO2_total must lie in [0, 1).")

    # Mass per 1 kg of total (CO2 + brine) mixture, before any flash:
    #   (1 - w_CO2_total) kg of CO2-free brine, which contains:
    #   m_water = brine_mass / (1 + m_NaCl * M_NaCl)   kg of H2O
    m_brine_in = 1.0 - w_CO2_total
    m_water = m_brine_in / (1.0 + m_NaCl * M_NACL)

    # If all CO2 were dissolved -> molality:
    n_CO2_total = w_CO2_total / M_CO2
    m_CO2_total_aq = n_CO2_total / m_water if m_water > 0 else np.inf

    # Solubility limit at local conditions:
    m_CO2_max = CO2_solubility_in_brine(T_K, P_Pa, m_NaCl)

    if m_CO2_total_aq <= m_CO2_max:
        # ---- Single phase: under-saturated, all CO2 dissolved ----
        m_CO2_dis = m_CO2_total_aq
        rho_liq = density_live_brine(T_K, P_Pa, m_NaCl, m_CO2_dis)
        mu_liq  = viscosity_brine(T_K, P_Pa, m_NaCl)
        sigma   = interfacial_tension_brine_CO2(T_K, P_Pa, m_NaCl)
        return dict(
            x_gas=0.0, x_liq=1.0,
            rho_gas=None, rho_liq=rho_liq,
            mu_gas=None,  mu_liq=mu_liq,
            sigma_gl=sigma,
            m_CO2_dis=m_CO2_dis,
            m_CO2_total_aq=m_CO2_total_aq,
            m_CO2_max=m_CO2_max,
            is_two_phase=False,
        )

    # ---- Two phase: liquid is at solubility limit, rest is free gas ----
    m_CO2_dis = m_CO2_max
    # Free CO2 mass per kg of total mixture
    m_CO2_free = (m_CO2_total_aq - m_CO2_max) * M_CO2 * m_water
    x_gas = m_CO2_free
    x_liq = 1.0 - x_gas

    rho_liq = density_live_brine(T_K, P_Pa, m_NaCl, m_CO2_dis)
    mu_liq  = viscosity_brine(T_K, P_Pa, m_NaCl)
    rho_gas = density_CO2(T_K, P_Pa, phase='gas')
    mu_gas  = viscosity_CO2(T_K, rho_gas)
    sigma   = interfacial_tension_brine_CO2(T_K, P_Pa, m_NaCl)

    return dict(
        x_gas=x_gas, x_liq=x_liq,
        rho_gas=rho_gas, rho_liq=rho_liq,
        mu_gas=mu_gas,   mu_liq=mu_liq,
        sigma_gl=sigma,
        m_CO2_dis=m_CO2_dis,
        m_CO2_total_aq=m_CO2_total_aq,
        m_CO2_max=m_CO2_max,
        is_two_phase=True,
    )


# =====================================================================
# 10. STANDARD-CONDITION GWR CONVERSIONS
# =====================================================================
def GWR_std_to_mass_fraction(GWR_std, m_NaCl=0.0,
                             T_std=T_STD, P_std=P_STD):
    """Convert standard-condition gas-water ratio (m^3 CO2 at standard
    conditions per m^3 of brine at standard conditions) to total CO2
    mass fraction in the produced stream."""
    rho_CO2_std   = density_CO2(T_std, P_std, phase='gas')
    rho_brine_std = density_brine_BW(T_std, P_std, m_NaCl)
    mass_CO2   = GWR_std * rho_CO2_std       # kg CO2 per m^3 brine
    mass_brine = rho_brine_std               # kg brine per m^3 brine
    return mass_CO2 / (mass_CO2 + mass_brine)


def mass_fraction_to_GWR_std(w_CO2, m_NaCl=0.0,
                             T_std=T_STD, P_std=P_STD):
    """Inverse of GWR_std_to_mass_fraction."""
    rho_CO2_std   = density_CO2(T_std, P_std, phase='gas')
    rho_brine_std = density_brine_BW(T_std, P_std, m_NaCl)
    if w_CO2 >= 1.0:
        return np.inf
    mass_CO2_per_kg_total   = w_CO2
    mass_brine_per_kg_total = 1.0 - w_CO2
    V_CO2   = mass_CO2_per_kg_total / rho_CO2_std
    V_brine = mass_brine_per_kg_total / rho_brine_std
    return V_CO2 / V_brine


# =====================================================================
# 11. SELF-TEST / VALIDATION BLOCK
# =====================================================================
def _print_row(label, computed, ref, units, tol_pct=10.0):
    err = abs(computed - ref) / abs(ref) * 100.0 if ref != 0 else 0.0
    flag = "OK" if err <= tol_pct else "CHECK"
    print(f"  {label:<48s}  {computed:12.4g}   "
          f"(ref {ref:10.4g}) {units:<10s}  err {err:5.1f}%  [{flag}]")


def _validate():
    print("=" * 78)
    print("PVT MODULE - VALIDATION AGAINST PUBLISHED EXPERIMENTAL DATA")
    print("=" * 78)
    print("Reference sources: NIST Webbook (CO2, water), IAPWS, Wiebe-Gaddy")
    print("(1939,1940), Bamberger et al. (2000), Tödheide-Franck (1963),")
    print("Kestin et al. (1981).")
    print("=" * 78)

    # ----- Pure water -----
    print("\n[Pure water]")
    rho = density_water_BW(C_to_K(25.0), bar_to_Pa(1.0))
    _print_row("rho_water(25 degC, 1 bar)", rho, 997.05, "kg/m^3", 1.0)
    rho = density_water_BW(C_to_K(150.0), bar_to_Pa(50.0))
    _print_row("rho_water(150 degC, 50 bar)", rho, 919.7, "kg/m^3", 2.0)
    mu = viscosity_water(C_to_K(25.0))
    _print_row("mu_water(25 degC)", mu * 1e3, 0.890, "mPa.s", 3.0)
    mu = viscosity_water(C_to_K(100.0))
    _print_row("mu_water(100 degC)", mu * 1e3, 0.282, "mPa.s", 5.0)
    psat = saturation_pressure_water(C_to_K(100.0))
    _print_row("Psat_water(100 degC)", psat * 1e-5, 1.0142, "bar", 1.0)
    psat = saturation_pressure_water(C_to_K(200.0))
    _print_row("Psat_water(200 degC)", psat * 1e-5, 15.55, "bar", 1.0)

    # ----- Brine -----
    print("\n[NaCl brine - Batzle-Wang density, Mao-Duan viscosity]")
    m1 = ppm_to_molality(50_000)        # 50 000 ppm = 0.879 mol/kg
    rho = density_brine_BW(C_to_K(25.0), bar_to_Pa(1.0), m1)
    _print_row("rho_brine(25 degC, 1 bar, 50,000 ppm)",
               rho, 1033.5, "kg/m^3", 2.0)
    mu = viscosity_brine(C_to_K(25.0), bar_to_Pa(1.0), m1)
    _print_row("mu_brine(25 degC, 1 bar, 50,000 ppm)",
               mu * 1e3, 0.97, "mPa.s", 5.0)
    m2 = ppm_to_molality(100_000)       # 100 000 ppm = 1.90 mol/kg
    rho = density_brine_BW(C_to_K(150.0), bar_to_Pa(100.0), m2)
    _print_row("rho_brine(150 degC, 100 bar, 100,000 ppm)",
               rho, 998.0, "kg/m^3", 3.0)

    # ----- CO2 density (Peng-Robinson EOS, NIST as reference) -----
    print("\n[Pure CO2 - PR-EOS density (NIST reference values)]")
    rho = density_CO2(C_to_K(25.0), bar_to_Pa(1.0))
    _print_row("rho_CO2(25 degC, 1 bar)", rho, 1.808, "kg/m^3", 5.0)
    rho = density_CO2(C_to_K(50.0), bar_to_Pa(50.0))
    _print_row("rho_CO2(50 degC, 50 bar) [dense gas]",
               rho, 119.3, "kg/m^3", 15.0)
    rho = density_CO2(C_to_K(100.0), bar_to_Pa(200.0))
    _print_row("rho_CO2(100 degC, 200 bar) [supercritical]",
               rho, 481.0, "kg/m^3", 10.0)
    mu = viscosity_CO2(C_to_K(25.0), 1.808)
    _print_row("mu_CO2(25 degC, ~1 bar)", mu * 1e6, 14.93, "uPa.s", 5.0)
    mu = viscosity_CO2(C_to_K(100.0), density_CO2(C_to_K(100.0),
                                                  bar_to_Pa(200.0)))
    _print_row("mu_CO2(100 degC, 200 bar)", mu * 1e6, 35.3, "uPa.s", 20.0)

    # ----- Duan-Sun CO2 solubility (vs. experimental data) -----
    print("\n[CO2 solubility - Duan & Sun (2003) with PR fugacity]")
    print("Reference values are EXPERIMENTAL DATA - tolerances reflect"
          " PR-EOS\nlimitations vs. the original Duan-EOS used by"
          " Duan & Sun.")
    s = CO2_solubility_in_brine(C_to_K(25.0), bar_to_Pa(50.0), 0.0)
    _print_row("m_CO2(25 degC, 50 bar, pure H2O) Bamberger",
               s, 1.054, "mol/kg", 15.0)
    s = CO2_solubility_in_brine(C_to_K(60.0), bar_to_Pa(100.0), 0.0)
    _print_row("m_CO2(60 degC, 100 bar, pure H2O) Wiebe-Gaddy",
               s, 1.000, "mol/kg", 10.0)
    s = CO2_solubility_in_brine(C_to_K(100.0), bar_to_Pa(100.0), 0.0)
    _print_row("m_CO2(100 degC, 100 bar, pure H2O) Wiebe-Gaddy",
               s, 0.875, "mol/kg", 15.0)
    s = CO2_solubility_in_brine(C_to_K(100.0), bar_to_Pa(200.0), 0.0)
    _print_row("m_CO2(100 degC, 200 bar, pure H2O) Wiebe-Gaddy",
               s, 1.292, "mol/kg", 20.0)
    s = CO2_solubility_in_brine(C_to_K(150.0), bar_to_Pa(200.0), 0.0)
    _print_row("m_CO2(150 degC, 200 bar, pure H2O) Toedheide",
               s, 1.135, "mol/kg", 10.0)
    # Salting-out: ratio at 1 mol/kg NaCl vs pure water should be ~0.85
    s_pw = CO2_solubility_in_brine(C_to_K(60.0), bar_to_Pa(100.0), 0.0)
    s_b1 = CO2_solubility_in_brine(C_to_K(60.0), bar_to_Pa(100.0), 1.0)
    print(f"  Salting-out @60C,100bar: m_CO2(1 m NaCl)/m_CO2(0) = "
          f"{s_b1/s_pw:.3f}  (Drummond exp. ~0.83)")

    # ----- Live brine density -----
    print("\n[Live brine (CO2 dissolved, Garcia 2001)]")
    msat = CO2_solubility_in_brine(C_to_K(60.0), bar_to_Pa(100.0), 1.0)
    rho0 = density_brine_BW(C_to_K(60.0), bar_to_Pa(100.0), 1.0)
    rho1 = density_live_brine(C_to_K(60.0), bar_to_Pa(100.0), 1.0, msat)
    print(f"  rho_brine(60 degC, 100 bar, 1m NaCl)        = "
          f"{rho0:8.2f} kg/m^3  (CO2-free)")
    print(f"  rho_live_brine (CO2 sat = {msat:.3f} mol/kg) = "
          f"{rho1:8.2f} kg/m^3")
    print(f"  Density increase from dissolved CO2:         "
          f"{rho1 - rho0:+.2f} kg/m^3  (~+5-10 kg/m^3 expected)")

    # ----- Surface tension -----
    print("\n[Brine-CO2 interfacial tension]")
    s = interfacial_tension_brine_CO2(C_to_K(50.0), bar_to_Pa(100.0), 0.0)
    print(f"  sigma(50 degC, 100 bar, pure H2O):       "
          f"{s * 1e3:5.1f} mN/m  (Bachu-Bennion ~30)")
    s = interfacial_tension_brine_CO2(C_to_K(100.0), bar_to_Pa(200.0), 1.0)
    print(f"  sigma(100 degC, 200 bar, 1 mol/kg NaCl): "
          f"{s * 1e3:5.1f} mN/m")

    # ----- Flash demonstration: profile down a hypothetical well -----
    print("\n[Flash demo - typical low-GWR geothermal well profile]")
    GWR = 10.0
    m_NaCl = ppm_to_molality(50_000)
    w_CO2 = GWR_std_to_mass_fraction(GWR, m_NaCl)
    print(f"  GWR_std = {GWR:.1f} Sm^3/Sm^3, NaCl 50,000 ppm")
    print(f"  => total CO2 mass fraction in produced stream: "
          f"{w_CO2 * 100:.3f} wt%\n")
    print(f"  {'T degC':>7} {'P bar':>6} | {'phase':<10} "
          f"{'x_gas wt%':>9} {'rho_L':>7} {'rho_G':>7} "
          f"{'mu_L mPas':>10} {'m_CO2_dis':>10}")
    print(f"  {'-'*7} {'-'*6} | {'-'*10} {'-'*9} {'-'*7} {'-'*7} "
          f"{'-'*10} {'-'*10}")
    for T_C, P_bar in [(200, 250), (180, 200), (170, 150), (160, 100),
                       (150,  60), (140,  30), (130,  10)]:
        f = flash_CO2_brine(C_to_K(T_C), bar_to_Pa(P_bar),
                            m_NaCl, w_CO2)
        phase = "2-phase" if f['is_two_phase'] else "1-phase L"
        rg = f"{f['rho_gas']:7.1f}" if f['rho_gas'] else "    -- "
        print(f"  {T_C:7d} {P_bar:6d} | {phase:<10} "
              f"{f['x_gas']*100:9.3f} {f['rho_liq']:7.1f} {rg} "
              f"{f['mu_liq']*1e3:10.4f} {f['m_CO2_dis']:10.3f}")

    print("\n" + "=" * 78)
    print("MODEL NOTES:")
    print("  * Density of pure water and brine: <1 % vs. NIST/IAPWS data.")
    print("  * Brine viscosity: ~3-5 % vs. Kestin et al. tabulations.")
    print("  * CO2 density (PR-EOS): ~5-15 % low in dense supercritical")
    print("    region. Acceptable for nodal analysis; upgrade to Span-")
    print("    Wagner only if needed.")
    print("  * CO2 solubility in brine: ~5-15 % vs. experiment. Larger")
    print("    errors (up to ~20 %) at intermediate T,P (~80-100 degC,")
    print("    100-200 bar) due to PR fugacity vs. Duan-EOS fugacity.")
    print("    Salting-out trend (Sechenov coefficient) is well captured.")
    print("  * Interfacial tension: engineering estimate (~15 %); enters")
    print("    H&B with weak (0.25) exponent so error is amplified <4 %.")
    print("=" * 78)


if __name__ == "__main__":
    _validate()
