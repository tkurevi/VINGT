# NOTE: prices/OPEX/CO2 in this file are NOT used by the economic model.
# ECONOMICS HANDLED IN economy*.py (this script supplies physical sizing only).
"""
=============================================================================
ESP (Electrical Submersible Pump) Power Calculator for Geothermal Brine Wells
Production Well  +  Injection Well  —  Full System Overview
=============================================================================
Assumes single-phase liquid water only (no flash / no gas phase).
All fluid properties are computed from temperature via established correlations.

PRODUCTION WELL:  ESP lifts brine from reservoir to surface.
                  TDH = h_static(lift) + h_friction + h_minor + h_backpressure

INJECTION WELL:   ESP (downhole or surface-equivalent) pushes cooled brine
                  into the formation against reservoir pressure.
                  TDH = h_reservoir - h_hydrostatic(gravity assist)
                        + h_friction + h_minor
                  If TDH <= 0 no pump is needed (gravity injection).

References:
  - IAPWS-IF97 simplified correlations for water properties
  - Colebrook-White equation for Darcy friction factor
  - Darcy-Weisbach for pipe friction losses
  - API / pump engineering standards for ESP sizing
=============================================================================
"""

import math

# ===========================================================================
# USER INPUTS  <- edit these values
# ===========================================================================

# ---------------------------------------------------------------------------
# PRODUCTION WELL
# ---------------------------------------------------------------------------
PROD_DEPTH_M          = 390  # [m]    TDH Total Dynamic Head (ESP depth/pressure minus intake pressure/depth)
PROD_FLOW_LS          = 34.7    # [L/s]  Required production flow rate
PROD_TEMP_C           = 94   # [deg C] Brine temperature at pump depth

# Production well efficiencies
# Typical ranges:
#   pump_eff  : 0.55-0.75  (centrifugal ESP stages; 0.65 good mid-range)
#   motor_eff : 0.85-0.95  (ESP motor; 0.90 industry standard)
#   cable_eff : 0.95-0.99  (power cable; 0.97 typical for deep wells)
#   vsd_eff   : 0.95-0.98  (Variable Speed Drive; set 1.0 if no VSD)
PROD_PUMP_EFF         = 0.70    # [-]  pump hydraulic efficiency
PROD_MOTOR_EFF        = 0.90    # [-]  motor efficiency
PROD_CABLE_EFF        = 0.97    # [-]  power cable efficiency
PROD_VSD_EFF          = 0.97    # [-]  variable speed drive efficiency

# Production well geometry
PROD_PIPE_ID_M        = 0.1016  # [m]   Tubing inner diameter (4" nominal = 0.1016 m)
PROD_PIPE_ROUGHNESS_M = 4.5e-5  # [m]   Absolute roughness (commercial steel = 4.5e-5 m)
PROD_BACKPRESSURE_BAR = 3.2     # [bar] Wellhead back-pressure above atmospheric
                                 #       e.g. 3.0 = closed district heating network

# ---------------------------------------------------------------------------
# INJECTION WELL
# ---------------------------------------------------------------------------
INJ_DEPTH_M           = 2700 # [m]    Depth to injection perforations / pump setting
INJ_FLOW_LS           = 35.3    # [L/s]  Injection flow rate (often = production flow)
INJ_TEMP_C            = 60    # [deg C] Temperature of fluid at injection pump
                                 #         = reinjection temperature after heat extraction.
                                 #         Used for both injection well hydraulics and
                                 #         geothermal thermal power calculation.
INJ_RESERVOIR_BAR     = 310    # [bar]  Required injection pressure at perforations
                                 #        (gauge). Pressure needed to push fluid into
                                 #        the formation. Typical: 5-30 bar gauge
                                 #        depending on reservoir permeability & depth.
                                 #        Hydrostatic column of fluid ASSISTS injection,
                                 #        so TDH = P_reservoir/rho/g - depth + friction.
                                 #        If TDH <= 0 gravity injection, no pump needed.
INJ_WELLHEAD_BAR      = 2.5     # [bar]  Additional surface back-pressure at the
                                 #        injection wellhead (gauge). Typically 0
                                 #        (open tank / direct pipe connection).
                                 #        Set > 0 if a surface control valve or
                                 #        pressurised manifold adds resistance on the
                                 #        discharge side of the injection pump.

# Injection well efficiencies
# Injection pumps are often surface centrifugal pumps; downhole ESP same ranges apply.
#   pump_eff  : 0.60-0.75  (injection pumps often higher eff. at lower head)
#   motor_eff : 0.85-0.95
#   cable_eff : 0.95-0.99  (shorter cable if surface pump -> higher efficiency)
#   vsd_eff   : 0.95-0.98  (set 1.0 if no VSD)
INJ_PUMP_EFF          = 0.70    # [-]  pump hydraulic efficiency
INJ_MOTOR_EFF         = 0.90    # [-]  motor efficiency
INJ_CABLE_EFF         = 0.98    # [-]  power cable efficiency (shorter cable for inj.)
INJ_VSD_EFF           = 0.97    # [-]  variable speed drive efficiency

# Injection well geometry
INJ_PIPE_ID_M         = 0.20  # [m]   Injection tubing inner diameter
INJ_PIPE_ROUGHNESS_M  = 4.5e-5  # [m]   Absolute roughness

# ---------------------------------------------------------------------------
# SHARED OPERATIONAL & ECONOMIC INPUTS
# ---------------------------------------------------------------------------
HOURS_PER_YEAR        = 4000    # [h/yr]      Operating hours per year (8760 = continuous)

PRICE_HEAT_EUR_MWhth  = 43    # [EUR/MWhth] Selling price of thermal energy
                                 #   Typical district heating range: 30-80 EUR/MWhth
PRICE_ELEC_EUR_MWhe   = 135.0    # [EUR/MWhe]  Electricity purchase price (ESP OPEX)
                                 #   Typical industrial tariff: 60-150 EUR/MWhe

# ---------------------------------------------------------------------------
# CO2 SAVINGS vs. NATURAL GAS REFERENCE PLANT
# ---------------------------------------------------------------------------
CO2_PRICE_EUR_TON     = 65    # [EUR/tCO2] CO2 market price (EU ETS or voluntary)
                                 #   EU ETS historical range: 25-100 EUR/tCO2
                                 #   Typical project assumption: 50-80 EUR/tCO2

GAS_BOILER_EFF        = 0.85    # [-]  Thermal efficiency of reference natural gas
                                 #      boiler / plant producing the same heat output
                                 #   Typical values:
                                 #     0.82-0.88 : conventional boiler (non-condensing)
                                 #     0.88-0.95 : condensing gas boiler
                                 #     0.85-0.92 : gas-fired CHP heat output side
                                 #
                                 # Emission factor used: 0.202 tCO2/MWhgas (LHV basis)
                                 # Source: IPCC 2006 Guidelines, natural gas = 56.1 kgCO2/GJ


# ===========================================================================
# WATER PROPERTY CORRELATIONS  (valid ~0-200 deg C, atmospheric to moderate P)
# ===========================================================================

def water_density(T_C):
    """
    Water density [kg/m3] vs temperature [deg C].
    5th-order polynomial fit to IAPWS-IF97 data (+/-0.3% over 0-200 deg C).
    """
    T = T_C
    return (999.842594
            + 6.793952e-2 * T
            - 9.095290e-3 * T**2
            + 1.001685e-4 * T**3
            - 1.120083e-6 * T**4
            + 6.536332e-9 * T**5)


def water_dynamic_viscosity(T_C):
    """
    Dynamic viscosity [Pa.s] vs temperature [deg C].
    Korosi & Fabuss (1968) Vogel-type equation (+/-2% over 0-200 deg C).
    """
    T_K = T_C + 273.15
    return 2.414e-5 * 10 ** (247.8 / (T_K - 140.0))


def water_vapor_pressure(T_C):
    """
    Saturation vapour pressure [Pa] - Buck equation.
    Used for NPSHa / cavitation check.
    """
    return 611.21 * math.exp((18.678 - T_C / 234.5) * (T_C / (257.14 + T_C)))


def water_specific_heat(T_C):
    """Isobaric specific heat [J/(kg.K)] - IAPWS polynomial fit."""
    T = T_C
    return (4217.6
            - 3.820    * T
            + 0.1103   * T**2
            - 9.609e-4 * T**3
            + 3.166e-6 * T**4)


# ===========================================================================
# HYDRAULIC HELPERS
# ===========================================================================

def friction_factor_colebrook(Re, relative_roughness):
    """
    Darcy friction factor via Colebrook-White (iterative, Haaland initial guess).
    Returns 64/Re for laminar flow (Re < 2300).
    """
    if Re < 2300:
        return 64.0 / Re
    eps_D = relative_roughness
    f = (-1.8 * math.log10((eps_D / 3.7) ** 1.11 + 6.9 / Re)) ** -2
    for _ in range(50):
        f_new = (-2.0 * math.log10(eps_D / 3.7 + 2.51 / (Re * math.sqrt(f)))) ** -2
        if abs(f_new - f) < 1e-10:
            break
        f = f_new
    return f_new


def pipe_friction_head(flow_m3s, pipe_id, pipe_roughness, rho, mu, length):
    """Darcy-Weisbach friction head [m]. Also returns velocity, Re, f."""
    area  = math.pi * pipe_id ** 2 / 4.0
    vel   = flow_m3s / area
    Re    = rho * vel * pipe_id / mu
    eps_D = pipe_roughness / pipe_id
    f     = friction_factor_colebrook(Re, eps_D)
    h_f   = f * (length / pipe_id) * vel ** 2 / (2.0 * 9.81)
    return h_f, vel, Re, f


# ===========================================================================
# GENERIC ESP WELL CALCULATOR
# ===========================================================================

def calc_well_pump(label, depth_m, flow_ls, temp_c,
                   pump_eff, motor_eff, cable_eff, vsd_eff,
                   pipe_id_m, pipe_roughness_m,
                   wellhead_pressure_bar,
                   mode,
                   reservoir_pressure_bar=0.0):
    """
    Calculate ESP power and energy for one well (production or injection).

    Production mode
    ---------------
      TDH = h_static(lift) + h_friction + h_minor + h_wellhead_backpressure

    Injection mode
    --------------
      The pump pushes fluid against the reservoir injection pressure.
      The hydrostatic column of fluid in the injection string ASSISTS the pump
      (gravity drives fluid downward), so it is subtracted from TDH:
        TDH = h_reservoir_pressure - h_hydrostatic + h_friction + h_minor
      If TDH <= 0 no pump is needed.
    """
    g = 9.81
    MINOR_LOSS_FACTOR = 0.08  # 8% of pipe friction head

    # --- Unit conversions ---
    flow_m3s     = flow_ls / 1000.0
    flow_m3h     = flow_m3s * 3600.0
    wellhead_pa  = wellhead_pressure_bar * 1e5
    reservoir_pa = reservoir_pressure_bar * 1e5

    # --- Fluid properties ---
    rho = water_density(temp_c)
    mu  = water_dynamic_viscosity(temp_c)
    nu  = mu / rho
    Pv  = water_vapor_pressure(temp_c)

    # --- Pipe friction ---
    h_fric, vel, Re, f_darcy = pipe_friction_head(
        flow_m3s, pipe_id_m, pipe_roughness_m, rho, mu, depth_m
    )
    h_minor = MINOR_LOSS_FACTOR * h_fric

    # --- TDH ---
    if mode == "production":
        h_static   = depth_m
        h_wellhead = wellhead_pa / (rho * g)
        TDH = h_static + h_fric + h_minor + h_wellhead
    else:  # injection
        h_reservoir   = reservoir_pa / (rho * g)
        h_hydrostatic = depth_m
        h_wellhead    = wellhead_pa / (rho * g)
        TDH = h_reservoir - h_hydrostatic + h_fric + h_minor + h_wellhead

    pump_needed = TDH > 0

    # --- NPSHa ---
    if mode == "production":
        P_intake = rho * g * depth_m
    else:
        P_intake = max(0.0, rho * g * depth_m - rho * g * (h_fric + h_minor))
    NPSHa = max(0.0, (P_intake - Pv) / (rho * g))

    # --- Specific speed ---
    omega = 2 * math.pi * 50  # 50 Hz (change to 60 for 60 Hz grids)
    Ns = omega * flow_m3s ** 0.5 / (g * TDH) ** 0.75 if pump_needed else 0.0

    # --- Power ---
    if pump_needed:
        P_hydraulic  = rho * g * flow_m3s * TDH
        P_shaft      = P_hydraulic / pump_eff
        P_motor      = P_shaft     / motor_eff
        P_surface_W  = P_motor     / (cable_eff * vsd_eff)
        P_surface_kW = P_surface_W / 1000.0
    else:
        P_hydraulic = P_shaft = P_motor = P_surface_W = P_surface_kW = 0.0

    E_annual_MWh = P_surface_kW * HOURS_PER_YEAR / 1000.0
    overall_eff  = pump_eff * motor_eff * cable_eff * vsd_eff

    # -----------------------------------------------------------------------
    # REPORT FOR THIS WELL
    # -----------------------------------------------------------------------
    sep        = "=" * 62
    mode_label = "PRODUCTION WELL" if mode == "production" else "INJECTION WELL"

    print(f"\n{sep}")
    print(f"  ESP CALCULATOR -- {mode_label} : {label}")
    print(sep)

    print("\n-- INPUT PARAMETERS -------------------------------------------")
    print(f"  Pump / perforation depth    : {depth_m:>10.1f}  m")
    print(f"  Flow rate                   : {flow_ls:>10.1f}  L/s  ({flow_m3h:.1f} m3/h)")
    print(f"  Fluid temperature at pump   : {temp_c:>10.1f}  deg C")
    if mode == "production":
        print(f"  Wellhead back-pressure      : {wellhead_pressure_bar:>10.2f}  bar (gauge)")
    else:
        print(f"  Reservoir injection pressure: {reservoir_pressure_bar:>10.2f}  bar (gauge)")
        print(f"  Additional wellhead pressure: {wellhead_pressure_bar:>10.2f}  bar (gauge)")
    print(f"  Tubing inner diameter       : {pipe_id_m*1000:>10.1f}  mm")
    print(f"  Pipe roughness              : {pipe_roughness_m*1e6:>10.1f}  um")

    print("\n-- EFFICIENCY ASSUMPTIONS -------------------------------------")
    print(f"  Pump hydraulic efficiency   : {pump_eff*100:>10.1f}  %")
    print(f"  Motor efficiency            : {motor_eff*100:>10.1f}  %")
    print(f"  Cable efficiency            : {cable_eff*100:>10.1f}  %")
    print(f"  VSD/drive efficiency        : {vsd_eff*100:>10.1f}  %")
    print(f"  Overall system efficiency   : {overall_eff*100:>10.1f}  %")

    print(f"\n-- FLUID PROPERTIES (at {temp_c:.0f} deg C) ----------------------------")
    print(f"  Density                     : {rho:>10.2f}  kg/m3")
    print(f"  Dynamic viscosity           : {mu*1e3:>10.4f}  mPa.s")
    print(f"  Kinematic viscosity         : {nu*1e6:>10.4f}  mm2/s (cSt)")
    print(f"  Vapour pressure             : {Pv/1e3:>10.2f}  kPa")

    print("\n-- HYDRAULIC ANALYSIS -----------------------------------------")
    print(f"  Flow velocity in tubing     : {vel:>10.3f}  m/s")
    regime = ("turbulent" if Re >= 4000 else
              "transitional" if Re >= 2300 else "laminar")
    print(f"  Reynolds number             : {Re:>10.0f}  ({regime})")
    print(f"  Darcy friction factor       : {f_darcy:>10.5f}")
    print(f"  Pipe friction head loss     : {h_fric:>10.2f}  m")
    print(f"  Minor losses head           : {h_minor:>10.2f}  m")
    if mode == "production":
        print(f"  Static lift head            : {depth_m:>10.2f}  m  (gravity to overcome)")
        print(f"  Wellhead back-pressure head : {wellhead_pa/(rho*g):>10.2f}  m")
    else:
        print(f"  Reservoir pressure head     : {reservoir_pa/(rho*g):>10.2f}  m  (to overcome)")
        print(f"  Hydrostatic assist head     : {depth_m:>10.2f}  m  (gravity helps, subtracted)")
        print(f"  Additional wellhead head    : {wellhead_pa/(rho*g):>10.2f}  m")
    print(f"  Total Dynamic Head (TDH)    : {TDH:>10.2f}  m")
    if not pump_needed:
        print("  >> TDH <= 0: gravity injection sufficient -- no pump required")

    print("\n-- CAVITATION CHECK -------------------------------------------")
    print(f"  NPSHa at pump intake        : {NPSHa:>10.2f}  m")
    if NPSHa > 10:
        print("  Status: OK  Adequate NPSH margin (fluid will not vaporise)")
    else:
        print("  Status: WARNING -- low NPSHa, verify pump NPSHr < NPSHa")

    print("\n-- POWER BREAKDOWN --------------------------------------------")
    if pump_needed:
        print(f"  Hydraulic power (fluid)     : {P_hydraulic/1e3:>10.2f}  kW")
        print(f"  Pump shaft power            : {P_shaft/1e3:>10.2f}  kW  (/ pump eff)")
        print(f"  Motor input power           : {P_motor/1e3:>10.2f}  kW  (/ motor eff)")
        print(f"  Surface electrical power    : {P_surface_kW:>10.2f}  kW  (/ cable & VSD eff)")
    else:
        print("  No pump required -- gravity drives injection (TDH <= 0)")

    print("\n-- ENERGY CONSUMPTION -----------------------------------------")
    print(f"  Operating hours / year      : {HOURS_PER_YEAR:>10.0f}  h/yr")
    print(f"  Annual energy consumption   : {E_annual_MWh:>10.2f}  MWh/yr")

    print("\n-- PUMP SIZING INDICATOR --------------------------------------")
    if pump_needed:
        print(f"  Dimensionless specific speed: {Ns:>10.4f}  (50 Hz base)")
        if Ns < 0.2:
            pump_type = "radial-flow -- consider multistage centrifugal"
        elif Ns < 1.0:
            pump_type = "mixed-flow -- typical ESP range"
        else:
            pump_type = "axial-flow -- consider propeller stages"
        print(f"  Suggested stage type        : {pump_type}")
    else:
        print("  N/A -- no pump required")

    return {
        "label":        label,
        "mode":         mode,
        "TDH_m":        TDH,
        "pump_needed":  pump_needed,
        "hydraulic_kW": P_hydraulic / 1e3,
        "surface_kW":   P_surface_kW,
        "annual_MWh":   E_annual_MWh,
        "rho":          rho,
        "mu_mPas":      mu * 1e3,
        "vel_m_s":      vel,
        "Re":           Re,
        "NPSHa_m":      NPSHa,
        "Ns":           Ns,
    }


# ===========================================================================
# MAIN CALCULATION
# ===========================================================================

def run_esp_calculation():

    # -----------------------------------------------------------------------
    # PRODUCTION WELL
    # -----------------------------------------------------------------------
    prod = calc_well_pump(
        label                 = "Production Well",
        depth_m               = PROD_DEPTH_M,
        flow_ls               = PROD_FLOW_LS,
        temp_c                = PROD_TEMP_C,
        pump_eff              = PROD_PUMP_EFF,
        motor_eff             = PROD_MOTOR_EFF,
        cable_eff             = PROD_CABLE_EFF,
        vsd_eff               = PROD_VSD_EFF,
        pipe_id_m             = PROD_PIPE_ID_M,
        pipe_roughness_m      = PROD_PIPE_ROUGHNESS_M,
        wellhead_pressure_bar = PROD_BACKPRESSURE_BAR,
        mode                  = "production",
    )

    # -----------------------------------------------------------------------
    # INJECTION WELL
    # -----------------------------------------------------------------------
    inj = calc_well_pump(
        label                 = "Injection Well",
        depth_m               = INJ_DEPTH_M,
        flow_ls               = INJ_FLOW_LS,
        temp_c                = INJ_TEMP_C,
        pump_eff              = INJ_PUMP_EFF,
        motor_eff             = INJ_MOTOR_EFF,
        cable_eff             = INJ_CABLE_EFF,
        vsd_eff               = INJ_VSD_EFF,
        pipe_id_m             = INJ_PIPE_ID_M,
        pipe_roughness_m      = INJ_PIPE_ROUGHNESS_M,
        wellhead_pressure_bar = INJ_WELLHEAD_BAR,
        mode                  = "injection",
        reservoir_pressure_bar= INJ_RESERVOIR_BAR,
    )

    # -----------------------------------------------------------------------
    # GEOTHERMAL THERMAL PERFORMANCE  (based on production well temperatures)
    # -----------------------------------------------------------------------
    if PROD_TEMP_C <= INJ_TEMP_C:
        raise ValueError("INJ_TEMP_C (reinjection temperature) must be lower than PROD_TEMP_C.")

    rho_avg   = water_density     ((PROD_TEMP_C + INJ_TEMP_C) / 2.0)
    cp_avg    = water_specific_heat((PROD_TEMP_C + INJ_TEMP_C) / 2.0)
    delta_T   = PROD_TEMP_C - INJ_TEMP_C
    mass_flow = (PROD_FLOW_LS / 1000.0) * rho_avg        # kg/s
    P_th_kW   = mass_flow * cp_avg * delta_T / 1000.0    # kW
    P_th_MW   = P_th_kW / 1000.0
    E_th_MWh  = P_th_kW * HOURS_PER_YEAR / 1000.0        # MWh/yr
    E_th_GWh  = E_th_MWh / 1000.0

    # Total electricity consumption (both pumps)
    E_total_MWh = prod["annual_MWh"] + inj["annual_MWh"]

    # SPF -- individual and system-wide
    SPF_prod   = E_th_MWh / prod["annual_MWh"] if prod["annual_MWh"] > 0 else float("inf")
    SPF_inj    = E_th_MWh / inj["annual_MWh"]  if inj["annual_MWh"]  > 0 else float("inf")
    SPF_system = E_th_MWh / E_total_MWh         if E_total_MWh        > 0 else float("inf")

    # Economics
    revenue_EUR    = E_th_MWh             * PRICE_HEAT_EUR_MWhth
    opex_prod_EUR  = prod["annual_MWh"]   * PRICE_ELEC_EUR_MWhe
    opex_inj_EUR   = inj["annual_MWh"]    * PRICE_ELEC_EUR_MWhe
    opex_total_EUR = opex_prod_EUR + opex_inj_EUR
    net_margin_EUR = revenue_EUR - opex_total_EUR

    # CO2 savings vs. natural gas reference plant
    # Natural gas emission factor (LHV basis, IPCC 2006): 0.202 tCO2/MWh_gas
    NG_EMISSION_FACTOR = 0.202             # tCO2 / MWh_gas (LHV)
    E_gas_MWh          = E_th_MWh / GAS_BOILER_EFF   # gas energy to produce same heat
    CO2_saved_ton      = E_gas_MWh * NG_EMISSION_FACTOR  # tCO2 avoided per year
    CO2_revenue_EUR    = CO2_saved_ton * CO2_PRICE_EUR_TON
    net_margin_with_CO2_EUR = net_margin_EUR + CO2_revenue_EUR

    # -----------------------------------------------------------------------
    # GEOTHERMAL THERMAL PERFORMANCE SECTION
    # -----------------------------------------------------------------------
    sep = "=" * 62
    print(f"\n{sep}")
    print("  GEOTHERMAL THERMAL PERFORMANCE")
    print(sep)
    print(f"  Brine temperature (produced): {PROD_TEMP_C:>10.1f}  deg C")
    print(f"  Reinjection temperature     : {INJ_TEMP_C:>10.1f}  deg C")
    print(f"  Temperature drop (dT)       : {delta_T:>10.1f}  K")
    print(f"  Average fluid density       : {rho_avg:>10.2f}  kg/m3")
    print(f"  Average specific heat       : {cp_avg:>10.1f}  J/(kg.K)")
    print(f"  Mass flow rate              : {mass_flow:>10.3f}  kg/s")
    print(f"  Installed thermal power     : {P_th_kW:>10.2f}  kW  ({P_th_MW:.3f} MW)")
    print(f"  Annual thermal energy       : {E_th_MWh:>10.2f}  MWh/yr  ({E_th_GWh:.3f} GWh/yr)")

    # -----------------------------------------------------------------------
    # PRODUCTION WELL OVERVIEW
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("  OVERVIEW -- PRODUCTION WELL")
    print(sep)
    print(f"  {'Parameter':<36} {'Value':>10}  Unit")
    print("-" * 62)
    print(f"  {'TDH':<36} {prod['TDH_m']:>10.2f}  m")
    print(f"  {'Surface electrical power':<36} {prod['surface_kW']:>10.2f}  kW")
    print(f"  {'Annual electricity':<36} {prod['annual_MWh']:>10.2f}  MWh/yr")
    print(f"  {'Electricity OPEX':<36} {opex_prod_EUR:>10,.0f}  EUR/yr")
    print(f"  {'Installed thermal power':<36} {P_th_kW:>10.2f}  kW")
    print(f"  {'Annual thermal energy':<36} {E_th_MWh:>10.2f}  MWh/yr")
    print(f"  {'SPF (thermal / prod. elec.)':<36} {SPF_prod:>10.2f}  --")

    # -----------------------------------------------------------------------
    # INJECTION WELL OVERVIEW
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("  OVERVIEW -- INJECTION WELL")
    print(sep)
    print(f"  {'Parameter':<36} {'Value':>10}  Unit")
    print("-" * 62)
    print(f"  {'TDH':<36} {inj['TDH_m']:>10.2f}  m")
    if inj["pump_needed"]:
        print(f"  {'Surface electrical power':<36} {inj['surface_kW']:>10.2f}  kW")
        print(f"  {'Annual electricity':<36} {inj['annual_MWh']:>10.2f}  MWh/yr")
        print(f"  {'Electricity OPEX':<36} {opex_inj_EUR:>10,.0f}  EUR/yr")
        print(f"  {'SPF (thermal / inj. elec.)':<36} {SPF_inj:>10.2f}  --")
    else:
        print("  No pump required -- gravity injection (TDH <= 0)")
        print(f"  {'Annual electricity':<36} {'0.00':>10}  MWh/yr")
        print(f"  {'Electricity OPEX':<36} {'0':>10}  EUR/yr")

    # -----------------------------------------------------------------------
    # FULL SYSTEM SUMMARY
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("  FULL SYSTEM SUMMARY")
    print(sep)
    print(f"  {'Parameter':<34} {'Prod.':>8}  {'Inj.':>8}  {'Total':>9}")
    print("-" * 62)
    print(f"  {'Elec. power [kW]':<34} "
          f"{prod['surface_kW']:>8.2f}  "
          f"{inj['surface_kW']:>8.2f}  "
          f"{prod['surface_kW']+inj['surface_kW']:>9.2f}")
    print(f"  {'Annual electricity [MWh/yr]':<34} "
          f"{prod['annual_MWh']:>8.2f}  "
          f"{inj['annual_MWh']:>8.2f}  "
          f"{E_total_MWh:>9.2f}")
    print(f"  {'Electricity OPEX [EUR/yr]':<34} "
          f"{opex_prod_EUR:>8,.0f}  "
          f"{opex_inj_EUR:>8,.0f}  "
          f"{opex_total_EUR:>9,.0f}")
    print("-" * 62)
    print(f"  {'Thermal power [kW]':<34} "
          f"{P_th_kW:>8.2f}  "
          f"{'--':>8}  "
          f"{P_th_kW:>9.2f}")
    print(f"  {'Annual thermal energy [MWh/yr]':<34} "
          f"{E_th_MWh:>8.2f}  "
          f"{'--':>8}  "
          f"{E_th_MWh:>9.2f}")
    print(f"  {'Heat revenue [EUR/yr]':<34} "
          f"{revenue_EUR:>8,.0f}  "
          f"{'--':>8}  "
          f"{revenue_EUR:>9,.0f}")
    print("-" * 62)
    print(f"  {'Net margin [EUR/yr]':<34} "
          f"{'':>8}  "
          f"{'':>8}  "
          f"{net_margin_EUR:>9,.0f}  "
          f"{'OK' if net_margin_EUR >= 0 else 'WARNING'}")
    print("-" * 62)
    print(f"  {'CO2 saved vs. gas boiler [tCO2/yr]':<34} "
          f"{'':>8}  "
          f"{'':>8}  "
          f"{CO2_saved_ton:>9,.1f}")
    print(f"  {'CO2 credit revenue [EUR/yr]':<34} "
          f"{'':>8}  "
          f"{'':>8}  "
          f"{CO2_revenue_EUR:>9,.0f}")
    print("-" * 62)
    print(f"  {'Net margin incl. CO2 [EUR/yr]':<34} "
          f"{'':>8}  "
          f"{'':>8}  "
          f"{net_margin_with_CO2_EUR:>9,.0f}  "
          f"{'OK' if net_margin_with_CO2_EUR >= 0 else 'WARNING'}")
    print("-" * 62)
    print(f"  {'SPF -- production pump only':<34} {SPF_prod:>8.2f}  {'':>8}  {'':>9}")
    print(f"  {'SPF -- injection pump only':<34} {'':>8}  {SPF_inj:>8.2f}  {'':>9}")
    print(f"  {'SPF -- full system':<34} {'':>8}  {'':>8}  {SPF_system:>9.2f}")
    print(sep)

    # -----------------------------------------------------------------------
    # CO2 SAVINGS DETAIL
    # -----------------------------------------------------------------------
    print(f"\n{sep}")
    print("  CO2 SAVINGS vs. NATURAL GAS REFERENCE PLANT")
    print(sep)
    print(f"  Reference plant type        : natural gas boiler / heat plant")
    print(f"  Gas boiler efficiency       : {GAS_BOILER_EFF*100:>10.1f}  %")
    print(f"  NG emission factor (LHV)    : {'0.202':>10}  tCO2/MWh_gas  (IPCC 2006)")
    print(f"  Annual geothermal heat      : {E_th_MWh:>10.2f}  MWh/yr")
    print(f"  Equiv. gas energy needed    : {E_gas_MWh:>10.2f}  MWh_gas/yr")
    print(f"  CO2 avoided (saved)         : {CO2_saved_ton:>10.1f}  tCO2/yr")
    print(f"  CO2 market price            : {CO2_PRICE_EUR_TON:>10.2f}  EUR/tCO2")
    print(f"  CO2 credit revenue          : {CO2_revenue_EUR:>10,.0f}  EUR/yr")
    print(f"  Net margin (excl. CO2)      : {net_margin_EUR:>10,.0f}  EUR/yr")
    print(f"  Net margin (incl. CO2)      : {net_margin_with_CO2_EUR:>10,.0f}  EUR/yr  "
          f"{'OK' if net_margin_with_CO2_EUR >= 0 else 'WARNING'}")
    print(sep)

    return {
        "production":              prod,
        "injection":               inj,
        "thermal_kW":              P_th_kW,
        "thermal_MWh_yr":          E_th_MWh,
        "total_elec_MWh":          E_total_MWh,
        "SPF_prod":                SPF_prod,
        "SPF_inj":                 SPF_inj,
        "SPF_system":              SPF_system,
        "revenue_EUR":             revenue_EUR,
        "opex_total_EUR":          opex_total_EUR,
        "net_margin_EUR":          net_margin_EUR,
        "CO2_saved_ton":           CO2_saved_ton,
        "CO2_revenue_EUR":         CO2_revenue_EUR,
        "net_margin_with_CO2_EUR": net_margin_with_CO2_EUR,
    }


# ===========================================================================
# ENTRY POINT
# ===========================================================================

if __name__ == "__main__":
    results = run_esp_calculation()
