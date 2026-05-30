"""Rock + brine thermal properties vs depth for the Sava / Drava sub-basins
(Pannonian Basin, Croatia).  Rock: Jelic (1987).  Used to ground the doublet
and Alves wellbore thermal properties instead of fixed defaults."""
import math
def rock_props(depth_m, subbasin="Sava"):
    """Return (rho_rock kg/m3, c_rock J/kgK, k_therm W/mK) at depth via Jelic 1987."""
    d_km = depth_m / 1000.0
    if subbasin.lower().startswith("s"):       # Sava
        g_depth = 2.198 + 0.0524 * d_km        # rock density g/cm3 (Jelic eq.2-style)
        c_rock = 0.602 * math.exp(-1.177 * d_km) + 0.898      # J/gK (eq.16)
    else:                                       # Drava
        g_depth = 2.214 + 0.0478 * d_km
        c_rock = 0.557 * math.exp(-1.460 * d_km) + 0.908      # J/gK (eq.17)
    rho_rock = g_depth * 1000.0
    k_therm = 0.142 * g_depth ** 2.86 + 0.623   # W/mK (eq.1 + 0.623)
    return rho_rock, c_rock * 1000.0, k_therm
def brine_props(T_C, salinity_ppm=5000.0):
    """Approx brine density (kg/m3) and cp (J/kgK) at temperature T."""
    S = salinity_ppm / 1e6
    rho = 1000.0*(1-((T_C+288.94)/(508929*(T_C+68.13)))*(T_C-3.99)**2) + 700.0*S
    cp = (4185.0 - 2.8*S*1e3)                   # mild salinity reduction
    return rho, cp
