"""
wellbore_T.py  -  Bottomhole -> wellhead brine temperature for a producing well.
Ramey (1962) conduction + Alves et al. (1992) elevation term + CO2 cooling
(exsolution + Joule-Thomson).  GWR is treated entirely as CO2 (conservative).

Wrapped from the uploaded wellbore_T_CO2.py into a single callable that returns
the wellhead temperature, used by economy*.py get_v3_results() to set year-0 T.
"""
import math
import numpy as np

R_GAS = 8.314
M_CO2, M_H2O = 0.04401, 0.01802
G = 9.81

def _p_sat_water(T_K):
    Tc, Pc = 647.096, 22.064
    if T_K >= Tc: return Pc
    t = 1.0 - T_K/Tc
    a = (-7.85951783*t + 1.84408259*t**1.5 - 11.7866497*t**3
         + 22.6807411*t**3.5 - 15.9618719*t**4 + 1.80122502*t**7.5)
    return Pc*math.exp(Tc/T_K*a)

def _henry_CO2(T_K):
    Tr = T_K/647.096; tau = 1.0 - Tr
    a, b, c = -9.4234, 4.0087, 10.3199
    ln_ratio = a/Tr + b*tau**0.355/Tr + c*math.exp(tau)*Tr**(-0.41)
    return _p_sat_water(T_K)*math.exp(ln_ratio)        # MPa

def _m_sat(T_K, P_MPa):
    kH = _henry_CO2(T_K)
    P_CO2 = max(0.0, P_MPa - _p_sat_water(T_K))
    x = min(P_CO2/kH, 0.3)
    mol = (x/(1.0-x))/M_H2O                              # mol CO2 / kg water
    return mol*M_CO2                                     # kg CO2 / kg water

def _dH_exsol(T_K, dT=2.0):
    return R_GAS*T_K**2*(math.log(_henry_CO2(T_K+dT))-math.log(_henry_CO2(T_K-dT)))/(2*dT)

def wellhead_temperature(L_m, q_ls, T_bh_C, geo_grad_K_m, T_surf_C=12.0,
                         k_e=2.5, alpha_e=1.0e-6, r_wb=0.125,
                         rho_f=975.0, cp_f=4185.0,
                         P_bh_bar=150.0, P_wh_bar=3.2,
                         GWR_m3m3=0.0, eta_JT=1.0, cp_CO2=900.0,
                         dyn_level_m=0.0, U_dry_frac=0.10,
                         prod_time_yr=1.0, n=200, return_detail=False):
    """Return wellhead temperature (degC).  RAMEY conduction + CO2/J-T (NO elevation term).
    Two-zone U: above the dynamic fluid level the tubing-casing annulus is gas-filled
    (insulating) so conduction loss is scaled by U_dry_frac; below it, full conduction."""
    m_dot = q_ls*1e-3*rho_f
    t = prod_time_yr*365.25*86400.0
    f_t = math.log(2.0*math.sqrt(alpha_e*t)/r_wb) - 0.29
    f_t = max(f_t, 0.5)
    A = (m_dot*cp_f)/(2.0*math.pi*k_e)*f_t               # Ramey relaxation length [m]
    w_tot = GWR_m3m3*1.842/999.0                          # kg CO2/kg water (GWR all-CO2)

    s = np.linspace(L_m, 0.0, n)                          # depth, bottom -> top
    P = np.linspace(P_bh_bar, P_wh_bar, n)/10.0           # MPa, linear
    T = T_bh_C
    diss_prev = _m_sat(T+273.15, P[0]) if w_tot > 0 else 0.0
    diss_prev = min(diss_prev, w_tot)
    dT_cond = dT_co2 = 0.0
    for i in range(1, n):
        dz = s[i-1]-s[i]                                  # >0, distance up
        T_e = T_surf_C + geo_grad_K_m*s[i]
        zone = U_dry_frac if s[i] < dyn_level_m else 1.0  # gas-insulated above fluid level
        d_cond = -(T - T_e)/A*dz*zone                     # Ramey conduction (two-zone U)
        d_co2 = 0.0
        if w_tot > 0:
            sat = min(_m_sat(T+273.15, P[i]), w_tot)
            d_exsol = max(diss_prev - sat, 0.0)           # mass exsolved this step
            free = max(w_tot - sat, 0.0)
            dHex = _dH_exsol(T+273.15)/M_CO2              # J/kg CO2
            d_co2 -= d_exsol*dHex/cp_f                    # endothermic cooling
            dP = (P[i-1]-P[i])                            # MPa drop
            d_co2 -= free*cp_CO2*eta_JT*dP/cp_f           # J-T cooling of free gas
            diss_prev = sat
        T += d_cond + d_co2
        dT_cond += d_cond; dT_co2 += d_co2
    if return_detail:
        return T, dict(cond=dT_cond, co2=dT_co2, A=A, total_drop=T_bh_C-T)
    return T


if __name__ == "__main__":
    for name, L, q, Tbh, g, GWR, dyn in [
        ("VinGT-1", 2080, 15, 100, 0.040, 0.40, 706),
        ("VGGT-1",  2211, 35, 101, 0.043, 0.89, 430),
        ("OsGT-1",  1515, 35, 88.5, 0.0466, 1.37, 549),
        ("ZapGT-1", 1617, 125, 71, 0.0365, 1.0, 120),
        ("ZapV2",   2047, 125, 86.7, 0.0365, 1.0, 120)]:
        Twh, d = wellhead_temperature(L, q, Tbh, g, GWR_m3m3=GWR, P_bh_bar=Tbh*1.5+30, dyn_level_m=dyn, return_detail=True)
        print(f"{name}: T_bh {Tbh} -> T_wh {Twh:.1f} C  (cond {d['cond']:.1f}, CO2 {d['co2']:.1f})")
