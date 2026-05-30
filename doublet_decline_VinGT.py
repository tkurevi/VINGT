"""
doublet_decline_VinGT.py - Gringarten & Sauty (1975) LAYERED doublet thermal model.
Compact, importable companion to gringarten_doublet-layered-vin.py.

Returns the commingled PRODUCER temperature decline T_mix(t) over the project
life for a producer/injector doublet, given an annual-AVERAGE circulation flow
and a well spacing that the USER controls (from FEFLOW / their own analysis).

Faithful to the source script:
    rho_A c_A,i = phi_i*rho_w*c_w + (1-phi_i)*rho_R*c_R           (per layer)
    retard_i    = rho_A c_A,i / (rho_w c_w)
    flow_frac_i = (k_i h_i) / sum(k_j h_j)
    lambda_i    = (rho_w c_w * rho_A c_A,i)/(K_imp * rho_R c_R,imp) * (q_i h_i / D^2)
    t_D,i(t)    = q_i t / (retard_i D^2 h_i)
    T_wD,i(t)   = sum_theta erfc( tau_D / sqrt(lambda_i (t_D - tau_D)) ) * dpsi
    T_layer,i(t)= T0_i - T_wD,i*(T0_i - T_inj)
    T_mix(t)    = sum_i flow_frac_i * T_layer,i(t)
"""
import math
import numpy as np
from scipy.special import erfc


def doublet_temperature_decline(avg_flow_ls, spacing_m, reinj_T_C, layers,
                                years=30.0,
                                rho_w=1000.0, cp_w=4184.0,
                                rock_rho=2589.9, rock_cp=931.8,
                                imp_K=2.5, imp_rho=2589.9, imp_cp=931.8,
                                n_t=361, n_theta=4000):
    """
    layers : list of dicts, each {'h':thickness_m, 'phi':porosity, 'k':perm_mD, 'T0':init_T_C}
    Returns (t_years[np.array], T_prod_C[np.array]).
    """
    h   = np.array([L['h']   for L in layers], float)
    phi = np.array([L['phi'] for L in layers], float)
    k   = np.array([L['k']   for L in layers], float)
    T0  = np.array([L['T0']  for L in layers], float)

    Q = avg_flow_ls * 1e-3                       # m3/s total
    rho_w_cw     = rho_w * cp_w
    rho_R_cR_imp = imp_rho * imp_cp
    rhoAcA = phi * rho_w_cw + (1.0 - phi) * (rock_rho * rock_cp)
    retard = rhoAcA / rho_w_cw

    frac = (k * h) / np.sum(k * h)               # transmissivity-weighted flow split
    q    = Q * frac                              # m3/s per layer
    lam  = (rho_w_cw * rhoAcA) / (imp_K * rho_R_cR_imp) * (q * h / spacing_m**2)

    theta = np.linspace(1e-5, math.pi - 1e-5, n_theta)
    tauD  = math.pi * (np.sin(theta) - theta * np.cos(theta)) / np.sin(theta)**3
    dpsi  = (theta[1] - theta[0]) / math.pi

    yr = 365.25 * 24.0 * 3600.0
    t_phys = np.linspace(0.0, years * yr, n_t)

    def TwD(i, tp):
        if tp <= 0:
            return 0.0
        tD = q[i] * tp / (retard[i] * spacing_m**2 * h[i])
        arr = tauD < tD
        if not np.any(arr):
            return 0.0
        if (not np.isfinite(lam[i])) or lam[i] > 1e4:
            return float(np.sum(arr) * dpsi)
        arg = tauD[arr] / np.sqrt(lam[i] * (tD - tauD[arr]))
        return float(min(np.sum(erfc(arg)) * dpsi, 1.0))

    Tmix = np.empty(n_t)
    for j, tp in enumerate(t_phys):
        Ti = np.array([T0[i] - TwD(i, tp) * (T0[i] - reinj_T_C)
                       for i in range(len(h))])
        Tmix[j] = float(np.sum(frac * Ti))

    return t_phys / yr, Tmix


# ===========================================================================
# Barends (2010) longitudinal-dispersion extension + layered GS/Barends view.
# Faithful compact port of gringarten_doublet-layered.py (the script that made
# the reference figure). Reuses the SAME analytic streamtube as the GS model
# above; Barends adds per-streamtube longitudinal dispersion (Eq. 4) with an
# analytic arc length L_psi = D*tauD/(pi/3) and cap+base bleeding.  Vectorised
# over streamtubes. Economics stay on the plain GS doublet_temperature_decline();
# this is a comparison view only, fed entirely with live model data.
# ===========================================================================
try:
    from numpy import trapezoid as _trap
except ImportError:                       # NumPy < 2.0
    _trap = np.trapz
YR = 365.25 * 24.0 * 3600.0


def decline_layered(avg_flow_ls, spacing_m, reinj_T_C, layers,
                    years=30.0,
                    rho_w=1000.0, cp_w=4184.0,
                    rock_rho=2589.9, rock_cp=931.8,
                    imp_K=2.5, imp_rho=2589.9, imp_cp=931.8,
                    a_L=5.0, include_underburden=True,
                    compute_barends=True,
                    n_theta=4000, n_theta_ba=1500):
    """Layered GS + (optional) Barends 2010 decline with per-layer diagnostics.
    Faithful compact port of gringarten_doublet-layered.py.  Economics use the
    plain GS doublet_temperature_decline(); this is a comparison view only.
    Returns a dict with everything the comparison figure and summary need.
    """
    h   = np.array([L['h']   for L in layers], float)
    phi = np.array([L['phi'] for L in layers], float)
    k   = np.array([L['k']   for L in layers], float)
    T0  = np.array([L['T0']  for L in layers], float)
    nL  = len(h)

    Q = avg_flow_ls * 1e-3
    rho_w_cw      = rho_w * cp_w
    rho_R_cR_perm = rock_rho * rock_cp
    rho_R_cR_imp  = imp_rho * imp_cp
    rhoAcA = phi * rho_w_cw + (1.0 - phi) * rho_R_cR_perm
    retard = rhoAcA / rho_w_cw
    frac   = (k * h) / np.sum(k * h)
    q      = Q * frac
    lam    = (rho_w_cw * rhoAcA) / (imp_K * rho_R_cR_imp) * (q * h / spacing_m**2)

    T0_bulk = float(np.sum(frac * T0))

    # per-layer advective breakthrough (pi/3 isochrone)
    t_bt_adv_s  = (math.pi / 3.0) * spacing_m**2 * h * retard / q
    t_bt_adv_yr = t_bt_adv_s / YR

    # shared streamtube
    theta = np.linspace(1e-5, math.pi - 1e-5, n_theta)
    tauD  = math.pi * (np.sin(theta) - theta * np.cos(theta)) / np.sin(theta)**3
    dpsi  = (theta[1] - theta[0]) / math.pi

    # ----- GS per-layer curves over a fine plot grid (vectorised over t) -----
    t_gs_yr = np.concatenate([np.linspace(0.0, 1.0, 50),
                              np.linspace(1.0, 100.0, 800)])
    t_gs_s  = t_gs_yr * YR

    def TwD_layer_GS_vec(i, t_s):
        tD = q[i] * t_s / (retard[i] * spacing_m**2 * h[i])           # (n_t,)
        M  = tauD[None, :] < tD[:, None]
        if lam[i] > 1e4 or not np.isfinite(lam[i]):
            return np.clip(np.sum(M, axis=1) * dpsi, 0.0, 1.0)
        denom = lam[i] * (tD[:, None] - tauD[None, :])
        ok    = M & (denom > 0)
        arg   = np.where(ok, tauD[None, :] / np.sqrt(np.where(denom > 0, denom, 1.0)),
                         np.inf)
        return np.clip(np.sum(erfc(arg), axis=1) * dpsi, 0.0, 1.0)

    T_layers_gs = np.array([T0[i] - TwD_layer_GS_vec(i, t_gs_s) * (T0[i] - reinj_T_C)
                            for i in range(nL)])
    T_mix_gs = (frac[:, None] * T_layers_gs).sum(axis=0)

    # ----- Barends per-layer T_mix (vectorised over streamtubes) -------------
    tauD_b = (np.linspace(1e-5, math.pi - 1e-5, n_theta_ba) if n_theta_ba else theta)
    tauD_b = math.pi * (np.sin(tauD_b) - tauD_b * np.cos(tauD_b)) / np.sin(tauD_b)**3 \
        if n_theta_ba else tauD
    dpsi_b = (math.pi / (n_theta_ba - 1)) / math.pi if n_theta_ba else dpsi
    D_prime = imp_K / rho_R_cR_imp

    def TwD_layer_BA_scalar(i, t_s):
        if t_s <= 0.0:
            return 0.0
        tD_arr = (q[i] / (retard[i] * spacing_m**2 * h[i])) * t_s
        arrived = tauD_b < tD_arr
        n_arr = int(arrived.sum())
        if n_arr == 0:
            return 0.0
        td = tauD_b[arrived]
        tau_phys = td * (retard[i] * spacing_m**2 * h[i] / q[i])
        x  = spacing_m * (td / (math.pi / 3.0))            # L_psi (n_arr,)
        v  = x / tau_phys                                  # constant per layer
        Dx = imp_K / rhoAcA[i] + a_L * v                   # (n_arr,) (v const -> const)
        P  = x * v / (4.0 * Dx)
        h_eff = (2.0 if include_underburden else 1.0) * (rho_R_cR_imp / rhoAcA[i])
        H  = h[i]
        sigma_min = x / (2.0 * np.sqrt(Dx * t_s))
        sig_lo = np.maximum(sigma_min, np.sqrt(P) / 50.0)
        sig_hi = np.maximum(6.0, np.sqrt(P) * 5.0)
        nM, nLft = 200, 100
        fM = np.linspace(0.0, 1.0, nM)
        sigma_M = sig_lo[:, None] + (sig_hi - sig_lo)[:, None] * fM[None, :]
        fL = np.linspace(0.0, 1.0, nLft, endpoint=False)
        sigma_L = sigma_min[:, None] + (sig_lo - sigma_min)[:, None] * fL[None, :]
        sigma = np.concatenate([sigma_L, sigma_M], axis=1)          # (n_arr, 300)
        log_outer = -(sigma - P[:, None] / sigma) ** 2
        t_resid = t_s - x[:, None] ** 2 / (4.0 * Dx[:, None] * sigma ** 2)
        t_resid = np.where(t_resid > 0, t_resid, 1e-300)
        inner = (x[:, None] ** 2 * h_eff * math.sqrt(D_prime)
                 / (8.0 * Dx[:, None] * H * sigma ** 2)) / np.sqrt(t_resid)
        integrand = np.exp(log_outer) * erfc(inner)
        Theta = (2.0 / math.sqrt(math.pi)) * _trap(integrand, sigma, axis=1)
        return float(np.clip(np.sum(Theta) * dpsi_b, 0.0, 1.0))

    def T_mix_BA(t_s):
        return float(np.sum([frac[i] * (T0[i] - TwD_layer_BA_scalar(i, t_s) * (T0[i] - reinj_T_C))
                             for i in range(nL)]))

    if compute_barends:
        t_ba_yr = np.concatenate([np.linspace(0.0, 5.0, 16), np.linspace(5.0, 100.0, 50)])
        T_mix_ba = np.array([T0_bulk if ty <= 0 else T_mix_BA(ty * YR) for ty in t_ba_yr])
    else:
        t_ba_yr = None
        T_mix_ba = None

    # ----- decline / breakthrough times (interpolated on a coarse grid) ------
    GS_PCT, BA_PCT, BT_PCT = 10.0, 10.0, 0.1
    tgt_GS = T0_bulk - GS_PCT / 100.0 * (T0_bulk - reinj_T_C)
    tgt_BA = T0_bulk - BA_PCT / 100.0 * (T0_bulk - reinj_T_C)
    tgt_BT = T0_bulk - BT_PCT / 100.0 * (T0_bulk - reinj_T_C)
    grid = np.array([1, 2, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100, 150, 200], float)
    Tgs_grid = np.interp(grid, t_gs_yr, T_mix_gs)
    Tba_grid = np.array([T_mix_BA(g * YR) for g in grid]) if compute_barends else None

    def cross(tg, Tg, target):
        below = Tg < target
        if not below.any():
            return float('nan')
        j = int(np.argmax(below))
        if j == 0:
            return float(tg[0])
        t1, t2, T1, T2 = tg[j-1], tg[j], Tg[j-1], Tg[j]
        return float(t1 + (target - T1) * (t2 - t1) / (T2 - T1))

    out = dict(
        t_gs_yr=t_gs_yr, T_layers_gs=T_layers_gs, T_mix_gs=T_mix_gs,
        t_ba_yr=t_ba_yr, T_mix_ba=T_mix_ba,
        lam=lam, flow_frac=frac, q_ls=q * 1000.0, h=h, T0=T0, retard=retard,
        t_bt_adv_yr=t_bt_adv_yr, T0_bulk=T0_bulk, reinj_T=reinj_T_C,
        a_L=a_L, include_underburden=include_underburden,
        gs_pct=GS_PCT, ba_pct=BA_PCT, bt_pct=BT_PCT,
        target_GS=tgt_GS, target_BA=tgt_BA, target_BT=tgt_BT,
        t_BT_GS=cross(grid, Tgs_grid, tgt_BT),
        t_BT_BA=(cross(grid, Tba_grid, tgt_BT) if compute_barends else float('nan')),
        t_GS_decline=cross(grid, Tgs_grid, tgt_GS),
        t_BA_decline=(cross(grid, Tba_grid, tgt_BA) if compute_barends else float('nan')),
        years=years, spacing_m=spacing_m, avg_flow_ls=avg_flow_ls,
    )
    return out


def plot_barends_comparison(d, png_path, well_label=""):
    """Render the 2-panel GS-vs-Barends layered figure (PNG-style)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nL = len(d['h'])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5),
                                   gridspec_kw={'width_ratios': [3, 2]})
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, nL))
    for i in range(nL):
        ax1.plot(d['t_gs_yr'], d['T_layers_gs'][i], color=colors[i], lw=1.2, alpha=0.7,
                 label=f"L{i+1}: h={d['h'][i]:.0f}m, q={d['q_ls'][i]:.1f}L/s")
    ax1.plot(d['t_gs_yr'], d['T_mix_gs'], color='#c0392b', lw=3.0, zorder=10,
             label='T_mix (mass-commingled, G&S)')
    if d['T_mix_ba'] is not None:
        ax1.plot(d['t_ba_yr'], d['T_mix_ba'], '--', color='#2c3e50', lw=2.0, zorder=11,
                 label='T_mix (Barends)')
    ax1.axhline(d['T0_bulk'], color='gray', ls=':', lw=1.0, label=f"T_0,bulk = {d['T0_bulk']:.1f} °C")
    ax1.axhline(d['reinj_T'], color='steelblue', ls=':', lw=1.0, label=f"T_i = {d['reinj_T']:.1f} °C")
    if not math.isnan(d['t_BT_GS']):
        ax1.axvline(d['t_BT_GS'], color='#27ae60', ls=':', lw=1.4, alpha=0.8,
                    label=f"T_mix breakthrough (G&S, {d['bt_pct']:.2f} %): t = {d['t_BT_GS']:.2f} yr")
    if d['T_mix_ba'] is not None and not math.isnan(d['t_BT_BA']):
        ax1.axvline(d['t_BT_BA'], color='#16a085', ls=':', lw=1.4, alpha=0.8,
                    label=f"T_mix breakthrough (Barends, {d['bt_pct']:.2f} %): t = {d['t_BT_BA']:.2f} yr")
    gs_lab = (f"G&S {d['gs_pct']:.1f} %-decline → {d['target_GS']:.1f} °C at t = {d['t_GS_decline']:.2f} yr"
              if not math.isnan(d['t_GS_decline']) else
              f"G&S {d['gs_pct']:.1f} %-decline → {d['target_GS']:.1f} °C (>200 yr)")
    ax1.axhline(d['target_GS'], color='orange', ls='-.', lw=0.8, alpha=0.6, label=gs_lab)
    if d['T_mix_ba'] is not None:
        ba_lab = (f"Barends {d['ba_pct']:.1f} %-decline → {d['target_BA']:.1f} °C at t = {d['t_BA_decline']:.2f} yr"
                  if not math.isnan(d['t_BA_decline']) else
                  f"Barends {d['ba_pct']:.1f} %-decline → {d['target_BA']:.1f} °C (>200 yr)")
        ax1.axhline(d['target_BA'], color='#8e44ad', ls='-.', lw=0.8, alpha=0.6, label=ba_lab)
    ax1.set_xlabel('Time [years]'); ax1.set_ylabel('Temperature [°C]')
    ax1.set_title('Per-layer T_layer,i(t) and commingled T_mix(t) at top of pay')
    ax1.legend(loc='lower left', fontsize=7, ncol=2, framealpha=0.9)
    ax1.grid(alpha=0.3); ax1.set_xlim([0, 100])

    y = np.arange(nL)
    ax2.barh(y, d['t_bt_adv_yr'], height=0.6, color=colors, edgecolor='black', lw=0.5)
    for i in range(nL):
        ax2.text(d['t_bt_adv_yr'][i] + 0.5, i,
                 f"{d['t_bt_adv_yr'][i]:.1f} yr  (λ={d['lam'][i]:.2f}, q={d['flow_frac'][i]*100:.1f}%)",
                 va='center', fontsize=8)
    ax2.axvline(d['years'], color='red', ls='--', lw=1.5, alpha=0.7, label=f"Δt = {d['years']:.0f} yr")
    ax2.set_yticks(y); ax2.set_yticklabels([f"L{i+1}" for i in range(nL)])
    ax2.set_xlabel('Advective breakthrough time [years]'); ax2.set_title('Per-layer breakthrough times')
    ax2.legend(loc='lower right', fontsize=9); ax2.grid(axis='x', alpha=0.3); ax2.invert_yaxis()
    ax2.set_xlim([0, max(d['t_bt_adv_yr'].max() * 1.4, d['years'] * 1.3)])
    a_txt = f"a_L = {d['a_L']:.1f} m, {'cap+base' if d['include_underburden'] else 'cap-only'} bleeding"
    ttl = (f"Layered-aquifer doublet{(' — ' + well_label) if well_label else ''} — {nL} permeable layers\n"
           f"D = {d['spacing_m']:.0f} m, Q_avg = {d['avg_flow_ls']:.1f} L/s, Δt = {d['years']:.0f} yr   ({a_txt})")
    fig.suptitle(ttl, fontsize=11, y=1.00)
    fig.tight_layout()
    fig.savefig(png_path, dpi=220, bbox_inches='tight')
    plt.close(fig)
    return png_path


if __name__ == "__main__":
    # VinGT-1 layers (from V3) - quick self-test
    layers = [dict(h=20, phi=0.092, k=27.8, T0=95.5),
              dict(h=15, phi=0.064, k=0.50, T0=101.0),
              dict(h=30, phi=0.079, k=3.37, T0=104.7),
              dict(h=35, phi=0.031, k=0.32, T0=117.7),
              dict(h=25, phi=0.064, k=0.64, T0=115.8),
              dict(h=50, phi=0.055, k=1.42, T0=125.6)]
    t, T = doublet_temperature_decline(11.0, 450.0, 60.0, layers, years=30)
    print("year   T_prod(C)")
    for yr in (0, 5, 10, 20, 30):
        idx = np.argmin(abs(t - yr))
        print(f"  {t[idx]:4.0f}   {T[idx]:6.2f}")
