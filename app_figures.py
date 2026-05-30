"""
app_figures.py  -  figure + table helpers for the web app.
Generates: well-construction schematic, DHS (PHE + pipelines) schematic,
and IPR/VLP curve + sampled table.  Uniform style via plot_style.
"""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mp
try:
    import plot_style
    from plot_style import NAVY, MAROON, TEAL, OLIVE, SLATE, GREY
except Exception:
    NAVY, MAROON, TEAL, OLIVE, SLATE, GREY = "#1f4e79","#a6324a","#2a7f7f","#7a6a1f","#42597a","#666666"


# --------------------------------------------------------------------------- 
def well_schematic(cfg, v3, layers, save_dir):
    """Vertical well-construction schematic from config/operating data."""
    zref = float(v3.get("inj_depth_m", 2000.0))
    dyn = float(v3.get("dynamic_level_m", 0.0))
    esp = float(v3.get("esp_depth_m", dyn + 50.0))
    bottom = zref * 1.05
    fig, ax = plt.subplots(figsize=(5.2, 8.6))
    # casing (two nominal strings + open/perf interval) — proportional schematic
    w_surf, w_prod = 0.42, 0.26
    surf_shoe = min(0.18*zref, 600.0)
    ax.add_patch(mp.Rectangle((-w_surf, 0), 2*w_surf, surf_shoe, fill=False, lw=2.0, ec=GREY))
    ax.add_patch(mp.Rectangle((-w_prod, 0), 2*w_prod, zref, fill=False, lw=2.0, ec=NAVY))
    # perforated / producing interval (most permeable layer span if available)
    tops = [L.get("top", None) for L in layers if isinstance(L, dict)]
    ax.add_patch(mp.Rectangle((-w_prod, zref*0.93), 2*w_prod, zref*0.12, color=MAROON, alpha=0.25))
    ax.text(w_prod+0.05, zref*0.99, "producing interval", va="center", fontsize=8, color=MAROON)
    # fluid level
    ax.axhline(dyn, color=TEAL, lw=1.4, ls="--")
    ax.text(-w_surf-0.05, dyn, f"dynamic fluid level {dyn:.0f} m", ha="right", va="center", fontsize=8, color=TEAL)
    # ESP
    ax.add_patch(mp.Rectangle((-0.10, esp), 0.20, zref*0.035, color=NAVY))
    ax.text(0.14, esp, f"ESP @ {esp:.0f} m", va="center", fontsize=8, color=NAVY)
    # annotations
    ax.text(0, -bottom*0.03, "wellhead", ha="center", fontsize=9, fontweight="bold")
    ax.text(0, zref*1.02, f"datum z_ref {zref:.0f} m", ha="center", va="top", fontsize=8)
    ax.set_ylim(bottom, -bottom*0.06); ax.set_xlim(-0.8, 0.8)
    ax.set_xticks([]); ax.set_ylabel("depth [m]")
    ax.set_title(f"Well construction — {v3.get('well_name','well')}")
    fig.tight_layout(); p = os.path.join(save_dir, "well_schematic.png")
    fig.savefig(p); plt.close(fig); return p


# --------------------------------------------------------------------------- 
def dhs_schematic(eng, v3, cfg, save_dir):
    """PHE + DH pipeline + injection line schematic with T / flow / dP labels."""
    R = eng.get("_dhs", {}) or {}
    phe = R.get("PHE", {}) if isinstance(R.get("PHE"), dict) else {}
    Tin = v3.get("wellhead_T_C", 0.0)
    Tout = eng.get("brine_out_C", 0.0)
    reinj = eng.get("reinj_T", 0.0)
    T_sup = R.get("T_HE_out", None)      # DH supply temperature (PHE outlet)
    T_ret = R.get("T_HE_in", reinj)      # DH return entering PHE
    q = v3.get("prod_flow_ls", 0.0)
    area = phe.get("A_req_m2", eng.get("phe_area_m2", 0.0))
    lmtd = phe.get("LMTD_K", float("nan"))
    U = phe.get("U_Wm2K", float("nan"))
    fig, ax = plt.subplots(figsize=(9.2, 5.2)); ax.axis("off")
    # PHE box
    ax.add_patch(mp.Rectangle((0.42, 0.30), 0.16, 0.40, fill=True, fc="#eef2f6", ec=NAVY, lw=1.8))
    ax.text(0.50, 0.72, "PHE", ha="center", fontsize=11, fontweight="bold", color=NAVY)
    ax.text(0.50, 0.135,
            f"{phe.get('plate_class','PHE')} · A={area:,.0f} m²\nLMTD={lmtd:.1f} K · U={U:,.0f} W/m²K",
            ha="center", va="center", fontsize=8.5)
    # brine (hot) side - red, left in / right out (to injection)
    ax.annotate("", xy=(0.42, 0.62), xytext=(0.16, 0.62), arrowprops=dict(arrowstyle="-|>", color=MAROON, lw=2.4))
    ax.text(0.16, 0.66, f"brine from well\n{Tin:.0f} °C · {q:.0f} L/s", fontsize=8.5, color=MAROON)
    ax.annotate("", xy=(0.84, 0.62), xytext=(0.58, 0.62), arrowprops=dict(arrowstyle="-|>", color=MAROON, lw=2.4))
    ax.text(0.85, 0.66, f"to injection\n{Tout:.0f} °C", fontsize=8.5, color=MAROON)
    # DH (cold) side - blue, counter-current
    ax.annotate("", xy=(0.16, 0.40), xytext=(0.42, 0.40), arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=2.4))
    _sup = f"{T_sup:.0f} °C" if T_sup is not None else ""
    ax.text(0.16, 0.30, f"DH supply {_sup}\n→ city ({cfg.city_distance_m/1000:.1f} km)", fontsize=8.5, color=NAVY)
    ax.annotate("", xy=(0.58, 0.40), xytext=(0.84, 0.40), arrowprops=dict(arrowstyle="-|>", color=NAVY, lw=2.4))
    ax.text(0.85, 0.30, f"DH return\n{reinj:.0f} °C", fontsize=8.5, color=NAVY)
    # pipeline / loss / circulator note
    dp = (R.get("dp_pump_Pa", 0.0) or 0.0)/1e5
    ax.text(0.5, 0.93,
            f"DH carrier DN {cfg.dh_pipe_DN_m*1000:.0f} mm · pipe loss {eng.get('pipe_loss_kW',0):.0f} kW · "
            f"circ {eng.get('circ_kW',0):.0f} kW · ΔP_circ {dp:.2f} bar",
            ha="center", fontsize=8.5, color=GREY)
    ax.text(0.5, 0.02, f"injection line DN {cfg.brine_pipe_DN_m*1000:.0f} mm · {cfg.injection_distance_m:.0f} m",
            ha="center", fontsize=8.5, color=GREY)
    ax.set_xlim(0, 1.05); ax.set_ylim(0, 1)
    ax.set_title("District-heating system — PHE, pipeline & injection line")
    fig.tight_layout(); p = os.path.join(save_dir, "dhs_schematic.png")
    fig.savefig(p); plt.close(fig); return p


# --------------------------------------------------------------------------- 
def ipr_vlp(comm, v3, cfg, save_dir, n=26, esp_pipe_id_m=0.1016):
    """Return (png_path, table_rows). IPR from reservoir; VLP simple tubing model."""
    from pvt import Pa_to_bar
    qmax = max(v3.get("prod_flow_ls", 35.0)*1.8, 10.0)
    qs = np.linspace(qmax/n, qmax, n)               # L/s
    ipr, vlp = [], []
    rho = 970.0; g = 9.81; ID = esp_pipe_id_m
    L = float(v3.get("dynamic_level_m", 300.0)); wh = 3.2
    for q in qs:
        try: pwf = Pa_to_bar(comm.Pwf_at_q(q/1000.0))
        except Exception: pwf = float("nan")
        ipr.append(pwf)
        v = (q/1000.0)/(math.pi*(ID/2)**2)
        dpf = 0.02*(L/ID)*0.5*rho*v*v/1e5            # friction (bar), f~0.02
        vlp.append(wh + rho*g*L/1e5 + dpf)           # tubing intake requirement
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(qs, ipr, color=NAVY, label="IPR (reservoir Pwf)")
    ax.plot(qs, vlp, color=MAROON, label="VLP (tubing intake req.)")
    ax.scatter([v3.get("prod_flow_ls")], [v3.get("Pwf_bar")], color=MAROON, zorder=6, label="design point")
    ax.set_xlabel("flow rate [L/s]"); ax.set_ylabel("pressure [bar]")
    ax.set_title("IPR/VLP ESP"); ax.legend()
    fig.tight_layout(); p = os.path.join(save_dir, "ipr_vlp.png")
    fig.savefig(p); plt.close(fig)
    rows = [dict(flow_Ls=round(float(q),1), IPR_Pwf_bar=round(float(a),1), VLP_intake_bar=round(float(b),1))
            for q, a, b in zip(qs, ipr, vlp)]
    return p, rows
