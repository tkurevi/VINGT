"""
doublet_viz.py  -  2D map-view visualization of a producer-injector doublet
(most-permeable layer).  Shows: streamlines, the advancing cold thermal front
at several times, and the swept (drainage) area.  Pattern 1+1 (collinear).

Faithful-but-compact companion to the full gringarten_doublet.py.  Uses
potential-flow superposition (source at injector, sink at producer) for the
velocity field and integrates tracer particles for the thermal front.
Styled by plot_style for uniform manual figures.
"""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import plot_style  # noqa: F401 (applies global style on import)
    from plot_style import NAVY, MAROON, TEAL, GREY
except Exception:
    NAVY, MAROON, TEAL, GREY = "#1f4e79", "#a6324a", "#2a7f7f", "#666666"


def visualize_doublet(layer, q_ls, spacing_m, reinj_T_C, T0_C,
                      phi_doublet=0.10, years=(5, 10, 20, 30),
                      well_name="well", save_dir="."):
    """layer: dict(h, k, name). Returns list of saved PNG paths."""
    h = layer["h"]; D = spacing_m
    Q = q_ls * 1e-3                              # m3/s
    phi = phi_doublet
    xi, yi = -D/2.0, 0.0                          # injector
    xp, yp = +D/2.0, 0.0                          # producer
    strength = Q / (2.0*math.pi*phi*h)            # source/sink strength per unit thickness

    def vel(x, y):
        dxi, dyi = x-xi, y-yi; ri2 = dxi*dxi+dyi*dyi+1e-6
        dxp, dyp = x-xp, y-yp; rp2 = dxp*dxp+dyp*dyp+1e-6
        u = strength*(dxi/ri2 - dxp/rp2)          # +source injector, -sink producer
        v = strength*(dyi/ri2 - dyp/rp2)
        return u, v

    # --- velocity grid for streamlines ---
    # Fixed field of view: x and y both span -1000 .. +1000 m.
    span = 1000.0
    gx = np.linspace(-span, span, 220)
    gy = np.linspace(-span, span, 220)
    GX, GY = np.meshgrid(gx, gy)
    U, V = vel(GX, GY)
    speed = np.sqrt(U*U+V*V)

    # --- trace thermal front: launch tracers from a ring around the injector,
    #     integrate position; front at time t = particle locations ---
    n_part = 240
    th = np.linspace(0, 2*math.pi, n_part, endpoint=False)
    r0 = 0.01*D
    px = xi + r0*np.cos(th); py = yi + r0*np.sin(th)
    yr = 365.25*86400.0
    tmax = max(years)*yr
    dt = tmax/4000.0
    retard = 1.0                                   # thermal retardation folded into phi here
    snaps = {}
    paths_x = [px.copy()]; paths_y = [py.copy()]
    t = 0.0; targets = sorted(set(years)); ti = 0
    while t < tmax and ti < len(targets):
        u, v = vel(px, py)
        cap = np.hypot(px-xp, py-yp) < 0.05*D   # captured by producer -> freeze
        u = np.where(cap, 0.0, u); v = np.where(cap, 0.0, v)
        px = px + u*dt; py = py + v*dt
        t += dt
        if t >= targets[ti]*yr:
            snaps[targets[ti]] = (px.copy(), py.copy()); ti += 1
        if int(t/dt) % 80 == 0:
            paths_x.append(px.copy()); paths_y.append(py.copy())

    # ============ FIGURE 1: streamlines + velocity magnitude ============
    fig1, ax = plt.subplots(figsize=(8.4, 7.4))
    lw = 1.0 + 1.6*np.tanh(speed/(np.nanmean(speed)+1e-9))
    ax.streamplot(gx, gy, U, V, color=GREY, density=1.3, linewidth=0.7, arrowsize=0.8)
    ax.scatter([xi], [yi], marker="v", s=140, color=NAVY, edgecolor="k", linewidth=0.5, zorder=6, label="injector")
    ax.scatter([xp], [yp], marker="^", s=140, color=MAROON, edgecolor="k", linewidth=0.5, zorder=6, label="producer")
    ax.set_aspect("equal"); ax.set_xlim(-span, span); ax.set_ylim(-span, span)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
    ax.set_title(f"{well_name} doublet — streamlines ({layer['name']}, k={layer['k']:.0f} mD, D={D:.0f} m)")
    ax.legend(loc="upper left")
    p1 = os.path.join(save_dir, f"{well_name}_viz_streamlines.png")
    fig1.tight_layout(); fig1.savefig(p1); plt.close(fig1)

    # ============ FIGURE 2: advancing cold front + swept area ============
    fig2, ax2 = plt.subplots(figsize=(8.4, 7.4))
    ax2.streamplot(gx, gy, U, V, color="#e2e2e2", density=1.0, linewidth=0.5, arrowsize=0.6)
    cols = [NAVY, TEAL, "#7a6a1f", MAROON]
    for k, yz in enumerate(targets):
        if yz not in snaps: continue
        fx, fy = snaps[yz]
        order = np.argsort(np.arctan2(fy-yi, fx-xi))
        fxs, fys = np.r_[fx[order], fx[order][0]], np.r_[fy[order], fy[order][0]]
        ax2.plot(fxs, fys, "-", color=cols[k % 4], lw=1.8, label=f"cold front, {yz} yr")
        if yz == targets[-1]:
            ax2.fill(fxs, fys, color=cols[k % 4], alpha=0.10)
    ax2.scatter([xi], [yi], marker="v", s=140, color=NAVY, edgecolor="k", linewidth=0.5, zorder=6, label="injector")
    ax2.scatter([xp], [yp], marker="^", s=140, color=MAROON, edgecolor="k", linewidth=0.5, zorder=6, label="producer")
    ax2.set_aspect("equal"); ax2.set_xlim(-span, span); ax2.set_ylim(-span, span)
    ax2.set_xlabel("x [m]"); ax2.set_ylabel("y [m]")
    ax2.set_title(f"{well_name} — cold-front advance & swept area  (T_inj={reinj_T_C:.0f} °C into {T0_C:.0f} °C, φ={phi:.2f})")
    ax2.legend(loc="upper left", fontsize=8)
    p2 = os.path.join(save_dir, f"{well_name}_viz_coldfront.png")
    fig2.tight_layout(); fig2.savefig(p2); plt.close(fig2)
    return [p1, p2]


if __name__ == "__main__":
    visualize_doublet(dict(h=161, k=9128, name="dolomite"), 90.0, 1400.0, 60.0, 71.0,
                      phi_doublet=0.10, well_name="ZapGT-1", save_dir=".")
    print("viz demo done")
