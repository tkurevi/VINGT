"""Uniform plot style shared by all well models (reservoir, economy, viz).
Sober, manual-ready: white bg, light grid, fixed muted palette, no flashy markers."""
import matplotlib as mpl
PALETTE = ["#1f4e79", "#a6324a", "#2a7f7f", "#7a6a1f", "#42597a", "#666666"]
NAVY, MAROON, TEAL, OLIVE, SLATE, GREY = PALETTE
OP_POINT = dict(marker="o", color="#a6324a", markersize=6, linestyle="none", zorder=6)
def apply():
    mpl.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white",
        "font.family": "DejaVu Sans", "font.size": 10,
        "axes.titlesize": 11, "axes.titleweight": "bold", "axes.labelsize": 10,
        "axes.edgecolor": "#333333", "axes.linewidth": 0.8, "axes.axisbelow": True,
        "axes.grid": True, "grid.color": "#d9d9d9", "grid.linewidth": 0.6, "grid.alpha": 0.9,
        "axes.prop_cycle": mpl.cycler(color=PALETTE),
        "lines.linewidth": 1.6, "lines.markersize": 5,
        "legend.frameon": True, "legend.framealpha": 0.92, "legend.fontsize": 9,
        "legend.edgecolor": "#cccccc",
        "figure.dpi": 140, "savefig.dpi": 150, "savefig.bbox": "tight",
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.labelsize": 9, "ytick.labelsize": 9,
    })
apply()
