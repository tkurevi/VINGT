"""
run_standalone.py  —  drive economy_VinGT.py WITHOUT Streamlit (e.g. in Spyder: F5).

Runs the model at the default Config, then sweeps the doublet average flow and
the DH return-temperature scenario, prints a comparison table, and writes a CSV.
Edit the BASE Config below or the sweep lists to taste; nothing here touches the
physics — it only calls economy_VinGT.run() repeatedly.
"""
import os, csv
import economy_VinGT as M

def metrics(cfg):
    R = M.run(cfg)
    v3, e, en, cf = R["v3"], R["eng"], R["ener"], R["cf"]
    return dict(
        scen=cfg.dh_scenario, avg=v3["avg_flow_ls"], peak=v3["peak_flow_ls"],
        Pwf=v3["Pwf_bar"], draw=v3["drawdown_bar"], whT=v3["wellhead_T_C"],
        reinj=e["reinj_T"], brine_out=e["brine_out_C"],
        deliv_MW=en["delivered_kWth"]/1000.0, deliv_MWh=en["delivered_MWh_y0"],
        esp_kW=e["esp_kW"], inj_kW=e["inj_kW"], circ_kW=e["circ_kW"], phe_m2=e["phe_area_m2"],
        NPV_keur=cf["proj_npv"]/1e3, IRR_pct=cf["proj_irr"]*100.0,
        LCOH=cf["lcoh"], payback=cf["disc_payback"],
    )

COLS = [("scen","%s"),("avg","%.1f"),("peak","%.1f"),("Pwf","%.0f"),("draw","%.0f"),
        ("whT","%.0f"),("reinj","%.0f"),("brine_out","%.1f"),("deliv_MW","%.2f"),
        ("deliv_MWh","%.0f"),("esp_kW","%.0f"),("inj_kW","%.0f"),("circ_kW","%.1f"),
        ("phe_m2","%.0f"),("NPV_keur","%.0f"),("IRR_pct","%.1f"),("LCOH","%.1f"),("payback","%s")]

def _cell(fmt, val):
    try: return (fmt % val)
    except Exception: return str(val)

def show(title, rows):
    print("\n" + title)
    print(" ".join(f"{c:>10}" for c, _ in COLS))
    for r in rows:
        print(" ".join(_cell(fmt, r[c]).rjust(10) for c, fmt in COLS))

if __name__ == "__main__":
    here = os.path.dirname(os.path.abspath(__file__))
    base = M.Config()

    print("=" * 64); print(" BASE CASE (default Config)"); print("=" * 64)
    M.print_report(M.run(base))

    a0 = base.doublet_avg_flow_ls
    flows = [round(a0 * f, 2) for f in (0.70, 0.85, 1.00, 1.15, 1.30)]
    frows = [metrics(M.Config(doublet_avg_flow_ls=a)) for a in flows]
    show("FLOW SWEEP  (dh_scenario=%s, operating_months=%d)" % (base.dh_scenario, base.operating_months_per_yr), frows)

    srows = [metrics(M.Config(dh_scenario=s)) for s in ("A", "B", "C")]
    show("SCENARIO SWEEP  (avg=%.1f L/s -> peak=%.2f L/s)" % (a0, M.run(base)["v3"]["peak_flow_ls"]), srows)

    out = os.path.join(here, "VinGT_standalone_sweep.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=[c for c, _ in COLS]); w.writeheader()
        for r in frows + srows:
            w.writerow({c: r[c] for c, _ in COLS})
    print("\nSaved CSV:", out)
