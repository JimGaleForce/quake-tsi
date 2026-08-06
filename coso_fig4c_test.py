"""Exact reproduction of Lu et al. (2025) Figure 4c per Weifan Lu's 2026-08-06 email.

Frozen protocol: COSO_FIG4C_PROTOCOL.md (hash in download_log.md). Bin lat [36.2, 36.6],
lon [-118.0, -117.6] (centered -117.8, 36.4); SCSN declustered M>=1.5; shear stress (tau)
resolved on the per-event YHS focal-mechanism plane; significance vs 3,000 orientation-
preserving time-randomized synthetics. Machinery reused verbatim from coso_fm_test.py;
only the box, seed, and a QTM M>=1.5 sensitivity cut differ.
"""
import json
import importlib.util
from pathlib import Path

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("fmtest", HERE / "coso_fm_test.py")
fmtest = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fmtest)

FIG4C_BOX = dict(lat=(36.2, 36.6), lon=(-118.0, -117.6))
fmtest.COSO = FIG4C_BOX
fmtest.SEED = 20260806

if __name__ == "__main__":
    out = {"box": {k: list(v) for k, v in FIG4C_BOX.items()},
           "protocol": "COSO_FIG4C_PROTOCOL.md", "seed": fmtest.SEED}

    # Primary (tau) + secondary (sigma_n), native catalogs (SCSN file is already M>=1.5)
    for comp in ["tau", "sigma_n"]:
        print(f"=== component: {comp} (native catalogs) ===")
        out[comp] = fmtest.run(comp)

    # Pre-declared sensitivity: QTM cut to M>=1.5 to match the SCSN threshold
    orig_load = fmtest.cc.load_declustered

    def load_m15(fname):
        df = orig_load(fname)
        return df[df.mag >= 1.5] if "QTM" in fname else df

    fmtest.cc.load_declustered = load_m15
    print("=== component: tau (QTM cut to M>=1.5) ===")
    out["tau_qtm_m1.5"] = {"QTM": fmtest.run("tau")["QTM"]}
    fmtest.cc.load_declustered = orig_load

    (HERE / "results_coso_fig4c.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    print("-> results_coso_fig4c.json")
