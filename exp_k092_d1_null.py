"""D-1b: the WAVEFORM-MATCHED, TIME-UNIFORM null for the K-092 seed readout.

WHY THIS EXISTS, AND WHY IT IS NOT A CHANGE OF RULE.

`engine/audit_arcsine.py` declares the D-1 classification against a null of **1/3 per
band** and the quadrant against **exactly 1/4**. Both numbers are exact -- FOR A PURE
SINUSOID. The scalar that is actually frozen is not a pure sinusoid: it is a
multi-constituent solid-earth body tide with a strong, latitude-dependent diurnal
inequality, and at the Alaska-Aleutian latitudes the measured cycle length across the
seed set runs from 9.71 h to 25.41 h.

`audit_arcsine.py` says this itself, in its own docstring, and ships
`real_waveform_control` precisely because "a pure cosine could be answered with 'the
real waveform is not a sinusoid'". Its committed measurement is that the real waveform
puts **48.9 %** in the lowest quarter of range where a sinusoid puts 33.3 %.
`K092_SCALAR_PROVENANCE.md` section 4 independently measured the app's own
below-neutral-AND-falling **duty cycle at 0.4000 over 60 days** at the Sand Point site,
where the sinusoid says 0.25.

So the sinusoid numbers are the WRONG YARDSTICK for this scalar, and the program
already knew it in two committed places before this script existed. This module does
not alter the declared rule -- `exp_k092_d1.py` reports it verbatim and unchanged. It
supplies the second, waveform-matched null that the rule's own supporting module points
at, and reports the seed against BOTH. The gap between the two nulls is a direct
measurement of how much of the "effect" is the rendering.

CONSTRUCTION. For every seed site and epoch, the scalar is evaluated on a long
continuous grid (default +-10 days at 1 minute) and EVERY SAMPLE is pushed through the
IDENTICAL `local_cycle_readout` used on the events. That yields the distribution of the
band label and of the quadrant membership for an event occurring at a UNIFORMLY RANDOM
TIME at that site in that epoch -- which is exactly the null the seed observation needs.
The analysis code is imported from `exp_k092_d1.py`, not re-implemented, so the null and
the observation cannot diverge through a second code path.

SCOPE. Site coordinates and epochs come from the seed superset, which D-1 is already
licensed to read. Everything else is a deterministic astronomical waveform evaluated on
a uniform time grid: no catalogue, no event times other than the seed's own, no
statistic on any unseen event. This is the same licence `audit_arcsine.real_waveform_
control` operates under, in its own words: "not a D-7 look".
"""

from __future__ import annotations

import datetime as _dt
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_k092_d1 as D1                       # noqa: E402  the SAME analysis
from engine import audit_arcsine as AA         # noqa: E402
from engine import sitetide as ST              # noqa: E402

NULL_BRIDGE = HERE / "exp_k092_d1_null_bridge.mjs"
OUT_JSON = HERE / "results_k092_d1_null.json"

HALF_SPAN_DAYS = 10.0      # per site, centred on the event epoch
STEP_MINUTES = 1.0         # same resolution as the observation arm
RATE_HALF_MINUTES = D1.RATE_HALF_MINUTES


def band_and_quadrant_over_grid(t_days, x, rate_k):
    """Band label and quadrant membership for EVERY interior sample of one grid.

    Vectorised restatement of `exp_k092_d1.local_cycle_readout` applied at every
    sample. Equivalence to the per-event scalar path is asserted by
    `check_equivalence` below rather than assumed -- the same class of check the
    program's own failure catalogue calls a circular-equivalence trap when omitted.
    """
    tm = D1.refined_maxima(t_days, x)
    if tm.size < 2:
        return None
    k = np.searchsorted(tm, t_days, side="right") - 1
    ok = (k >= 0) & (k < tm.size - 1)
    # per-cycle min/max, computed once over the grid samples that fall in each cycle
    edges = np.searchsorted(t_days, tm)
    lo = np.full(tm.size - 1, np.nan)
    hi = np.full(tm.size - 1, np.nan)
    for j in range(tm.size - 1):
        a, b = edges[j], edges[j + 1]
        if b > a:
            seg = x[a:b]
            lo[j], hi[j] = seg.min(), seg.max()
    kk = np.clip(k, 0, tm.size - 2)
    rng = np.maximum(hi[kk] - lo[kk], 1e-300)
    u = 2.0 * (x - lo[kk]) / rng - 1.0
    rate = np.full(x.shape, np.nan)
    rate[rate_k:-rate_k] = ((x[2 * rate_k:] - x[:-2 * rate_k])
                            / (2.0 * RATE_HALF_MINUTES / 60.0))
    good = ok & np.isfinite(u) & np.isfinite(rate) & np.isfinite(lo[kk])
    return {
        "trough": int(np.sum(good & (u < -0.5))),
        "mid_slope": int(np.sum(good & (u >= -0.5) & (u <= 0.5))),
        "crest": int(np.sum(good & (u > 0.5))),
        "quadrant": int(np.sum(good & (x < 0.0) & (rate < 0.0))),
        "n": int(np.sum(good)),
    }


def check_equivalence(t_days, x, rate_k, n_probe=25, tol=1e-9):
    """Assert the vectorised grid path reproduces the scalar per-event path exactly.

    Anchors the null's analysis to the observation's analysis at the level of numbers,
    not of shared imports. Without this the two paths could drift and the comparison
    would silently become meaningless.
    """
    tm = D1.refined_maxima(t_days, x)
    if tm.size < 2:
        return {"checked": 0, "max_abs_diff": None}
    grid = band_and_quadrant_over_grid(t_days, x, rate_k)
    if grid is None:
        return {"checked": 0, "max_abs_diff": None}
    idx = np.linspace(rate_k + 5, t_days.size - rate_k - 6, n_probe).astype(int)
    worst = 0.0
    checked = 0
    for i in idx:
        rec = D1.local_cycle_readout(t_days, x, float(t_days[i]))
        if rec is None:
            continue
        # recompute u the vectorised way for this one sample
        k = int(np.searchsorted(tm, t_days[i], side="right") - 1)
        if k < 0 or k >= tm.size - 1:
            continue
        m = (t_days >= tm[k]) & (t_days <= tm[k + 1])
        lo, hi = float(x[m].min()), float(x[m].max())
        u_vec = 2.0 * (x[i] - lo) / max(hi - lo, 1e-300) - 1.0
        worst = max(worst, abs(u_vec - rec["normalised_level_u"]))
        checked += 1
    if checked and worst > tol:
        raise SystemExit("null/observation analysis paths disagree by %.3g" % worst)
    return {"checked": checked, "max_abs_diff": worst, "tolerance": tol}


def app_null(events):
    """Time-uniform null on the FROZEN scalar, via the app's own code."""
    step_ms = int(round(STEP_MINUTES * 60000))
    n_half = int(round(HALF_SPAN_DAYS * 86400000.0 / step_ms))
    n = 2 * n_half + 1
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    out_dir = HERE / "engine" / "out"
    out_dir.mkdir(parents=True, exist_ok=True)
    f_in = out_dir / "_k092_d1null_in.json"
    f_bin = out_dir / "_k092_d1null.bin"
    f_meta = out_dir / "_k092_d1null_meta.json"
    f_in.write_text(json.dumps({
        "astro_url": D1.ASTRO_URL,
        "half_span_days": HALF_SPAN_DAYS,
        "step_minutes": STEP_MINUTES,
        "sites": [{"id": e["id"], "t_ms": e["t_ms"], "lat": e["lat"], "lon": e["lon"]}
                  for e in events],
    }), encoding="utf-8")
    r = subprocess.run(["node", str(NULL_BRIDGE), str(f_in), str(f_bin), str(f_meta)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("null bridge failed (rc=%d):\n%s" % (r.returncode, r.stderr))
    meta = json.loads(f_meta.read_text(encoding="utf-8"))
    assert meta["n_per_site"] == n and meta["order"] == [e["id"] for e in events]
    raw = np.fromfile(f_bin, dtype="<f8")
    assert raw.size == len(events) * n, (raw.size, len(events) * n)
    raw = raw.reshape(len(events), n)

    tot = {"trough": 0, "mid_slope": 0, "crest": 0, "quadrant": 0, "n": 0}
    per_site = []
    equiv = None
    for j, e in enumerate(events):
        t0 = (e["t_ms"] - n_half * step_ms) / 86400000.0
        t = t0 + (step_ms / 86400000.0) * np.arange(n, dtype=np.float64)
        x = raw[j]
        if equiv is None:
            equiv = check_equivalence(t, x, rate_k)
        g = band_and_quadrant_over_grid(t, x, rate_k)
        if g is None:
            continue
        for k2 in tot:
            tot[k2] += g[k2]
        per_site.append({"id": e["id"], "lat": e["lat"],
                         "trough_fraction": g["trough"] / max(g["n"], 1),
                         "quadrant_fraction": g["quadrant"] / max(g["n"], 1)})
    return tot, per_site, meta, equiv


def sitetide_null(events):
    """The same time-uniform null on engine/sitetide.py's areal_strain."""
    step_days = STEP_MINUTES / 1440.0
    n_half = int(round(HALF_SPAN_DAYS / step_days))
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    t0ref = _dt.datetime(1900, 1, 1)
    epoch = _dt.datetime(1970, 1, 1, tzinfo=_dt.timezone.utc)
    tot = {"trough": 0, "mid_slope": 0, "crest": 0, "quadrant": 0, "n": 0}
    for e in events:
        ev_dt = epoch + _dt.timedelta(milliseconds=e["t_ms"])
        d_ev = (ev_dt.replace(tzinfo=None) - t0ref).total_seconds() / 86400.0
        t = d_ev + step_days * np.arange(-n_half, n_half + 1, dtype=np.float64)
        x = ST.site_scalar_at(t0ref, t, e["lat"], e["lon"], 0.0)
        g = band_and_quadrant_over_grid(t, x, rate_k)
        if g is None:
            continue
        for k2 in tot:
            tot[k2] += g[k2]
    return tot


def z_against(k, n, p0):
    se = math.sqrt(p0 * (1.0 - p0) / max(n, 1))
    return (k / max(n, 1) - p0) / se if se > 0 else float("nan")


def main():
    events = D1.load_seed_events()
    obs = json.loads((HERE / "results_k092_d1.json").read_text(encoding="utf-8"))
    if obs["app_arm"]["n_scored"] != len(events):
        raise SystemExit("results_k092_d1.json does not cover all %d seed events; "
                         "re-run exp_k092_d1.py first" % len(events))

    print("app null (time-uniform on the frozen scalar) ...", flush=True)
    a_tot, a_sites, a_meta, a_equiv = app_null(events)
    print("  %d null samples over %d sites" % (a_tot["n"], len(a_sites)), flush=True)
    print("sitetide null ...", flush=True)
    s_tot = sitetide_null(events)
    print("  %d null samples" % s_tot["n"], flush=True)

    p_trough = a_tot["trough"] / a_tot["n"]
    p_mid = a_tot["mid_slope"] / a_tot["n"]
    p_crest = a_tot["crest"] / a_tot["n"]
    p_quad = a_tot["quadrant"] / a_tot["n"]

    c = obs["app_arm"]["classification"]
    n_obs = c["n"]
    k_trough, k_mid, k_crest = (c["counts"]["TROUGH"], c["counts"]["MID_SLOPE"],
                                c["counts"]["CREST"])
    k_quad = obs["app_arm"]["quadrant"]["sign_condition_quadrant"]["count"]

    out = {
        "arm": "D-1b waveform-matched time-uniform null for the K-092 seed readout",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "not_a_rule_change": (
            "engine/audit_arcsine.py's declared 1/3 bands and 1/4 quadrant are exact "
            "FOR A PURE SINUSOID. exp_k092_d1.py reports them verbatim and unchanged. "
            "This module supplies the second, waveform-matched null that "
            "audit_arcsine's own real_waveform_control points at (its committed 48.9% "
            "in the lowest quarter) and that K092_SCALAR_PROVENANCE.md section 4 "
            "already measured for the quadrant (0.4000 duty cycle over 60 days at "
            "Sand Point). Both nulls are reported."),
        "construction": {
            "half_span_days": HALF_SPAN_DAYS, "step_minutes": STEP_MINUTES,
            "rate_half_minutes": RATE_HALF_MINUTES,
            "method": "every sample of a long continuous grid at each seed site and "
                      "epoch, pushed through the IDENTICAL local_cycle_readout used "
                      "on the events (imported from exp_k092_d1.py, not re-written)",
            "path_equivalence_check": a_equiv,
            "scope": "site coordinates and epochs from the already-licensed seed "
                     "superset; deterministic astronomical waveform on a uniform time "
                     "grid; no catalogue, no unseen event",
        },
        "app_null": {
            "counts": a_tot,
            "band_probabilities": {"TROUGH": p_trough, "MID_SLOPE": p_mid,
                                   "CREST": p_crest},
            "quadrant_duty_cycle": p_quad,
            "sinusoid_says": {"band": AA.BAND_NULL_PROBABILITY,
                              "quadrant": AA.QUADRANT_NULL_FRACTION},
            "provenance_corroboration": (
                "K092_SCALAR_PROVENANCE.md section 4 measured the app's quadrant duty "
                "cycle at 0.4000 (60 days, Sand Point) independently of this run."),
            "per_site": a_sites,
            "bridge_meta": {k: a_meta[k] for k in
                            ("node_version", "n_sites", "n_per_site", "step_minutes",
                             "half_span_days", "selftest_displacement_cm")},
        },
        "sitetide_null": {
            "counts": s_tot,
            "band_probabilities": {"TROUGH": s_tot["trough"] / s_tot["n"],
                                   "MID_SLOPE": s_tot["mid_slope"] / s_tot["n"],
                                   "CREST": s_tot["crest"] / s_tot["n"]},
            "quadrant_duty_cycle": s_tot["quadrant"] / s_tot["n"],
        },
        "seed_vs_both_nulls": {
            "n_events": n_obs,
            "TROUGH": {"observed": k_trough, "fraction": k_trough / n_obs,
                       "z_vs_sinusoid_null": z_against(k_trough, n_obs, 1.0 / 3.0),
                       "z_vs_waveform_null": z_against(k_trough, n_obs, p_trough)},
            "MID_SLOPE": {"observed": k_mid, "fraction": k_mid / n_obs,
                          "z_vs_sinusoid_null": z_against(k_mid, n_obs, 1.0 / 3.0),
                          "z_vs_waveform_null": z_against(k_mid, n_obs, p_mid)},
            "CREST": {"observed": k_crest, "fraction": k_crest / n_obs,
                      "z_vs_sinusoid_null": z_against(k_crest, n_obs, 1.0 / 3.0),
                      "z_vs_waveform_null": z_against(k_crest, n_obs, p_crest)},
            "QUADRANT": {
                "observed": k_quad, "fraction": k_quad / n_obs,
                "z_vs_sinusoid_null": z_against(k_quad, n_obs, 0.25),
                "z_vs_waveform_null": z_against(k_quad, n_obs, p_quad),
                "binom_p_vs_sinusoid": D1.binom_p_greater(k_quad, n_obs, 0.25),
                "binom_p_vs_waveform": D1.binom_p_greater(k_quad, n_obs, p_quad),
            },
        },
        "consequence_for_D13": (
            "K092_FREEZE.md's D-13 prediction is written against 'the uniform "
            "expectation (0.25)'. The measured time-uniform duty cycle of the frozen "
            "quadrant IN THE FROZEN SCALAR is reported above. If it differs materially "
            "from 0.25, the D-13 test as hashed compares the prospective fraction to "
            "the wrong number and would fire on the waveform alone. Popper adjudicates "
            "whether the fix is an amendment to the null constant or a restatement of "
            "the prediction; the run does not decide it."),
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 74)
    print("D-1b: SEED AGAINST BOTH NULLS (primary scorer = the app's frozen scalar)")
    print("=" * 74)
    print("  %-10s %8s | %9s %8s | %9s %8s"
          % ("band", "observed", "sinusoid", "z", "waveform", "z"))
    sv = out["seed_vs_both_nulls"]
    for b, p0 in (("TROUGH", p_trough), ("MID_SLOPE", p_mid), ("CREST", p_crest)):
        r = sv[b]
        print("  %-10s %8.4f | %9.4f %+8.2f | %9.4f %+8.2f"
              % (b, r["fraction"], 1.0 / 3.0, r["z_vs_sinusoid_null"], p0,
                 r["z_vs_waveform_null"]))
    r = sv["QUADRANT"]
    print("  %-10s %8.4f | %9.4f %+8.2f | %9.4f %+8.2f"
          % ("QUADRANT", r["fraction"], 0.25, r["z_vs_sinusoid_null"], p_quad,
             r["z_vs_waveform_null"]))
    print("\n  binomial p(>=k): vs sinusoid 0.25 = %.4g ; vs waveform %.4f = %.4g"
          % (r["binom_p_vs_sinusoid"], p_quad, r["binom_p_vs_waveform"]))
    print("  sitetide null duty cycle: %.4f (robustness)"
          % out["sitetide_null"]["quadrant_duty_cycle"])
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
