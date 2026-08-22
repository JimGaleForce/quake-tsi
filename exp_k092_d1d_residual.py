"""D-1d: THE ORTHOGONAL COMPLEMENT of the D-1 bound, on the K-092 seed set. Price 0.

WHY THIS ARM EXISTS, AND WHY IT IS NOT A RESCUE.

K-108 established that D-1c's resolved-Coulomb statistics are ~89 % variance-shared
with the scalar D-1 had already bounded (`corr = -0.9506`, verified at machine
precision), so D-1c was not independent evidence and its apparent second signal was the
first one in anti-correlated coordinates. The obvious response is to discard the arm.

**That response throws away the 11 %, and the 11 % is the interesting part.** A bound on
one direction is not a bound on the space. The right operation is PROJECTION, not
rejection: split every candidate into the part the bounded axis explains and the part it
cannot, and test only the second. This arm is that test.

The measurement that makes it worth running is Kepler's, re-derived by the supervisor:
the D-1 ceiling covers the body tide's isotropic content at `r = 1.000000` (an
identity), normal stress at every bearing at |r| = 0.949 to 0.952, and resolved thrust
Coulomb at r = -0.9506 -- but it covers stressing rate at r = -0.000022 (orthogonal by
the theorem that the derivative of a stationary process is uncorrelated with the
process at lag zero) and horizontal shear at r = +0.0098. **The directions D-1 did not
bound are exactly the directions the "rope" lives in, and that was true before the seed
was read.**

---------------------------------------------------------------------------
THE DECLARATION. Fixed BEFORE any number is computed.
---------------------------------------------------------------------------

FAULT GEOMETRY: strike 250, dip 20, rake 90, friction 0.4 -- the same declared
megathrust assumption as D-1c. NOTE, per K-111: at the free surface friction and dip
are algebraically inert, so no grid over them appears here. That is a correction of
D-1c's design, not a new choice.

THE BOUNDED AXIS that is projected out of every candidate: `sitetide.areal_strain`, the
program's declared SCALAR_FOR_PHASE, which is a positive multiple of the same W2 as the
frozen app scalar and therefore shares its quadrant exactly.

THE BATTERY, eight statistics, all chosen for ORTHOGONALITY rather than for physical
appeal, which is the inversion K-098 argues for:

  R1 RESID_CFS_TROUGH  residual Coulomb in the lowest third of its own local cycle
  R2 RESID_CFS_CREST   ... highest third
  R3 RATE_RISING       d(scalar)/dt > 0 at the event. Orthogonal BY IDENTITY, and the
                       exact quadrature complement of the level statistic D-1 bounded
  R4 RATE_HIGH         |d(scalar)/dt| in the top third of its own local cycle -- the
                       "maximum unloading stress rate" of §P7-23(C), as a CONTINUOUS
                       covariate rather than a sign, which Popper recorded as his own
                       omission at §P7-25(5) K-096 item 3
  R5 SHEAR_TROUGH      horizontal shear along the fault-normal bearing, lowest third
  R6 RAKE0_CREST       Coulomb at rake 0 -- a near-orthogonal receiver (K-112), highest
                       third
  R7 BEARING_RESID_R2  axial concentration of the bearing AFTER the scalar-explained
                       part is removed, on the doubled angle
  R8 ROPE_RESID_EASING residual fault-normal traction positive and falling: the rope,
                       on the part of the traction the bound does not cover

THE NULL: waveform-matched and time-uniform, identical in construction to D-1b and
D-1c. Every sample of a long grid at each seed site and epoch is pushed through the
same feature code; the event is the centre sample; the rest are that site's null pool.
A residual is a DIFFERENT statistic from the raw one and gets its OWN null -- scoring a
residual against the raw statistic's null would be a fresh instance of exactly the
error D-1 was killed by.

MULTIPLICITY: max-statistic calibration over the whole battery via
`engine/dwell_null.calibrate`, with the battery's EFFECTIVE RANK reported beside its
nominal count per Kepler's proposed SP-9. Power floor asserted before the run.

NOT EVIDENCE, AND THE ASYMMETRY IS UNCHANGED. The seed set was selected by having been
looked at. A null here extends the bound into the orthogonal directions, which is worth
having. A hit here is NOT a finding; it is a reason to build the lattice and pay for it.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_k092_d1 as D1                        # scope, loader, cycle machinery
from engine import dwell_null as DN
from engine import residual as RS
from engine import sitetide as ST
from engine import tidal_tensor as TT

OUT_JSON = HERE / "results_k092_d1d_residual.json"

FAULT = {"strike_deg": 250.0, "dip_deg": 20.0, "rake_deg": 90.0}
FRICTION = 0.4
HALF_SPAN_DAYS = 10.0
STEP_MINUTES = 2.0
RATE_HALF_MINUTES = 10.0
N_NULL_REPLICATES = 2000
RNG_SEED = 20260822

STAT_NAMES = ("RESID_CFS_TROUGH", "RESID_CFS_CREST", "RATE_RISING", "RATE_HIGH",
              "SHEAR_TROUGH", "RAKE0_CREST", "BEARING_RESID_R2",
              "ROPE_RESID_EASING")
SIMPLE = tuple(n for n in STAT_NAMES if n != "BEARING_RESID_R2")


def local_cycle_u(t, x):
    """u within each bracketing-maxima cycle. The D-1b construction, reused."""
    tm = D1.refined_maxima(t, x)
    n = x.size
    if tm.size < 2:
        return np.full(n, np.nan)
    k = np.searchsorted(tm, t, side="right") - 1
    ok = (k >= 0) & (k < tm.size - 1)
    edges = np.searchsorted(t, tm)
    lo = np.full(tm.size - 1, np.nan)
    hi = np.full(tm.size - 1, np.nan)
    for j in range(tm.size - 1):
        a, b = edges[j], edges[j + 1]
        if b > a:
            lo[j], hi[j] = x[a:b].min(), x[a:b].max()
    kk = np.clip(k, 0, tm.size - 2)
    u = 2.0 * (x - lo[kk]) / np.maximum(hi[kk] - lo[kk], 1e-300) - 1.0
    return np.where(ok, u, np.nan)


def site_features(t, lat, lon, rate_k):
    """Every declared feature at every sample, plus the orthogonality table."""
    jd = t + 2440587.5
    scal = ST.site_scalar(jd, lat, lon, 0.0)
    st = TT.stress_tensor(jd, lat, lon, 0.0)
    strike = FAULT["strike_deg"]

    cfs = TT.coulomb(st, strike, FAULT["dip_deg"], FAULT["rake_deg"],
                     FRICTION)["coulomb_pa"]
    cfs0 = TT.coulomb(st, strike, FAULT["dip_deg"], 0.0, FRICTION)["coulomb_pa"]
    shear = TT.shear_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], strike + 90.0)
    rope = TT.normal_along_bearing(st["s_NN"], st["s_EE"], st["s_NE"], strike + 90.0)
    bearing = TT.principal_bearing(st["e_NN"], st["e_EE"], st["e_NE"])

    # PROJECT OUT THE BOUNDED AXIS. This is the whole point of the arm.
    resid_cfs, sh_cfs = RS.project_out(cfs, [scal])
    resid_rope, sh_rope = RS.project_out(rope, [scal])
    resid_shear, sh_shear = RS.project_out(shear, [scal])
    resid_cfs0, sh_cfs0 = RS.project_out(cfs0, [scal])
    ax = RS.axial_shared_variance(bearing, [scal])

    n = scal.size
    rate = np.full(n, np.nan)
    rate[rate_k:-rate_k] = ((scal[2 * rate_k:] - scal[:-2 * rate_k])
                            / (2.0 * RATE_HALF_MINUTES / 60.0))
    d_rope = np.full(n, np.nan)
    d_rope[rate_k:-rate_k] = ((resid_rope[2 * rate_k:] - resid_rope[:-2 * rate_k])
                              / (2.0 * RATE_HALF_MINUTES / 60.0))

    u_cfs = local_cycle_u(t, resid_cfs)
    u_shear = local_cycle_u(t, resid_shear)
    u_cfs0 = local_cycle_u(t, resid_cfs0)
    u_rate = local_cycle_u(t, np.nan_to_num(np.abs(rate), nan=0.0))

    ok = np.zeros(n, dtype=bool)
    margin = rate_k + int(round(1.5 * 1440.0 / STEP_MINUTES))
    ok[margin:n - margin] = True
    ok &= np.isfinite(u_cfs) & np.isfinite(u_shear) & np.isfinite(u_cfs0)
    ok &= np.isfinite(u_rate) & np.isfinite(rate) & np.isfinite(d_rope)

    # residual bearing, reconstructed from the residual doubled-angle components
    b_resid = np.mod(0.5 * np.degrees(np.arctan2(ax["residual_sin"],
                                                 ax["residual_cos"])), 180.0)
    return {
        "ok": ok,
        "RESID_CFS_TROUGH": u_cfs < -0.5,
        "RESID_CFS_CREST": u_cfs > 0.5,
        "RATE_RISING": rate > 0.0,
        "RATE_HIGH": u_rate > 0.5,
        "SHEAR_TROUGH": u_shear < -0.5,
        "RAKE0_CREST": u_cfs0 > 0.5,
        "ROPE_RESID_EASING": (resid_rope > 0.0) & (d_rope < 0.0),
        "_ANG2": TT.axial_to_circular(b_resid),
        "_shared": {"cfs_thrust": sh_cfs, "rope_traction": sh_rope,
                    "horizontal_shear": sh_shear, "cfs_rake0": sh_cfs0,
                    "bearing_axial": ax["shared_variance"]},
        "_raw_for_rank": {"cfs": cfs, "cfs0": cfs0, "shear": shear, "rope": rope,
                          "rate": np.nan_to_num(rate), "scalar": scal},
    }


def battery(mats, idx):
    rows = np.arange(mats["RATE_RISING"].shape[0])
    if idx.ndim == 1:
        take = lambda m: m[rows, idx]
        red = lambda v: float(v.mean())
    else:
        take = lambda m: m[rows[:, None], idx]
        red = lambda v: v.mean(axis=0)
    out = {k: red(take(mats[k])) for k in SIMPLE}
    c = red(take(mats["_COS2"]))
    s = red(take(mats["_SIN2"]))
    out["BEARING_RESID_R2"] = (float(np.hypot(c, s)) if idx.ndim == 1
                               else np.hypot(c, s))
    return out


def main():
    floor = DN.assert_power_floor(len(STAT_NAMES), N_NULL_REPLICATES, alpha=0.05)
    events = D1.load_seed_events()
    rng = np.random.default_rng(RNG_SEED)
    step_days = STEP_MINUTES / 1440.0
    n_half = int(round(HALF_SPAN_DAYS / step_days))
    rate_k = int(round(RATE_HALF_MINUTES / STEP_MINUTES))
    print("D-1d: %d seed events, power floor %.4f" % (len(events), floor["floor"]),
          flush=True)

    feats, shared_rows, rank_rows = [], [], None
    for e in events:
        t = (e["t_ms"] / 86400000.0
             + step_days * np.arange(-n_half, n_half + 1, dtype=np.float64))
        f = site_features(t, e["lat"], e["lon"], rate_k)
        if not f["ok"][n_half]:
            raise SystemExit("event %s not scoreable at the grid centre" % e["id"])
        shared_rows.append(f.pop("_shared"))
        raw = f.pop("_raw_for_rank")
        if rank_rows is None:
            rank_rows = RS.orthogonality_report(
                {k: v for k, v in raw.items() if k != "scalar"}, [raw["scalar"]])
        feats.append(f)
    print("  all %d events scoreable" % len(feats), flush=True)

    mats = {k: np.stack([f[k].astype(np.float32) for f in feats]) for k in SIMPLE}
    ang = np.stack([f["_ANG2"] for f in feats])
    mats["_COS2"] = np.cos(ang).astype(np.float32)
    mats["_SIN2"] = np.sin(ang).astype(np.float32)

    n_sites = len(feats)
    obs = battery(mats, np.full(n_sites, n_half, dtype=np.int64))
    pool = DN.time_uniform_pool(feats[0]["ok"].size,
                                rate_k + int(round(1.5 * 1440.0 / STEP_MINUTES)))
    idx = DN.draw_replicates(pool, n_sites, N_NULL_REPLICATES, rng)
    bad = sum(int((~f["ok"][idx[s]]).sum()) for s, f in enumerate(feats))
    if bad:
        raise SystemExit("null pool hit %d invalid samples" % bad)
    nb = battery(mats, idx)
    cal = DN.calibrate(obs, {k: np.asarray(nb[k], float) for k in STAT_NAMES},
                       STAT_NAMES)

    pooled_shared = {k: float(np.mean([r[k] for r in shared_rows]))
                     for k in shared_rows[0]}
    out = {
        "arm": "D-1d orthogonal-complement arm on the K-092 seed set",
        "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "priced_tests": 0,
        "scope": {"input": "K092_seed_exclusion_superset.csv ONLY",
                  "sha256_verified": D1.SEED_CSV_SHA256, "n_events": len(events),
                  "licence": "K-096 first-move licence, second view; already-seen set, "
                             "excluded from D-12 by construction; D-7 remains STRUCK"},
        "why_not_a_rescue": (
            "K-108 showed D-1c was ~89% variance-shared with the axis D-1 bounded. "
            "Discarding it throws away the 11%. This arm PROJECTS OUT the bounded axis "
            "and tests only what remains, which is a graded answer instead of a binary "
            "one. The orthogonality of the rate axis is a theorem, not a discovery "
            "made after the seed was read."),
        "not_evidence": (
            "the seed set was selected by having been looked at. A null here EXTENDS "
            "the bound into the orthogonal directions. A hit here is NOT a finding."),
        "declaration": {
            "fault": FAULT, "friction": FRICTION,
            "no_friction_or_dip_grid_because": (
                "K-111: at the free surface friction and dip are algebraically inert, "
                "so D-1c's 45-cell grid was 5 statistics. Gridding them again would "
                "repeat a known-vacuous design."),
            "bounded_axis_projected_out": "sitetide.areal_strain (SCALAR_FOR_PHASE)",
            "statistics": list(STAT_NAMES),
            "half_span_days": HALF_SPAN_DAYS, "step_minutes": STEP_MINUTES,
            "rate_half_minutes": RATE_HALF_MINUTES,
            "n_null_replicates": N_NULL_REPLICATES, "rng_seed": RNG_SEED,
            "power_floor": floor,
            "null": "waveform-matched time-uniform, engine/dwell_null.py; the RESIDUAL "
                    "gets its OWN null and does not inherit the raw statistic's",
            "scope_flags": TT.TENSOR_SCOPE_FLAGS,
        },
        "orthogonality": {
            "pooled_shared_variance_with_bounded_axis": pooled_shared,
            "pooled_new_information_fraction": {k: 1.0 - v
                                                for k, v in pooled_shared.items()},
            "per_site_spread": {
                k: {"min": float(np.min([r[k] for r in shared_rows])),
                    "max": float(np.max([r[k] for r in shared_rows]))}
                for k in shared_rows[0]},
            "battery_effective_rank": rank_rows["battery_effective_rank"],
            "ranked_by_novelty": rank_rows["ranked_by_novelty"],
            "sp9_note": rank_rows["rule"],
        },
        "calibration": cal,
    }
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")

    print("\n" + "=" * 76)
    print("D-1d: THE ORTHOGONAL COMPLEMENT OF THE D-1 BOUND")
    print("=" * 76)
    print("  new-information fraction vs the bounded axis (1 = fully orthogonal):")
    for k, v in sorted(pooled_shared.items(), key=lambda kv: kv[1]):
        print("    %-18s %.4f" % (k, 1.0 - v))
    print("  battery effective rank: %d of %d nominal"
          % (rank_rows["battery_effective_rank"]["effective_rank"],
             rank_rows["battery_effective_rank"]["n_nominal"]))
    print()
    for k in STAT_NAMES:
        r = cal["per_statistic"][k]
        print("  %-18s obs %8.4f  null %8.4f +- %7.4f   z %+6.2f"
              % (k, r["observed"], r["null_mean"], r["null_sd"], r["z"]))
    print("\n  max|z| = %.3f ; null 95th pct %.3f ; MAX-STATISTIC p = %.4f"
          % (cal["max_abs_z_observed"], cal["max_abs_z_null_p95"],
             cal["max_statistic_p"]))
    print("\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
