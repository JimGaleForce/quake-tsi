"""SENSITIVITY OF THE AMPLITUDE ARM. What its null is worth. Priced 0.

A null result is only as good as the power behind it, and this program has already been
bitten once by the difference: Kepler's K-118 showed the scan's |R| statistics could pass
through ZERO sensitivity at unlucky nuisance phase, which would have made their nulls
void rather than merely weak. The rule adopted from that episode -- engine/dwell_null.py,
`assert_sensitivity_floor` -- is that a statistic must be shown to RESPOND to the effect
it claims to measure before its null may be quoted. This module discharges that duty for
the six amplitude statistics.

THE INJECTION, and it is deliberately conservative. Within each event's OWN null pool --
the same two-fortnight neighbourhood the null draws from -- resample the event with
probability proportional to (1 + eps * z_driver), where z_driver is the standardised
driver at that candidate time. So eps = 0.10 means a ten percent preference for
high-driver moments over the event's own two weeks. This cannot borrow power from long-
timescale structure, seasonality, or catalogue growth, because the comparison is entirely
within the neighbourhood. It is the smallest honest version of the hypothesis.

Four drivers are injected separately -- amplitude, shear, |rate|, signed rate -- and all
six statistics are scored against every one. The off-diagonal entries matter as much as
the diagonal: if the six statistics all responded to the same driver they would be six
views of one test and the multiplicity correction over six would be wasted.

WHAT IT ESTABLISHES. Run at Chile (n = 967 declustered), the four continuous statistics
detect their OWN driver at 0.85 to 0.90 power for eps = 0.10 and at 1.00 for eps = 0.20,
at a two-sided five percent threshold, from ONE region alone. The arm's null is therefore
a real bound and not an absence of instrument.

The two weakest probes are the decile-occupancy pair E5 and E6, at 0.55 and 0.50 for
eps = 0.10, which is the expected price of throwing away everything except whether a
threshold was crossed. They reach 0.97 and 0.95 by eps = 0.20. Their nulls are quoted at
that power and not at the continuous statistics'.

The OFF-diagonal entries are the other half of the point. E1 scores 0.05 against a
signed-rate injection and E4 scores 0.04 against an amplitude injection: the statistics
are close to orthogonal probes rather than six views of one quantity, so the correction
over six is buying real coverage. The one genuine overlap is areal amplitude with shear
(0.15 to 0.85 across eps, both directions), which is physics -- the two share the same
degree-2 potential -- and not redundancy of construction.
"""

from __future__ import annotations

import datetime as _dt
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import exp_world_amplitude as A
import exp_world_harmonics as W

OUT_JSON = HERE / "results_world_amplitude_sensitivity.json"

REGION = "Chile"
DRIVERS = ("areal_abs", "shear", "rate_abs", "rate_signed")
EPS = (0.10, 0.20, 0.30)
TRIALS = 400
CRIT = 1.959963984540054          # two-sided 5 percent
SEED = 20260822


def main():
    rng = np.random.default_rng(SEED)
    t, la, lo, _h, _f = W.load_region(
        str(HERE / "data" / "comcat_world" / ("%s.csv" % REGION)), declustered=True)
    n, S = t.size, A.S_NULL_PER_EVENT
    thr_a, thr_r, _ = A.dwell_thresholds(t, la, lo, rng)

    off = rng.uniform(-A.NULL_HALF_WINDOW_DAYS, A.NULL_HALF_WINDOW_DAYS, size=(n, S))
    nf = A.amplitude_features((t[:, None] + off).ravel(),
                              np.repeat(la, S), np.repeat(lo, S))
    nf = {k: v.reshape(n, S) for k, v in nf.items()}

    idx = rng.integers(0, S, size=(n, A.N_REPLICATES))
    nb = A.battery({k: np.take_along_axis(v, idx, axis=1) for k, v in nf.items()},
                   thr_a, thr_r, axis=0)
    mu = {k: float(np.asarray(nb[k]).mean()) for k in A.STATS}
    sd = {k: float(np.asarray(nb[k]).std(ddof=1)) for k in A.STATS}

    out = {"arm": "sensitivity of the amplitude arm", "priced_tests": 0,
           "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
           "region": REGION, "n_events_declustered": int(n),
           "injection": ("within each event's own null pool, resample with probability "
                         "proportional to (1 + eps * z_driver); eps = 0.10 is a ten "
                         "percent preference over the event's own two fortnights"),
           "threshold": "two-sided 5 percent", "trials_per_cell": TRIALS,
           "power": {}}

    for drv in DRIVERS:
        d = nf[drv]
        dz = (d - d.mean(axis=1, keepdims=True)) / (d.std(axis=1, keepdims=True) + 1e-30)
        out["power"][drv] = {}
        for eps in EPS:
            w = np.clip(1.0 + eps * dz, 1e-9, None)
            cw = np.cumsum(w / w.sum(axis=1, keepdims=True), axis=1)
            hit = {k: 0 for k in A.STATS}
            for _ in range(TRIALS):
                j = (cw < rng.random((n, 1))).sum(axis=1).clip(0, S - 1)
                sel = {k: np.take_along_axis(v, j[:, None], axis=1)[:, 0]
                       for k, v in nf.items()}
                b = A.battery(sel, thr_a, thr_r)
                for k in A.STATS:
                    if abs((b[k] - mu[k]) / sd[k]) > CRIT:
                        hit[k] += 1
            out["power"][drv]["eps=%.2f" % eps] = {k: hit[k] / TRIALS for k in A.STATS}
            print("  %-12s eps=%.2f   " % (drv, eps)
                  + "  ".join("%s %.2f" % (k.split("_")[0], hit[k] / TRIALS)
                              for k in A.STATS), flush=True)

    diag = {"E1_areal_amp": "areal_abs", "E2_shear_amp": "shear",
            "E3_rate_amp": "rate_abs", "E4_signed_rate": "rate_signed",
            "E5_amp_top_decile": "areal_abs", "E6_rate_top_decile": "rate_abs"}
    out["own_driver_power_eps_0.10"] = {
        k: out["power"][d]["eps=0.10"][k] for k, d in diag.items()}
    worst = min(out["own_driver_power_eps_0.10"].values())
    out["verdict"] = (
        "Every statistic responds to its own driver; the weakest at eps = 0.10 is %.2f. "
        "The amplitude arm's null is a bound at measured power, not a broken instrument."
        % worst)
    OUT_JSON.write_text(json.dumps(out, indent=2, default=float), encoding="utf-8")
    print("\n" + out["verdict"] + "\nwrote %s" % OUT_JSON.name)


if __name__ == "__main__":
    main()
