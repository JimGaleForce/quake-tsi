"""F8-15 `clock_random_reference` -- THE NULL CLOCK, and the refusal that makes it binding.

MINING_CATALOG F8-15, in full because it is short and every clause is operative:

    Def: a random monotone time warp with the same coarse statistics as the real
         clocks.
    Src: closed form.
    Hunch: the NULL CLOCK. If six clocks are tried, the improvement of the best over
         a random warp is the only honest measure of whether the clock discovery is
         real.
    Pit: none; it is the control that makes the whole family interpretable, and
         without it a clock scan is a p-hacking machine.
    S15: not applicable.
    Price: 1 control per clock scan. **Mandatory if any clock scan is run.**

§P7-2(b) assigns it: *"F8-15 random-clock control | new | Build in A, USE in C.
Mandatory the moment any clock scan runs."* §P7-4 freezes its definition NOW, before
B's results are seen, because *"a clock family reshaped after seeing which statistics
survived is the most efficient forking-paths machine this program could build, and
freezing it costs nothing today."* That is why this module exists in the Tranche B
build even though its consumer is Tranche C: the definition has to be frozen while
nobody has an interest in its shape.

WHAT "THE SAME COARSE STATISTICS" MEANS HERE, made specific so it cannot drift
-----------------------------------------------------------------------------
A clock is a monotone map `t -> tau(t)` from calendar days to clock units. A real
clock (natural time, ETAS-rescaled time, aftershock-free time, completeness-corrected
time, reverse time) is a cumulative sum of NON-NEGATIVE, HIGHLY UNEVEN daily
increments: quiet stretches barely advance, a big sequence advances a year in a week.
A random warp that matched only the total duration would be a straw man -- it would
be beaten by any real clock for the trivial reason that real clocks are lumpy.

So `random_monotone_clock` matches, by construction:
  * the total elapsed clock time (endpoint identity),
  * the MEAN and the COEFFICIENT OF VARIATION of the daily increments, and
  * strict monotonicity (increments > 0, so the map is invertible and no two days
    collapse onto one clock instant).

It does NOT match the increments' autocorrelation, and that omission is DECLARED
rather than hidden: a real clock's lumpiness is temporally clustered (an aftershock
sequence is contiguous) and a shuffled-increment warp is not. `random_monotone_clock`
therefore comes in two declared modes -- `iid` (increments drawn independently) and
`block` (increments block-permuted from the reference clock, which preserves the
short-range clustering exactly). Declare which one before the scan; report both if
the answer depends on it, because "the answer depended on the control's construction"
is itself the result.
"""

from __future__ import annotations

import numpy as np

F8_15_RULE_ID = "F8-15-v1"
F8_15_MODES = ("iid", "block")
F8_15_DEFAULT_BLOCK_DAYS = 90

MANDATORY_NOTE = (
    "MINING_CATALOG F8-15: the random-clock control is MANDATORY if any clock scan "
    "is run -- 1 control per clock scan. Without it a clock scan is a p-hacking "
    "machine: the best of k warps beats the identity warp by construction, and the "
    "only interpretable quantity is the best REAL clock's margin over a RANDOM one.")


class RandomClockControlMissing(AssertionError):
    """A clock scan was declared without its mandatory F8-15 control."""


def clock_increments(tau):
    """Daily increments of a monotone clock, with the monotonicity checked not assumed."""
    tau = np.asarray(tau, dtype=np.float64)
    d = np.diff(tau)
    if d.size and (d <= 0).any():
        raise ValueError("clock is not strictly monotone: %d non-positive increments"
                         % int((d <= 0).sum()))
    return d


def clock_summary(tau):
    """The coarse statistics F8-15 requires a random warp to match."""
    d = clock_increments(tau)
    m = float(d.mean()) if d.size else 0.0
    sd = float(d.std(ddof=0)) if d.size else 0.0
    return {"n_days": int(np.asarray(tau).size), "total": float(tau[-1] - tau[0]),
            "mean_increment": m, "sd_increment": sd,
            "cv_increment": float(sd / m) if m > 0 else 0.0}


def random_monotone_clock(reference_tau, rng, mode="iid",
                          block_days=F8_15_DEFAULT_BLOCK_DAYS):
    """F8-15: a random monotone warp matching `reference_tau`'s coarse statistics.

    `iid`   -- gamma increments with the reference's mean and CV, drawn independently
               and then rescaled so the endpoint matches EXACTLY. Gamma because the
               increments must be strictly positive and right-skewed, which is what a
               clock built from event counts actually looks like; a normal draw would
               need truncation and the truncation would move the CV it was matching.
    `block` -- the reference's OWN increments, permuted in contiguous blocks of
               `block_days`. Matches mean, CV and short-range clustering exactly and
               destroys only the long-range alignment, which is the property a clock
               claim is about.
    """
    if mode not in F8_15_MODES:
        raise ValueError("F8-15 mode must be one of %r, got %r" % (F8_15_MODES, mode))
    ref = np.asarray(reference_tau, dtype=np.float64)
    d = clock_increments(ref)
    if d.size == 0:
        return ref.copy()
    if mode == "iid":
        m, cv = float(d.mean()), float(d.std(ddof=0) / max(d.mean(), 1e-300))
        shape = 1.0 / max(cv * cv, 1e-12)
        draw = rng.gamma(shape, m / shape, size=d.size)
    else:
        L = int(max(1, min(int(block_days), d.size)))
        nb = int(np.ceil(d.size / L))
        starts = rng.permutation(nb) * L
        draw = np.concatenate([d[s:s + L] for s in starts])[:d.size]
    draw = draw * (float(d.sum()) / float(draw.sum()))     # endpoint identity
    return np.concatenate([[float(ref[0])], float(ref[0]) + np.cumsum(draw)])


def random_clock_control(reference_tau, rng, mode="iid",
                         block_days=F8_15_DEFAULT_BLOCK_DAYS, label=""):
    """The control clock plus the audit line proving it matched what it claims to."""
    tau = random_monotone_clock(reference_tau, rng, mode, block_days)
    ref_s, new_s = clock_summary(reference_tau), clock_summary(tau)
    return {
        "rule_id": F8_15_RULE_ID, "mode": mode, "label": label,
        "block_days": (int(block_days) if mode == "block" else None),
        "tau": tau,
        "reference_summary": ref_s, "control_summary": new_s,
        "matched": {
            "total": abs(new_s["total"] - ref_s["total"]) < 1e-6 * max(
                abs(ref_s["total"]), 1.0),
            "mean_increment": abs(new_s["mean_increment"]
                                  - ref_s["mean_increment"]) < 1e-9 * max(
                abs(ref_s["mean_increment"]), 1.0),
            "cv_within_20pct": abs(new_s["cv_increment"] - ref_s["cv_increment"])
            <= 0.20 * max(ref_s["cv_increment"], 1e-12),
        },
        "declared_omission": (
            "the `iid` mode does NOT match the increments' autocorrelation; a real "
            "clock's lumpiness is temporally clustered and an iid warp's is not. Use "
            "`block` when that matters, declare which before the scan, and report "
            "both if the answer depends on it."),
        "mandatory_note": MANDATORY_NOTE,
    }


def assert_random_clock_control(clock_names, control_names):
    """Refuse a clock scan that has not declared its F8-15 control -- one per scan."""
    clocks = [c for c in clock_names]
    controls = [c for c in control_names]
    if clocks and len(controls) < 1:
        raise RandomClockControlMissing(
            "clock scan declares %d clock(s) %r and ZERO F8-15 random-clock "
            "controls. %s" % (len(clocks), clocks, MANDATORY_NOTE))
    return {"n_clocks": len(clocks), "n_controls": len(controls),
            "rule_id": F8_15_RULE_ID, "satisfied": True,
            "mandatory_note": MANDATORY_NOTE}
