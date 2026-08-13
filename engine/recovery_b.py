"""S-17 DEMONSTRATED RECOVERY for Tranche B's new statistics. SIMULATION ONLY.

§P7-3 attaches a recovery demand to every one of B's new estimators, and §P7-15(b)
restates the list so nothing is inferred:

  * **F9-01 second moment** -- two-lobed positive control AND **the sinusoid negative
    control** (the one Kepler did not name).
  * **F9-04 Kuiper/Watson** -- narrow-arc positive and day-lattice negative,
    demonstrated **band by band**.
  * **F9-10 mark path** -- **G-M1 arm (i) re-run on the mark path specifically**,
    since it is not notched and arm (i) there becomes a live falsification rather
    than a formality.

G-M1's logic applies to ESTIMATORS, not only to pipelines (§P7-3): *"before any
result from it means anything."* This module is the harness that discharges those
demands. It runs on SIMULATED catalogues only -- B's run on real data is gated
(§P7-14(d)) and nothing here touches `data/`, `engine/EXPLORE_COUNT.jsonl`,
`engine/out/mine` or any holdout hash.

WHY THE PLANTS ARE THE SIZE THEY ARE
------------------------------------
§P7-8(d), enforced by `engine/floors.py`: *plant at >= 2x the operative floor*, or a
failure to recover is a POWER verdict wearing an INSTRUMENT verdict's clothes. Every
plant below goes through `floors.assert_plant_above_floor` (count path, rate-
modulation units) or `floors.assert_mark_plant_above_floor` (mark path, correlation
units, §P7-10(c)) and the harness REFUSES to run a non-compliant plant.

The floor used is the declared Tranche B floor at `alpha = 1.0e-4` under the measured
VIF = 24.08, and the simulator is built to MATCH that VIF rather than to be easier
than it (see `simulate_phase_catalog`). That matters more than it sounds: plants sized
for a VIF-24 world and then evaluated in a VIF-1 world would recover trivially, and
every recovery rate in this table would be an artifact of the simulator being kind.

AND THE UNITS THE PLANTS ARE SIZED IN, because getting this wrong is the exact
category error §P7-8(d) exists to catch. The floor is quoted in SINUSOIDAL amplitude.
A concentrated arc is not a sinusoid, so it is converted first:
`circstat.kuiper_equivalent_amplitude` (CDF-excursion units -> equivalent sinusoid)
for the omnibus arm, and `floors.rho_min` (correlation units, §P7-10(c)) for the mark
arm. Sizing an arc by its raw height against a sinusoidal floor would compare a CDF
excursion to a Fourier coefficient.

THE SIMULATOR CARRIES DISPERSION ON PURPOSE
-------------------------------------------
`simulate_phase_catalog` puts a slow AR(1) log-rate field into the intensity and does
NOT give it to the baseline offset, plus a broadband gamma mixture that sets the VIF
at every frequency. That is not a nicety: a clean Poisson world has a null so narrow
that every statistic fires on everything, and a "Kuiper found what Rayleigh missed"
demonstration on a Poisson world would be a demonstration about Poisson. The excess
low-frequency residual power is also exactly what §P7-10(d) named as the leading
explanation of the real VIF, so the simulated world fails in the direction the real
one does.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os

import numpy as np

from . import circstat, floors, marks_ext, mine as M

RECOVERY_ID = "TRANCHE-B-S17-RECOVERY"
RESULTS_JSON = "results_recovery_b.json"
OUT_ROOT = os.path.join("engine", "out", "recovery_b")

SIM_ONLY_BANNER = (
    "SIMULATION ONLY. This harness runs on synthetic catalogues generated in this "
    "process. It reads no real catalogue, writes no exploration-ledger line, spends "
    "no holdout hash and makes no declaration. Its product is a RECOVERY TABLE -- an "
    "instrument calibration -- and it may not be entered for or against any ledger "
    "entry (§P7-14(d): B's BUILD is authorized, B's RUN is separately gated).")

# Declared bands for the band-by-band demonstration §P7-3(2) requires. Fortnightly is
# where v1 actually claims; monthly and annual bracket it; the diurnal band is
# ABSENT ON PURPOSE -- it is structurally notched on the count path and demonstrating
# there would be demonstrating against a zero (§P5-6's band-matching).
DECLARED_BANDS = (("fortnightly", 14.765), ("monthly", 29.530588),
                  ("annual", 365.2422))

# Frozen before the run (S-9). One declared value each.
ALPHA = floors.ALPHA_TRANCHE_B                    # 1.0e-4
DETECT_ALPHA = 0.01          # the declared detection threshold for a recovery count
N_SURR = 1000
N_REPS = 10
NARROW_DUTY = 0.10
N_ARCS = 3                   # 3 equally spaced arcs kill harmonics 1 AND 2 exactly
DUTY_EACH = 0.05             # duty of EACH arc in the fundamental-suppressed plant
RECOVERY_BAR = 0.80          # G-M1's own bar: >= 80% recovery

# THE SIMULATED WORLD IS SIZED TO THE REAL ONE, and that is a correctness
# requirement rather than realism for its own sake. A recovery demonstration is a
# statement about THIS instrument at THIS operating point: the floor scales as
# 1/sqrt(N), and the block-bootstrap null's width at a band scales with how many
# BLOCKS of `2 x period` fit in the record. Demonstrate the annual band on a 3,000 d
# record and the bootstrap has four blocks to work with and the band fails -- which
# is a true statement about a 3,000 d record and a false one about the 7,716 d window
# Tranche B would actually run on. So: 7,716 days at 6.04 events/day = 46,600 events,
# which is the v2 exploration window and N = 46,585 to three significant figures.
SIM_N_DAYS = 7716
SIM_RATE = 6.04
GM1_ARM_I_COUNT_PATH_BAR = 0.10   # arm (i)'s v1 requirement: A_hat / A_ref < 0.1
GM1_ARM_II_RANGE = (0.8, 1.2)     # the amplitude-recovery tolerance


# ============================================================== simulators ====
SIM_TARGET_VIF = floors.MEASURED_VIF_DF2_PHASE     # 24.08, the measured F4-58 value


def simulate_phase_catalog(n_days, period_days, modulation=None, seed=0,
                           rate=6.0, od_phi=0.995, od_sd=0.45,
                           target_vif=SIM_TARGET_VIF):
    """A day-binned catalogue with a planted phase modulation and real dispersion.

    Returns (theta, counts, offset, meta). `offset` is the baseline the miner would
    have: the MEAN rate, with no knowledge of the AR(1) log-rate field -- so the
    residual carries excess low-frequency power and the block-bootstrap null has the
    width a real null has.

    TWO DISPERSION TERMS, AND THE SECOND ONE IS NOT COSMETIC
    -------------------------------------------------------
    The AR(1) log-rate field is LOW-FREQUENCY: at phi = 0.995 its correlation time is
    ~200 d, so it inflates the variance of a decadal statistic and leaves a
    fortnightly one almost untouched. A world with only that term has an effective
    VIF near 1 at the bands this tranche tests, and every statistic fires on
    everything -- which would make "Kuiper found what the 2-df form missed" a
    statement about Poisson data rather than about shape.

    So a second, BROADBAND term is added: the daily rate is gamma-mixed
    (negative-binomial counts) with shape `k = rate / (target_vif - 1)`, which gives
    `Var = mu * target_vif` at EVERY frequency. `target_vif` defaults to the measured
    F4-58 value 24.08 -- the same number `engine/floors.py` builds the declared floor
    from -- so the simulated world's null width and the floor the plants are sized
    against are the SAME arithmetic rather than two unrelated ones. Without this the
    plant discipline would be internally inconsistent: plants sized for a VIF-24
    world, evaluated in a VIF-1 world.
    """
    rng = np.random.default_rng(int(seed))
    t = np.arange(int(n_days), dtype=np.float64)
    th = np.mod(2 * np.pi * t / float(period_days), 2 * np.pi)
    u = np.zeros(int(n_days))
    for i in range(1, int(n_days)):
        u[i] = od_phi * u[i - 1] + np.sqrt(1 - od_phi ** 2) * rng.normal()
    u = od_sd * u
    m = np.ones(int(n_days)) if modulation is None else np.asarray(
        modulation, dtype=np.float64)
    lam = float(rate) * np.exp(u - 0.5 * od_sd ** 2) * m
    if target_vif and float(target_vif) > 1.0:
        k = float(rate) / (float(target_vif) - 1.0)
        lam = lam * rng.gamma(k, 1.0 / k, size=int(n_days))
    counts = rng.poisson(lam).astype(np.float64)
    offset = np.full(int(n_days), float(rate))
    return th, counts, offset, {
        "n_days": int(n_days), "period_days": float(period_days),
        "n_events": float(counts.sum()), "rate": float(rate),
        "od_phi": float(od_phi), "od_sd": float(od_sd),
        "target_vif": (float(target_vif) if target_vif else 1.0),
        "dispersion_note": ("the AR(1) log-rate field is IN the intensity and NOT "
                            "in the offset, so the residual carries excess "
                            "low-frequency power -- §P7-10(d)'s own candidate "
                            "explanation of the measured VIF -- and a broadband "
                            "gamma mixture supplies the dispersion the declared "
                            "floor assumes at the bands actually tested"),
    }


def sinusoid_modulation(theta, amplitude, phase=0.0):
    """`1 + A cos(theta - phi)` -- the pure first-harmonic alternative."""
    return 1.0 + float(amplitude) * np.cos(np.asarray(theta) - float(phase))


def antipodal_modulation(theta, amplitude, phase=0.0):
    """`1 + A cos(2(theta - phi))` -- two lobes 180 deg apart, FIRST MOMENT ZERO.

    §P7-3(1)'s positive control: *"a planted antipodal two-lobed phase distribution
    with first moment ~ 0 -- the statistic must recover it."* This is K-088's shape:
    unlock at phase X, release at X + pi, and the first moment cancels exactly while
    the second survives at full strength.
    """
    return 1.0 + float(amplitude) * np.cos(2.0 * (np.asarray(theta) - float(phase)))


# ============================== F9-01: positive control AND the sinusoid one ===
def _second_moment_and_rayleigh(th, counts, offset, period_days, n_surr, seed):
    """Both statistics on the same data under the same nulls -- the only fair form."""
    rng = np.random.default_rng(int(seed))
    block = float(np.clip(2.0 * period_days, 30.0, 800.0))
    m2 = circstat.second_moment_test(th, counts, offset, n_surr, rng,
                                     block_days=block, periodic=True)
    rng2 = np.random.default_rng(int(seed) + 1)
    X1 = np.column_stack([np.sin(th), np.cos(th)])
    S = M.score_stat_all_shifts(X1, counts, offset)
    Sb = M.score_stat_block_bootstrap(X1, counts, offset, int(n_surr), rng2,
                                      mean_block=block)
    p1 = M.bootstrap_p(S[0], Sb)
    return m2, {"test": "rayleigh_form_2df_score", "statistic": float(S[0]),
                "p_block_bootstrap": float(p1), "p_raw": float(p1)}


def f9_01_recovery(n_days=SIM_N_DAYS, rate=SIM_RATE, n_reps=N_REPS, n_surr=N_SURR,
                   alpha=DETECT_ALPHA, bands=DECLARED_BANDS, seed0=901,
                   verbose=True):
    """F9-01's two S-17 controls, band by band. Returns the recovery table.

    POSITIVE: an antipodal two-lobed plant (first moment ~ 0). The second-moment
              statistic must recover it.
    NEGATIVE: a pure sinusoid at the SAME amplitude. The second-moment statistic must
              NOT fire beyond its nominal size -- §P7-3(1), the control without which
              a second-moment detection is unfalsifiably confounded with the first
              moment this programme has already bounded.
    """
    out = []
    for band, P in bands:
        n_ev = rate * n_days
        amp = floors.PLANT_FACTOR * floors.a_min(floors.MEASURED_VIF_DF2_PHASE,
                                                 ALPHA, n_ev)
        plant = floors.assert_plant_above_floor(
            amp, n_ev, ALPHA, what="F9-01 %s band" % band)
        pos_hits = neg_hits = 0
        pos_rows, neg_rows = [], []
        for r in range(int(n_reps)):
            s = int(seed0) + 1000 * len(out) + r
            th = np.mod(2 * np.pi * np.arange(n_days) / P, 2 * np.pi)
            # POSITIVE: antipodal
            _th, c, o, meta = simulate_phase_catalog(
                n_days, P, antipodal_modulation(th, amp), seed=s, rate=rate)
            m2, r1 = _second_moment_and_rayleigh(th, c, o, P, n_surr, s)
            pos_hits += int(m2["p_raw"] <= alpha)
            pos_rows.append({"rep": r, "p_second_moment": m2["p_raw"],
                             "p_rayleigh_form": r1["p_raw"],
                             "R1": m2["R1"], "R2": m2["R2"],
                             "reading": m2["reading"]})
            # NEGATIVE: pure sinusoid, SAME amplitude
            _th, c2, o2, _m = simulate_phase_catalog(
                n_days, P, sinusoid_modulation(th, amp), seed=s + 500000, rate=rate)
            m2n, r1n = _second_moment_and_rayleigh(th, c2, o2, P, n_surr,
                                                   s + 500000)
            neg_hits += int(m2n["p_raw"] <= alpha)
            neg_rows.append({"rep": r, "p_second_moment": m2n["p_raw"],
                             "p_rayleigh_form": r1n["p_raw"],
                             "R1": m2n["R1"], "R2": m2n["R2"]})
        row = {
            "band": band, "period_days": P,
            "planted_amplitude": amp, "plant_report": plant,
            "n_reps": int(n_reps), "detect_alpha": float(alpha),
            "positive_two_lobed_recovery_rate": pos_hits / float(n_reps),
            "negative_sinusoid_false_positive_rate": neg_hits / float(n_reps),
            "nominal_size": float(alpha),
            "positive_rows": pos_rows, "negative_rows": neg_rows,
            "positive_pass": pos_hits / float(n_reps) >= RECOVERY_BAR,
            "negative_pass": neg_hits / float(n_reps) <= max(5.0 * alpha, 0.10),
            "rayleigh_on_two_lobed_median_p": float(np.median(
                [x["p_rayleigh_form"] for x in pos_rows])),
            "why_the_negative_control_is_the_one_that_matters": (
                "§P7-3(1): without the SINUSOID negative control a 'second-moment "
                "detection' is unfalsifiably confounded with the first moment the "
                "programme has already bounded. The positive control alone cannot "
                "distinguish a new estimator from a re-labelled old one."),
        }
        if verbose:
            print("  F9-01 %-12s P=%8.3f d  plant %.3f  two-lobed recovery %.0f%%  "
                  "sinusoid false-positive %.0f%% (nominal %.0f%%)"
                  % (band, P, amp, 100 * row["positive_two_lobed_recovery_rate"],
                     100 * row["negative_sinusoid_false_positive_rate"],
                     100 * alpha))
        out.append(row)
    return {"statistic": "F9-01 second circular moment", "rows": out,
            "bar": RECOVERY_BAR,
            "pass": all(r["positive_pass"] and r["negative_pass"] for r in out)}


# ================== F9-04: narrow-arc positive, day-lattice negative, and THE ONE ==
def _omnibus_and_rayleigh(th, counts, offset, period_days, n_surr, seed):
    rng = np.random.default_rng(int(seed))
    block = float(np.clip(2.0 * period_days, 30.0, 800.0))
    rows = circstat.omnibus_test(th, counts, offset, n_surr, rng,
                                 block_days=block, periodic=True)
    rng2 = np.random.default_rng(int(seed) + 1)
    X1 = np.column_stack([np.sin(th), np.cos(th)])
    S = M.score_stat_all_shifts(X1, counts, offset)
    Sb = M.score_stat_block_bootstrap(X1, counts, offset, int(n_surr), rng2,
                                      mean_block=block)
    return rows, float(M.bootstrap_p(S[0], Sb))


def f9_04_recovery(n_days=SIM_N_DAYS, rate=SIM_RATE, n_reps=N_REPS, n_surr=N_SURR,
                   alpha=DETECT_ALPHA, bands=DECLARED_BANDS, duty=NARROW_DUTY,
                   n_arcs=N_ARCS, duty_each=DUTY_EACH, seed0=904, verbose=True):
    """F9-04's controls band by band, and the money demonstration inside them.

    THREE ARMS, and the second is the one §P7-3(2) is really asking for.

    (A) SINGLE NARROW ARC, `duty` duty cycle -- the catalog's own positive control,
        sized so its KUIPER-equivalent amplitude is >= 2x the floor. Its rates are
        reported for BOTH statistics and they come out close, which is the honest
        measured answer: a single mean-preserving boxcar can buy Kuiper AT MOST a
        factor pi/2 = 1.571 in equivalent amplitude
        (`circstat.single_arc_kuiper_edge`, derived not asserted). MINING_CATALOG
        F9-04's "nearly invisible to it and obvious to Kuiper" OVERSTATES this case,
        and the harness says so with the number rather than reproducing the rhetoric.

    (B) THE DECISIVE PLANT: `n_arcs` equally spaced arcs of `duty_each` duty each. At
        n_arcs = 3 the first AND second harmonics vanish EXACTLY, so the 2-df
        Rayleigh-form statistic and the F9-01 second-moment statistic are blind BY
        CONSTRUCTION while Kuiper and Watson read the whole CDF. This is
        "concentrated-phase structure recovered where the Rayleigh-form test misses
        it", with the miss guaranteed by algebra stated in advance instead of by a
        tuned amplitude -- and the F9-01 column in the same table shows the two new
        statistics are COMPLEMENTARY, not two names for one repair.

    (C) NEGATIVE: the day-binning lattice itself. A feature whose period divides the
        1-day lattice has ONE distinct phase, so V and U^2 are identically zero and
        neither statistic may fire. §P7-3(2)'s negative control, and the reason
        `circstat.phase_group_ends` exists.
    """
    out = []
    for band, P in bands:
        n_ev = rate * n_days
        fl = floors.a_min(floors.MEASURED_VIF_DF2_PHASE, ALPHA, n_ev)
        a_eq = floors.PLANT_FACTOR * fl        # the plant, in KUIPER-equivalent units

        amp_single = circstat.arc_amplitude_for_kuiper_equivalent(a_eq, duty)
        fund_single = circstat.fundamental_coefficient(amp_single, duty)
        plant_single = floors.assert_plant_above_floor(
            a_eq, n_ev, ALPHA,
            what="F9-04 %s single narrow arc (Kuiper-equivalent)" % band)
        amp_kfold = circstat.arc_amplitude_for_kuiper_equivalent(a_eq, duty_each)
        plant_kfold = floors.assert_plant_above_floor(
            a_eq, n_ev, ALPHA,
            what="F9-04 %s %d-fold arc (Kuiper-equivalent)" % (band, n_arcs))

        th = np.mod(2 * np.pi * np.arange(n_days) / P, 2 * np.pi)
        single = {"k": 0, "w": 0, "r": 0, "rows": []}
        kfold = {"k": 0, "w": 0, "r": 0, "m2": 0, "rows": []}
        for r in range(int(n_reps)):
            s = int(seed0) + 1000 * len(out) + r
            _t, c, o, _m = simulate_phase_catalog(
                n_days, P, circstat.narrow_arc_intensity(th, amp_single, duty=duty),
                seed=s, rate=rate)
            omni, p_ray = _omnibus_and_rayleigh(th, c, o, P, n_surr, s)
            pk, pw = float(omni[0]["p_raw"]), float(omni[1]["p_raw"])
            single["k"] += int(pk <= alpha)
            single["w"] += int(pw <= alpha)
            single["r"] += int(p_ray <= alpha)
            single["rows"].append({"rep": r, "p_kuiper": pk, "p_watson": pw,
                                   "p_rayleigh_form": p_ray})

            s2 = s + 40000
            _t, c2, o2, _m2 = simulate_phase_catalog(
                n_days, P,
                circstat.k_fold_arc_intensity(th, amp_kfold, n_arcs, duty_each),
                seed=s2, rate=rate)
            omni2, p_ray2 = _omnibus_and_rayleigh(th, c2, o2, P, n_surr, s2)
            m2row, _r1 = _second_moment_and_rayleigh(th, c2, o2, P, n_surr, s2)
            pk2, pw2 = float(omni2[0]["p_raw"]), float(omni2[1]["p_raw"])
            kfold["k"] += int(pk2 <= alpha)
            kfold["w"] += int(pw2 <= alpha)
            kfold["r"] += int(p_ray2 <= alpha)
            kfold["m2"] += int(float(m2row["p_raw"]) <= alpha)
            kfold["rows"].append({"rep": r, "p_kuiper": pk2, "p_watson": pw2,
                                  "p_rayleigh_form": p_ray2,
                                  "p_second_moment": float(m2row["p_raw"])})

        neg_hits, neg_rows = 0, []
        for r in range(int(n_reps)):
            s = int(seed0) + 90000 + 1000 * len(out) + r
            th_lat = circstat.day_lattice_phase(n_days, period_days=1.0)
            _t, c, o, _m = simulate_phase_catalog(n_days, 1.0, None, seed=s,
                                                  rate=rate)
            om = circstat.omnibus_test(th_lat, c, o, n_surr,
                                       np.random.default_rng(s),
                                       block_days=30.0, periodic=True)
            p_lat = float(min(om[0]["p_raw"], om[1]["p_raw"]))
            neg_hits += int(p_lat <= alpha)
            neg_rows.append({"rep": r, "p_min_omnibus": p_lat,
                             "V_star": float(om[0]["statistic"]),
                             "U2_star": float(om[1]["statistic"]),
                             "n_distinct_phases": int(
                                 circstat.phase_group_ends(np.sort(th_lat)).size)})

        nr = float(n_reps)
        row = {
            "band": band, "period_days": P,
            "operative_floor_A_min": fl,
            "planted_kuiper_equivalent_amplitude": a_eq,
            "plant_factor_over_floor": a_eq / fl,
            "n_reps": int(n_reps), "detect_alpha": float(alpha),

            "A_single_arc": {
                "duty_cycle": float(duty),
                "arc_height": amp_single,
                "fundamental_presented_to_2df": fund_single,
                "fundamental_over_floor": fund_single / fl,
                "kuiper_edge_ratio_ceiling": circstat.single_arc_kuiper_edge(duty),
                "kuiper_recovery_rate": single["k"] / nr,
                "watson_recovery_rate": single["w"] / nr,
                "rayleigh_recovery_rate": single["r"] / nr,
                "reading": ("a single mean-preserving boxcar buys Kuiper at most a "
                            "factor %.3f in equivalent amplitude (derived), so BOTH "
                            "statistics are expected to see a plant sized at 2x the "
                            "floor in Kuiper-equivalent units. MINING_CATALOG "
                            "F9-04's 'nearly invisible / obvious' overstates the "
                            "single-arc case; arm B is the decisive one."
                            % circstat.single_arc_kuiper_edge(duty)),
                "rows": single["rows"],
            },

            "B_fundamental_suppressed": {
                "n_arcs": int(n_arcs), "duty_each": float(duty_each),
                "arc_height": amp_kfold,
                "harmonics_suppressed": [m for m in (1, 2) if m % int(n_arcs)],
                "kuiper_recovery_rate": kfold["k"] / nr,
                "watson_recovery_rate": kfold["w"] / nr,
                "rayleigh_recovery_rate": kfold["r"] / nr,
                "second_moment_recovery_rate": kfold["m2"] / nr,
                "kuiper_sees_and_rayleigh_misses_rate": (
                    sum(1 for x in kfold["rows"]
                        if min(x["p_kuiper"], x["p_watson"]) <= alpha
                        and x["p_rayleigh_form"] > alpha) / nr),
                "reading": ("%d equally spaced arcs put ZERO power in harmonics 1 "
                            "and 2, so the 2-df Rayleigh-form statistic AND the "
                            "F9-01 second moment are blind by algebra, not by luck. "
                            "This is the money demonstration and it also shows the "
                            "two new statistics are complementary." % n_arcs),
                "rows": kfold["rows"],
            },

            "C_day_lattice_negative": {
                "false_positive_rate": neg_hits / nr,
                "nominal_size": float(alpha),
                "rows": neg_rows,
            },

            "plant_reports": {"single_arc": plant_single, "k_fold": plant_kfold},
            "positive_pass": bool(kfold["k"] / nr >= RECOVERY_BAR),
            "decisive_pass": bool(kfold["k"] / nr >= RECOVERY_BAR
                                  and kfold["r"] / nr <= max(5.0 * alpha, 0.10)),
            "negative_pass": bool(neg_hits / nr <= max(5.0 * alpha, 0.10)),
            "band_note": circstat.BAND_NOTE,
        }
        if verbose:
            print("  F9-04 %-12s P=%8.3f d  a_eq %.3f (2x floor %.3f)"
                  % (band, P, a_eq, fl))
            print("        (A) single arc  d=%.2f: Kuiper %.0f%%  Watson %.0f%%  "
                  "Rayleigh-form %.0f%%   [edge ceiling x%.3f]"
                  % (duty, 100 * single["k"] / nr, 100 * single["w"] / nr,
                     100 * single["r"] / nr,
                     circstat.single_arc_kuiper_edge(duty)))
            print("        (B) %d-fold arc  : Kuiper %.0f%%  Watson %.0f%%  "
                  "Rayleigh-form %.0f%%  2nd-moment %.0f%%   <-- THE DEMONSTRATION"
                  % (n_arcs, 100 * kfold["k"] / nr, 100 * kfold["w"] / nr,
                     100 * kfold["r"] / nr, 100 * kfold["m2"] / nr))
            print("        (C) day lattice : false positives %.0f%% (nominal %.0f%%)"
                  % (100 * neg_hits / nr, 100 * alpha))
        out.append(row)
    return {"statistic": "F9-04 Kuiper V / Watson U^2", "rows": out,
            "bar": RECOVERY_BAR,
            "pass": all(r["decisive_pass"] and r["negative_pass"] for r in out)}


# ========================= F9-10: G-M1 ARM (i) RE-RUN ON THE MARK PATH ==========
def simulate_marked_catalog(n_days=1500, n_per_day=12.0, mag_amplitude=0.20,
                            seed=910, mc=4.5, beta=1.0):
    """A catalogue whose MAGNITUDES carry a local-solar-hour modulation.

    This is the S1/S2 detection artifact in the only form a mark test can see it: a
    time-of-day-dependent completeness threshold does not change the daily COUNT
    much (it is a global sum over all longitudes, which is §P5-1's second exact
    zero), but it does change WHICH events are in the catalogue -- and the mark axis
    sees exactly that (F9-10's own Hunch, restated as a plant).

    The planted quantity is the amplitude of the first harmonic of mean magnitude
    against LOCAL solar hour. It is recoverable in the mark's own units, which is
    what makes arm (i) a live falsification here rather than a formality.
    """
    rng = np.random.default_rng(int(seed))
    n = int(n_days * n_per_day)
    day_float = np.sort(rng.uniform(0.0, float(n_days), size=n))
    lon = rng.uniform(-180.0, 180.0, size=n)
    lat = rng.uniform(-60.0, 60.0, size=n)
    th_loc = marks_ext.local_solar_hour_phase(day_float, lon)
    mag = mc + rng.exponential(1.0 / (beta * np.log(10.0)), size=n)
    mag = mag + float(mag_amplitude) * np.cos(th_loc)
    depth = rng.uniform(1.0, 200.0, size=n)
    marks = {"day": np.floor(day_float).astype(np.int64), "day_float": day_float,
             "mag": mag, "depth": depth, "lat": lat, "lon": lon}
    return marks, {"n_events": int(n), "n_days": int(n_days),
                   "planted_mag_amplitude": float(mag_amplitude),
                   "mc": float(mc), "beta": float(beta)}


def gm1_arm_i_mark_path(n_days=1500, n_per_day=12.0, mag_amplitude=0.20,
                        n_reps=5, n_surr=N_SURR, seed0=910, verbose=True):
    """G-M1 arm (i), re-run ON THE MARK PATH -- §P7-3(3)'s live falsification test.

    On the COUNT path arm (i) demanded `A_hat / A_ref < 0.1` and could only destroy:
    v1 is *designed* blind there by two exact zeros. On the MARK path there is no
    zero, so the demand INVERTS: the mark path must RECOVER the planted local-solar-
    hour modulation to G-M1 arm (ii)'s own tolerance [0.8, 1.2]. A mark path that
    cannot recover a known artifact of known amplitude has no business reporting a
    sub-daily mark result -- and one that recovers it has just measured its own
    sub-daily sensitivity, which is the number the whole sub-daily arm rests on.

    The count path is run on the SAME catalogues in the same call, so the pair of
    numbers -- mark recovery ~1, count recovery ~0 -- is the escape from the sinc
    MEASURED rather than asserted.
    """
    rows = []
    ok = 0
    for r in range(int(n_reps)):
        s = int(seed0) + r
        marks, meta = simulate_marked_catalog(n_days, n_per_day, mag_amplitude,
                                              seed=s)
        th_loc = marks_ext.local_solar_hour_phase(marks["day_float"], marks["lon"])
        a_hat, phi_hat = marks_ext.harmonic_amplitude(th_loc, marks["mag"])
        ratio = a_hat / float(mag_amplitude)

        # the mark test itself, under the engine's unmodified statistic
        rng = np.random.default_rng(s + 7)
        mt = M.mark_test(th_loc, marks["mag"], "phase", int(n_surr), rng)
        n_ev = int(marks["mag"].size)
        mfl = floors.mark_floor_report(n_ev, None, feature="local_solar_hour x mag")

        # THE COUNT PATH ON THE SAME CATALOGUE. A daily bin has exactly one value of
        # the diurnal phase, so the design column is constant: the recovery is zero
        # BY CONSTRUCTION and the engine refuses to build the design at all. That
        # refusal IS §P5-1's first exact zero, caught rather than described.
        counts = np.bincount(marks["day"], minlength=int(n_days)).astype(float)
        th_day = circstat.day_lattice_phase(int(n_days), period_days=1.0,
                                            phase0=np.pi)
        design_sd = float(np.std(np.cos(th_day)))
        if design_sd < 1e-12:
            # THE EXACT ZERO, caught rather than described. The design column is
            # CONSTANT, so the count path cannot form the test at all: the recovered
            # amplitude is 0 identically, not "small". A least-squares fit on a
            # rank-deficient design would still RETURN a number -- and that number is
            # meaningless noise from the pseudo-inverse, which is exactly how a
            # structural zero gets accidentally reported as a 2.5x over-recovery.
            count_ratio = 0.0
            count_note = ("STRUCTURAL ZERO (§P5-1's first exact zero): a daily bin "
                          "has exactly one value of the diurnal phase, so the "
                          "design column is constant (sd = %.3g) and the count path "
                          "cannot form the test. A_hat/A_ref = 0 identically, not "
                          "by measurement." % design_sd)
        else:
            a_c, _ = marks_ext.harmonic_amplitude(
                th_day, counts / max(counts.mean(), 1e-12))
            count_ratio = float(a_c) / float(mag_amplitude)
            count_note = "measured on a non-degenerate design"

        in_range = GM1_ARM_II_RANGE[0] <= ratio <= GM1_ARM_II_RANGE[1]
        ok += int(in_range)
        rows.append({
            "rep": r, "n_events": n_ev,
            "planted_mag_amplitude": float(mag_amplitude),
            "recovered_mag_amplitude": float(a_hat),
            "recovery_ratio_mark_path": float(ratio),
            "recovered_phase_rad": float(phi_hat),
            "phase_error_deg": float(np.degrees(np.angle(np.exp(1j * phi_hat)))),
            "in_gm1_range": bool(in_range),
            "mark_statistic": float(mt["statistic"]),
            "mark_p_raw": float(mt["p_raw"]),
            "rho_min": mfl["rho_min"],
            "mark_statistic_over_floor": float(mt["statistic"]) / mfl["rho_min"],
            "count_path_recovery_ratio": count_ratio,
            "count_path_design_sd_of_cos_phase": design_sd,
            "count_path_note": count_note,
        })
        if verbose:
            print("  G-M1 arm (i) mark path rep %d: A_hat/A = %.3f (mark), "
                  "%.3g (count path), mark p = %.3g, stat/floor = %.1f"
                  % (r, ratio, count_ratio, mt["p_raw"],
                     rows[-1]["mark_statistic_over_floor"]))
    rate = ok / float(n_reps)
    count_max = max(abs(x["count_path_recovery_ratio"]) for x in rows)
    return {
        "statistic": "F9-10 mark path -- G-M1 arm (i) re-run (§P7-3(3))",
        "n_reps": int(n_reps),
        "mark_path_recovery_rate_in_[0.8,1.2]": rate,
        "mark_path_median_ratio": float(np.median(
            [x["recovery_ratio_mark_path"] for x in rows])),
        "count_path_max_abs_ratio": float(count_max),
        "count_path_bar": GM1_ARM_I_COUNT_PATH_BAR,
        "count_path_structurally_blind": bool(count_max
                                              < GM1_ARM_I_COUNT_PATH_BAR),
        "bar": RECOVERY_BAR,
        "pass": bool(rate >= RECOVERY_BAR
                     and count_max < GM1_ARM_I_COUNT_PATH_BAR),
        "rows": rows,
        "inversion_note": (
            "arm (i) INVERTS on the mark path. On the count path it demanded "
            "A_hat/A_ref < 0.1 and could only destroy (v1 is designed blind there). "
            "On the mark path there is no notch, so the demand is RECOVERY to "
            "[0.8, 1.2] -- and the count-path column in this same table is the "
            "notch, measured on the identical catalogues."),
        "what_it_licenses": (
            "the mark path's sub-daily SENSITIVITY, measured. It does NOT license "
            "any sub-daily mark result on real data to be read as Earth: §P7-3(3) "
            "gates that arm on the F7-01/02/03 observer controls, and this "
            "harness's own plant is an OBSERVER artifact by construction -- which "
            "is precisely why the gate exists."),
    }


# ================================================================== driver ====
def run_all(quick=False, verbose=True):
    n_reps = 3 if quick else N_REPS
    n_surr = 200 if quick else N_SURR
    n_days = 1200 if quick else SIM_N_DAYS
    bands = DECLARED_BANDS[:1] if quick else DECLARED_BANDS
    t_open = _dt.datetime.now(_dt.timezone.utc).isoformat()
    if verbose:
        print("=" * 78)
        print("TRANCHE B -- S-17 DEMONSTRATED RECOVERY (%s)"
              % ("QUICK" if quick else "FULL"))
        print("=" * 78)
        print(SIM_ONLY_BANNER)
        print("")
    rate = 8.0 if quick else SIM_RATE
    f901 = f9_01_recovery(n_days=n_days, rate=rate, n_reps=n_reps, n_surr=n_surr,
                          bands=bands, verbose=verbose)
    f904 = f9_04_recovery(n_days=n_days, rate=rate, n_reps=n_reps, n_surr=n_surr,
                          bands=bands, verbose=verbose)
    f910 = gm1_arm_i_mark_path(n_reps=max(3, n_reps // 2), n_surr=n_surr,
                               verbose=verbose)
    payload = {
        "id": RECOVERY_ID,
        "banner": SIM_ONLY_BANNER,
        "opened_utc": t_open,
        "closed_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(),
        "quick": bool(quick),
        "declared": {"alpha_tranche_b": ALPHA, "detect_alpha": DETECT_ALPHA,
                     "n_reps": n_reps, "n_surrogates": n_surr,
                     "sim_n_days": n_days, "sim_rate": rate,
                     "recovery_bar": RECOVERY_BAR, "bands": list(bands),
                     "narrow_arc_duty": NARROW_DUTY,
                     "plant_factor": floors.PLANT_FACTOR,
                     "vif_count_path": floors.MEASURED_VIF_DF2_PHASE,
                     "vif_mark_fallback": floors.VIF_MARK_FALLBACK},
        "F9-01": f901, "F9-04": f904, "F9-10": f910,
        "all_pass": bool(f901["pass"] and f904["pass"] and f910["pass"]),
        "what_none_of_this_licenses": (
            "B's RUN. §P7-14(d) authorized B's BUILD and gated its run separately; "
            "§P7-15(b) confirms the run is gated on per-statistic G-M1 clearance, "
            "which is what this table is EVIDENCE FOR and not a substitute for. "
            "Adjudication is the Popper seat's."),
    }
    return payload


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--json", default=RESULTS_JSON)
    a = ap.parse_args(argv)
    os.makedirs(OUT_ROOT, exist_ok=True)
    rep = run_all(quick=a.quick)
    with open(a.json, "w", encoding="utf-8") as fh:
        json.dump(rep, fh, indent=1, default=float)
    print("\nwrote %s   all_pass = %s" % (a.json, rep["all_pass"]))
    return rep


if __name__ == "__main__":
    main()
