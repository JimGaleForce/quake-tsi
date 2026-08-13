"""§P7-17 -- THE DISPOSITION TAXONOMY. Every executed row is one of four things.

WHAT THIS RESOLVES, AND WHY NONE OF THE OBVIOUS ANSWERS WAS RIGHT
------------------------------------------------------------------
The Tranche B build flagged a real hole: the engine executes the count-path GLM
sweep on every science feature in every session, and B's partition declares strata
only for its three new statistics. The flag posed three options -- suppress the GLM
axis, declare and price a GLM stratum, or leave it broken. §P7-17 took none of them,
because all three assume every executed row is the same KIND of thing. They are not:

  **COMPONENT-OF <parent>.** A row that exists only as part of another row's claim.
  The first-moment GLM on a cyclic feature is *inside* F9-01's second-moment claim --
  `R2` is uninterpretable without `R1` beside it (F9-01's own Pit clause), and the
  row reporting `R1` is not a second hypothesis, it is half of one. Likewise the
  Rayleigh-form statistic beside Kuiper in F9-04: **the comparison IS the claim.**
  A component owes no multiplicity because it makes no independent rejection -- and
  it must never be able to make one, which is what the enforcement below is for.

  **REPLICATION-OF-DECLARED -- CROSS-SESSION ONLY (§P7-18 clarifying §P7-17(2)).**
  A re-execution, IN A LATER SESSION, of a row already declared in an earlier one
  (Tranche A's 550, the period scan). Not re-priced -- paying twice for the same
  hypothesis inflates the denominator without adding a hypothesis. And it comes with
  a dividend §P7-17 makes mandatory rather than optional: **where the re-execution
  shares a test-key digest with the prior session, its p MUST be bitwise identical,
  and that is a free §P6-5 determinism check** on data we are running anyway.

  **A WITHIN-SESSION same-shape row is NOT a replication. It is a DUPLICATE, and it
  is SUPPRESSED (§P7-18).** The distinction is not bookkeeping. A cross-session
  replication is a second observation of the same hypothesis made under a different
  declaration, and its agreement with the first is INFORMATION (that is the whole
  determinism dividend). A within-session duplicate is the same statistic on the
  same data in the same table twice: it carries no information, and it reads as
  phantom replication -- one result appearing to be confirmed twice, which is the
  hazard §P7-16 de-duplicated `log_moment` for. Suppression is proved, not asserted:
  see `same_shape` and `suppression_map`.

  **DECLARED.** A genuinely new independent row. These, and only these, are `m`.

**THE TAXONOMY IS THREE TAGS (§P7-18).** COMPONENT-OF, REPLICATION-OF-DECLARED,
DECLARED -- and every EXECUTED, PRICEABLE row carries exactly one. Two further
classes exist and are deliberately NOT dispositions, because calling them
dispositions would blur what a disposition is for:

  **UNPRICED-CONTROL** (§P7-16) is a PRICING class, not a disposition: a control that
  can only calibrate a reference or condemn our own instrument owes no multiplicity,
  so the question "which kind of test is this" never arises for it.

  **SUPPRESSED** (§P7-18) is not a property of an executed row at all -- a suppressed
  row is one that is NOT EMITTED. It has no disposition because it has no standing.

THE ENFORCEMENT IS THE POINT, NOT THE LABEL
-------------------------------------------
A taxonomy that only annotates is a taxonomy that will be ignored the first time a
component row looks interesting. §P7-17 ruled that a COMPONENT-OF row appearing
standalone is an ERROR, so:

  * `attach_components` moves component rows INSIDE their parent's record;
  * `assert_no_component_standalone` RAISES if one reaches `stubs.json` or any other
    standalone list, and it is wired into `write_stubs`.

The failure it prevents is specific and it is the one that would actually happen: a
first-moment GLM row surviving BH on its own, being written up as a finding, and
nobody noticing that the multiplicity it faced was priced for a different test.
"""

from __future__ import annotations

import numpy as np

from . import mine as M

COMPONENT = "COMPONENT-OF"
REPLICATION = "REPLICATION-OF-DECLARED"
DECLARED = "DECLARED"

# §P7-18: THREE tags. This tuple is the taxonomy and nothing else belongs in it.
DISPOSITIONS = (COMPONENT, REPLICATION, DECLARED)

# Two further classes that are NOT dispositions -- see the module docstring for why
# the distinction is kept rather than collapsed for tidiness.
UNPRICED_CONTROL = "UNPRICED-CONTROL"     # §P7-16: a pricing class
SUPPRESSED = "SUPPRESSED"                 # §P7-18: not emitted, so not a row
ROW_CLASSES = DISPOSITIONS + (UNPRICED_CONTROL, SUPPRESSED)

DISPOSITION_RULE = {
    COMPONENT: ("exists only as part of another row's claim (F9-01's R1 beside R2; "
                "F9-04's Rayleigh beside Kuiper -- the comparison IS the claim). "
                "Owes no multiplicity because it makes no independent rejection, "
                "and MUST NOT be able to make one: it appears only inside its "
                "parent's record and never standalone."),
    REPLICATION: ("a re-execution of a row already declared elsewhere. Not "
                  "re-priced. Where the test-key digest matches the prior session, "
                  "p-invariance is a FREE §P6-5 determinism check and is wired."),
    DECLARED: ("a genuinely new independent row. These, and only these, are m."),
}

CLASS_RULE = dict(DISPOSITION_RULE)
CLASS_RULE[UNPRICED_CONTROL] = (
    "§P7-16: a control that can only calibrate a reference or condemn our own "
    "instrument. A PRICING class, not a disposition. Declared at m_s = 0.")
CLASS_RULE[SUPPRESSED] = (
    "§P7-18: a WITHIN-SESSION same-shape duplicate. NOT EMITTED, so it has no "
    "disposition -- it has no standing. Suppression is proved same-shape (feature, "
    "mark, statistic, time base) and mapped, and it LIFTS automatically if the "
    "shapes stop matching (e.g. the sub-daily arm turning on).")

# The parent statistics a first-moment GLM row is a COMPONENT of.
COMPONENT_PARENT_TESTS = ("second_circular_moment_score", "kuiper_V", "watson_U2")
DECLARED_TESTS = ("second_circular_moment_score", "kuiper_V", "watson_U2")


class ComponentRowStandalone(AssertionError):
    """A COMPONENT-OF row reached a standalone list. §P7-17 makes this an error."""


class MultipleDispositions(AssertionError):
    """A row carries zero or more than one disposition tag."""


class DeterminismViolation(AssertionError):
    """A REPLICATION row with a matching test-key digest returned a different p."""


# ------------------------------------------------------------------ tagging ---
def tag_rows(tests, prior_keys=None, cyclic_features=None):
    """Assign exactly one disposition to every row. Returns the rows, annotated.

    `prior_keys` is the set of `(test, feature, lag, mark)` tuples already declared
    in an earlier session; `cyclic_features` the names carrying a declared cycle (the
    ones F9-01/F9-04 run on, and therefore the ones whose first-moment GLM is a
    component rather than a hypothesis).
    """
    prior = set(prior_keys or ())
    cyc = set(cyclic_features or ())
    for t in tests:
        test = t.get("test")
        key = (test, t.get("feature"), t.get("lag"), t.get("mark"))
        if test in DECLARED_TESTS or t.get("mark_axis") == "F9-10":
            tag, parent, why = DECLARED, None, "new independent row; priced in m"
        elif t.get("control") or t.get("family") == 7:
            tag, parent, why = (UNPRICED_CONTROL, None,
                                "F7 observer control (§P7-16), m_s = 0")
        elif (test == "glm_poisson_offset_etas" and t.get("feature") in cyc
              and int(t.get("lag") or 0) == 0):
            tag = COMPONENT
            parent = "F9-01/F9-04 on %s" % t.get("feature")
            why = ("the first moment inside the second-moment claim, and the "
                   "Rayleigh form beside Kuiper -- the comparison IS the claim")
        elif key in prior:
            tag, parent, why = (REPLICATION, None,
                                "already declared in a PRIOR SESSION; not re-priced. "
                                "Cross-session only (§P7-18): a within-session "
                                "same-shape row is a DUPLICATE and is suppressed, "
                                "not tagged")
        else:
            tag, parent, why = (DECLARED, None,
                                "not a component, not a prior declaration -- "
                                "GENUINELY NEW and the integer must move for it")
        t["disposition"] = tag
        t["disposition_parent"] = parent
        t["disposition_reason"] = why
    return tests


def assert_one_disposition(tests):
    """Exactly one class per row, drawn from the three tags or the pricing class.

    SUPPRESSED never appears here: a suppressed row is not emitted, so if one
    reaches this function something upstream failed to suppress it.
    """
    leaked = [t for t in tests if t.get("disposition") == SUPPRESSED]
    if leaked:
        raise MultipleDispositions(
            "§P7-18: %d SUPPRESSED row(s) reached the emitted test list. A suppressed "
            "row is not emitted at all; if it is here, `suppression_map` was not "
            "applied. First: %r"
            % (len(leaked), {k: leaked[0].get(k)
                             for k in ("test", "feature", "mark")}))
    bad = [t for t in tests
           if t.get("disposition") not in (DISPOSITIONS + (UNPRICED_CONTROL,))]
    if bad:
        raise MultipleDispositions(
            "§P7-17: %d row(s) carry no valid disposition tag (first: %r). Every "
            "executed priceable row carries exactly one of %r; a control carries "
            "%r, which is a PRICING class and not a disposition (§P7-18)."
            % (len(bad), {k: bad[0].get(k) for k in ("test", "feature", "lag",
                                                     "mark", "disposition")},
               list(DISPOSITIONS), UNPRICED_CONTROL))
    return {"n_rows": len(tests), "ok": True}


def counts_by_disposition(tests):
    out = {d: 0 for d in ROW_CLASSES}
    for t in tests:
        out[t["disposition"]] = out.get(t["disposition"], 0) + 1
    return out


# ------------------------- §P7-18: within-session duplicates, PROVED and mapped -
# The declared shape of a mark row. Four coordinates, and the fourth is the one that
# does the work: two rows testing the same feature against the same mark with the
# same statistic are STILL different tests if one reads the feature on the day
# lattice and the other at event times. That is the entire content of F9-10's
# escape from the sinc, so it is a shape coordinate and not a footnote.
SHAPE_FIELDS = ("feature", "mark", "statistic", "time_base")
DAY_BINNED = "day-binned"
EVENT_TIMES = "event-times (sub-daily)"

SUPPRESSION_RULE = (
    "§P7-18: the v1 two-mark rows (`marks:<f>` on mag and depth) and the F9-10 mark "
    "axis rows for the same two marks are the SAME SHAPE while the sub-daily arm is "
    "OFF -- same feature, same mark, same statistic, same time base -- so the v1 "
    "rows are SUPPRESSED, not tagged. Tagging them REPLICATION would have been "
    "wrong twice over: a replication is CROSS-SESSION (§P7-17(2) as clarified), and "
    "a within-session same-shape row carries no information while reading as "
    "phantom replication. THE SUPPRESSION IS CONDITIONAL AND IT LIFTS ITSELF: turn "
    "the sub-daily arm on and the F9-10 row's time base becomes `event-times`, the "
    "shapes stop matching, and the v1 row is no longer a duplicate of anything.")


def row_shape(t, subdaily_default=False):
    """The declared shape of a row. Declare-then-prove: this IS the declaration."""
    if t.get("mark_axis") == "F9-10":
        tb = EVENT_TIMES if t.get("subdaily", subdaily_default) else DAY_BINNED
    else:
        # the v1 mark axis reads the feature at its DAY value, always
        tb = DAY_BINNED
    return {"feature": t.get("feature"), "mark": t.get("mark"),
            "statistic": t.get("test"), "time_base": tb}


def same_shape(a, b, subdaily_default=False):
    """True when two rows are the same test. The PROOF half of declare-then-prove."""
    sa, sb = row_shape(a, subdaily_default), row_shape(b, subdaily_default)
    return all(sa[f] == sb[f] for f in SHAPE_FIELDS), sa, sb


def suppression_map(tests, subdaily=False):
    """Suppress within-session same-shape duplicates. Returns (kept, suppressed, map).

    The SURVIVOR is always the F9-10 row and the suppressed one always the v1 row,
    and the direction is not arbitrary: the F9-10 row is the one this tranche
    DECLARED and priced, and the v1 row is the one the engine emits by inheritance.
    Suppressing the declared row and keeping the inherited one would leave the
    tranche's own priced slot unfilled.
    """
    declared = {}
    for t in tests:
        if t.get("mark_axis") == "F9-10":
            sh = row_shape(t, subdaily)
            declared[tuple(sh[f] for f in SHAPE_FIELDS)] = t
    kept, suppressed, mapping = [], [], []
    for t in tests:
        if t.get("mark_axis") == "F9-10" or t.get("mark") is None:
            kept.append(t)
            continue
        sh = row_shape(t, subdaily)
        key = tuple(sh[f] for f in SHAPE_FIELDS)
        parent = declared.get(key)
        if parent is None:
            kept.append(t)
            continue
        ok, sa, sb = same_shape(t, parent, subdaily)
        assert ok, "suppression_map matched two rows of different shape"
        t["disposition"] = SUPPRESSED
        suppressed.append(t)
        mapping.append({
            "suppressed": {"source": "v1 mark axis", "shape": sa},
            "survivor": {"source": "F9-10 mark axis (declared and priced)",
                         "shape": sb},
            "proof": ("all %d shape coordinates equal: %s"
                      % (len(SHAPE_FIELDS),
                         ", ".join("%s=%r" % (f, sa[f]) for f in SHAPE_FIELDS))),
            "lifts_if": ("the sub-daily arm turns on: the survivor's time_base "
                         "becomes %r, the shapes differ, and this suppression is "
                         "void" % EVENT_TIMES),
        })
    return kept, suppressed, {
        "rule": SUPPRESSION_RULE,
        "subdaily_arm": bool(subdaily),
        "active": (not subdaily),
        "n_suppressed": len(suppressed),
        "shape_fields": list(SHAPE_FIELDS),
        "entries": mapping,
    }


# -------------------------------------------------------------- enforcement ---
def attach_components(tests):
    """Move COMPONENT-OF rows INSIDE their parent's record. Returns (kept, moved).

    The parent is the set of DECLARED rows on the same feature. A component is
    attached to every parent it is a component of -- F9-01's second moment and
    F9-04's two omnibus rows all need the first moment beside them, and giving it to
    one of them arbitrarily would make the other two's reading depend on which
    statistic happened to be listed first.
    """
    comps = [t for t in tests if t.get("disposition") == COMPONENT]
    kept = [t for t in tests if t.get("disposition") != COMPONENT]
    by_feature = {}
    for c in comps:
        by_feature.setdefault(c.get("feature"), []).append(c)
    for t in kept:
        if t.get("test") not in COMPONENT_PARENT_TESTS:
            continue
        for c in by_feature.get(t.get("feature"), []):
            t.setdefault("components", []).append({
                "test": c.get("test"), "feature": c.get("feature"),
                "lag": c.get("lag"), "p_raw": c.get("p_raw"),
                "amplitude_log_rate": c.get("amplitude_log_rate"),
                "chi2_score": c.get("chi2_score"),
                "disposition": COMPONENT,
                "why": ("F9-01's Pit: R2 is uninterpretable without R1 beside it. "
                        "F9-04: the comparison against the 2-df form IS the claim. "
                        "This row is half of the parent's claim, not a claim."),
            })
    return kept, comps


def assert_no_component_standalone(rows, context):
    """§P7-17: a COMPONENT-OF row in a standalone list is an ERROR, not a warning."""
    bad = [r for r in rows if r.get("disposition") == COMPONENT]
    if bad:
        raise ComponentRowStandalone(
            "§P7-17 VIOLATED: %d COMPONENT-OF row(s) reached %s, which is a "
            "STANDALONE list. A component makes no independent rejection and owes no "
            "multiplicity precisely because it cannot be read alone; a component "
            "surviving BH on its own and being written up as a finding would face a "
            "threshold priced for a different test. First offender: %r. Attach it to "
            "its parent (`attach_components`) or exclude it."
            % (len(bad), context,
               {k: bad[0].get(k) for k in ("test", "feature", "lag",
                                           "disposition_parent")}))
    return {"context": context, "n_rows": len(rows), "ok": True}


# ---------------------------------------- the free §P6-5 determinism dividend --
def _digest(seed, t):
    kind = {"glm_poisson_offset_etas": "glm"}.get(t.get("test"), t.get("test"))
    if t.get("mark"):
        kind = "mark_%s" % t["mark"]
    return M.test_key_digest(M.test_key(seed, t.get("feature"), kind,
                                        lag=t.get("lag"),
                                        null_type="block_bootstrap"))


def replication_invariance(tests, prior_tests, seed, prior_seed):
    """§P7-17: p-invariance of REPLICATION rows against the prior session.

    THE CHECK IS CONDITIONAL AND THE CONDITION IS THE INTERESTING PART. A test's
    random stream is addressed by its COMPLETE key, and `master_seed` is a field of
    that key -- so two sessions at different master seeds draw different surrogates
    for the same hypothesis and their p-values are not supposed to match. Comparing
    them and reporting a difference would be reporting the seed.

    So rows are matched on test-key DIGEST:

      * digest MATCHES  -> the two runs addressed the identical stream, and the p
        MUST be bitwise identical. A difference here is a determinism failure and
        raises: it means the same declared key produced two different answers, which
        is the §P6-5 invariant the whole reproducibility argument rests on.
      * digest DIFFERS  -> NOT COMPARABLE, with the differing field named. Reported,
        never silently counted as agreement and never counted as a failure.
    """
    prior_by_digest = {}
    prior_by_name = {}
    for p in prior_tests:
        prior_by_digest[_digest(prior_seed, p)] = p
        prior_by_name[(p.get("test"), p.get("feature"), p.get("lag"),
                       p.get("mark"))] = p
    checked, not_comparable, mismatches = [], [], []
    for t in tests:
        if t.get("disposition") != REPLICATION:
            continue
        d = _digest(seed, t)
        name = (t.get("test"), t.get("feature"), t.get("lag"), t.get("mark"))
        if d in prior_by_digest:
            p_now = float(t.get("p_raw", float("nan")))
            p_then = float(prior_by_digest[d].get("p_raw", float("nan")))
            same = (p_now == p_then)
            rec = {"key": name, "digest": d, "p_now": p_now, "p_prior": p_then,
                   "bitwise_identical": bool(same)}
            checked.append(rec)
            if not same:
                mismatches.append(rec)
        elif name in prior_by_name:
            not_comparable.append({
                "key": name,
                "reason": ("test-key digest differs -- the two sessions declare "
                           "different master seeds (%s vs %s), so they address "
                           "different streams and their p-values are not supposed "
                           "to agree. NOT COMPARABLE, not a failure."
                           % (seed, prior_seed))})
    if mismatches:
        raise DeterminismViolation(
            "§P6-5 VIOLATED: %d REPLICATION row(s) share a test-key digest with the "
            "prior session but returned a DIFFERENT p. The same declared key must "
            "address the same stream and produce the same answer. First: %r"
            % (len(mismatches), mismatches[0]))
    return {
        "n_replication_rows": sum(1 for t in tests
                                  if t.get("disposition") == REPLICATION),
        "n_digest_matched_and_checked": len(checked),
        "n_bitwise_identical": sum(1 for c in checked if c["bitwise_identical"]),
        "n_not_comparable_different_seed": len(not_comparable),
        "rows": checked[:200],
        "not_comparable": not_comparable[:50],
        "verdict": ("EXERCISED -- %d rows re-derived bitwise" % len(checked)
                    if checked else
                    "NOT EXERCISED under this config: the two sessions declare "
                    "different master seeds, so no digest matches. The check costs "
                    "nothing and would become live the moment a replication session "
                    "declares the prior session's master seed."),
        "rule": CLASS_RULE[REPLICATION],
    }


# ------------------------------------------------------------- the enumeration -
def enumerate_session_rows(cyclic_features, other_science_features,
                           control_features, n_period_peaks, other_lags,
                           mark_axis_features, scored_marks, v1_marks=("mag",
                                                                       "depth"),
                           prior_declared=None):
    """Every row a Tranche B session would execute, categorised. NO DATA REQUIRED.

    Symbolic on purpose: B's run is gated, so the enumeration has to be derivable
    from the declared feature names and the frozen config alone. Everything below is
    a count of rows the engine WOULD emit, not a count of rows it did.
    """
    rows = []

    def add(n, test, tag, parent, why, detail):
        rows.append({"n": int(n), "test": test, "disposition": tag,
                     "parent": parent, "why": why, "detail": detail})

    n_cyc = len(cyclic_features)
    n_other = len(other_science_features)
    n_ctl = len(control_features)

    # ---- the three DECLARED statistics --------------------------------------
    add(n_cyc, "second_circular_moment_score", DECLARED, None,
        "F9-01, priced", "%d cyclic features x 1 statistic" % n_cyc)
    add(2 * n_cyc, "kuiper_V / watson_U2", DECLARED, None,
        "F9-04, priced", "%d cyclic features x 2 statistics" % n_cyc)
    add(len(scored_marks) * mark_axis_features, "markx", DECLARED, None,
        "F9-10, priced",
        "%d scored marks x %d features" % (len(scored_marks), mark_axis_features))

    # ---- COMPONENT-OF -------------------------------------------------------
    add(n_cyc, "glm_poisson_offset_etas", COMPONENT,
        "F9-01/F9-04 on the same feature",
        "the first moment inside the second-moment claim; the Rayleigh form beside "
        "Kuiper -- the comparison IS the claim",
        "%d cyclic features x lag 0" % n_cyc)

    # ---- REPLICATION-OF-DECLARED --------------------------------------------
    add(n_other * other_lags, "glm_poisson_offset_etas", REPLICATION, None,
        "already declared (K-089-R tranche 1 / Tranche A's 550); not re-priced",
        "%d non-cyclic science features x %d lags" % (n_other, other_lags))
    add((n_cyc + n_other) * len(v1_marks), "spearman / circular-linear", SUPPRESSED,
        None,
        "§P7-18: SAME SHAPE as the F9-10 `mag`/`depth` rows in this same session "
        "while the sub-daily arm is off (same feature, mark, statistic, time base). "
        "A within-session same-shape row is a DUPLICATE, not a replication: it is "
        "SUPPRESSED with a proved shape map, and the suppression lifts by itself if "
        "the sub-daily arm turns on.",
        "%d features x %d v1 marks" % (n_cyc + n_other, len(v1_marks)))
    add(n_period_peaks, "lomb_scargle_peak", REPLICATION, None,
        "the period scan runs in every session and is already declared",
        "%d declared peaks" % n_period_peaks)

    # ---- UNPRICED-CONTROL ---------------------------------------------------
    add(n_ctl, "glm_poisson_offset_etas", UNPRICED_CONTROL, None,
        "F7 observer controls (§P7-16): can only calibrate a reference or condemn "
        "our own instrument", "%d count-path control features x lag 0" % n_ctl)

    totals = {d: sum(r["n"] for r in rows if r["disposition"] == d)
              for d in ROW_CLASSES}
    totals["TOTAL"] = sum(r["n"] for r in rows)
    return {"rows": rows, "totals": totals,
            "rule": dict(CLASS_RULE),
            "n_genuinely_new_beyond_the_declared_statistics": 0,
            "genuinely_new_note": (
                "Every executed row is a COMPONENT-OF, a (cross-session) "
                "REPLICATION-OF-DECLARED, an UNPRICED-CONTROL, a SUPPRESSED "
                "within-session duplicate, or one of the DECLARED statistics. NO row "
                "falls outside the taxonomy, so no genuinely-new independent "
                "hypothesis is executed beyond the priced ones and the integer does "
                "not move on this account.")}
