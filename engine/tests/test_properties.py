"""B1 / §S2: the per-event property join, and the four refusals it enforces in code.

  * `property_class` is REQUIRED (§P7-24 SP-2 / §S11(4)) and a column without one
    cannot be built;
  * a property whose class has no null layer ATTACHED may rank and may not promote
    (SP-2's last sentence, as a query rather than a review-time judgement);
  * a LEVEL property cannot exist without its dwell-time occupancy measure, and the
    PIT transform is checked against §P7-23(C)'s ANALYTIC arcsine case -- 0.318 at
    mid-level, 1.019 at |x| = 0.95, and exactly 1/3 of the measure in the lowest
    quarter of range;
  * provenance is a FIELD (§S2.2(3)): a column missing any of the five fields raises.

Plus the property that separates this module from the count path at all: every value
is evaluated at the event's own sub-daily time (§S2.2(1)), so two events on the same
calendar day get DIFFERENT phases -- which is exactly what `obs_utc_hour_phase` cannot
do on the count path by construction.
"""

from __future__ import annotations

import datetime as _dt

import numpy as np
import pytest

from engine import properties as P

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

T0 = _dt.datetime(2000, 1, 1)
SITE = dict(site_lat=38.0, site_lon=142.0)


def _cols(**kw):
    kw.setdefault("families", ("solid_tide", "ephemeris", "human_clock", "season"))
    return P.build_property_matrix(T0, 0.0, 400.0, sample_minutes=60.0,
                                   **SITE, **kw)


# ------------------------------------------------------- SP-2's required class --
def test_a_property_without_a_class_cannot_be_built():
    with pytest.raises(P.PropertyClassMissing):
        P.PropertyColumn("x", "season", "phase", None, lambda t: t,
                         {k: "" for k in P.PROVENANCE_FIELDS})


def test_a_property_with_an_unknown_class_cannot_be_built():
    with pytest.raises(P.UnknownPropertyClass):
        P.PropertyColumn("x", "season", "phase", "vibes", lambda t: t,
                         {k: "" for k in P.PROVENANCE_FIELDS})


def test_the_class_table_is_sp2s_table():
    assert set(P.PROPERTY_CLASSES) == {
        "human-schedule", "level-waveform-phase", "clustering-derived",
        "count-path-periodic", "catalogue-endogenous"}
    for cls, meta in P.PROPERTY_CLASSES.items():
        assert meta["layer"] and meta["attach"], cls


def test_human_schedule_cannot_promote_until_the_f7_layer_is_attached():
    cols, _a = _cols()
    dow = [c for c in cols if c.name == "day_of_week"][0]
    ok, why = dow.may_promote()
    assert not ok
    assert "human-schedule" in why and "not attached" in why
    dow.attach("observer_features")
    ok2, why2 = dow.may_promote()
    assert ok2 and "attached" in why2


def test_the_dwell_layer_is_auto_attached_and_the_measured_layers_are_not():
    """Only the layer a column carries BY CONSTRUCTION is auto-attached."""
    cols, _a = _cols()
    for c in cols:
        if c.pclass == "level-waveform-phase":
            assert "dwell" in c.attached, c.name
        else:
            assert not c.attached, c.name


# --------------------------------------------------------- provenance is a field --
def test_a_column_missing_provenance_fields_raises():
    with pytest.raises(P.ProvenanceMissing):
        P.PropertyColumn("x", "season", "phase", "level-waveform-phase",
                         lambda t: t, {"source": "s"})


def test_every_built_column_carries_all_five_provenance_fields():
    cols, audit = _cols()
    for rec in audit["columns"]:
        for f in P.PROVENANCE_FIELDS:
            assert rec["provenance"].get(f), (rec["property"], f)


def test_the_tide_columns_carry_the_frozen_scalar_and_convention():
    """§P7-23(D): the quadrant is convention-free ONLY once the scalar is fixed."""
    cols, _a = _cols(families=("solid_tide",))
    q = [c for c in cols if c.name == "tide_quadrant"][0]
    assert q.provenance["convention"].startswith("D0-tanaka-stressmax-v1")
    assert q.provenance["scalar_definition"] == "areal_strain"


# ---------------------------------------------------------- the dwell-time layer --
def test_a_level_property_cannot_exist_without_its_dwell_measure():
    with pytest.raises(P.DwellDensityMissing):
        P.PropertyColumn("x", "season", "level", "level-waveform-phase",
                         lambda t: t, {k: "v" for k in P.PROVENANCE_FIELDS})


def test_arcsine_cdf_reproduces_p7_23_C_exactly():
    """The three numbers §P7-23(C) verified: the density and the lowest quarter."""
    # density f(x) = dF/dx = 1/(pi sqrt(1-x^2))
    for x, want in ((0.0, 0.3183), (0.95, 1.0191)):
        h = 1e-6
        f = (P.arcsine_cdf(x + h) - P.arcsine_cdf(x - h)) / (2 * h)
        assert f == pytest.approx(want, abs=2e-3), x
    # exactly 1/3 of the measure lies in the lowest quarter of the range [-1, 1]
    assert P.arcsine_cdf(-0.5) == pytest.approx(1.0 / 3.0, abs=1e-9)


def test_the_empirical_dwell_of_a_sinusoid_is_the_arcsine_measure():
    """The transform is checked against the case whose answer is known in closed form."""
    theta = np.linspace(0.0, 200.0 * 2 * np.pi, 400001)
    dwell = P.empirical_dwell(np.sin(theta))
    assert P.dwell_lowest_quarter_fraction(dwell) == pytest.approx(1.0 / 3.0,
                                                                  abs=5e-3)
    for x in (-0.9, -0.3, 0.0, 0.4, 0.95):
        emp = float(np.searchsorted(dwell, x) / dwell.size)
        assert emp == pytest.approx(float(P.arcsine_cdf(x)), abs=5e-3), x


def test_dwell_pit_makes_uniform_phase_events_uniform_in_the_transform():
    """The whole point: an UNCORRECTED level piles events at the extremes; the PIT
    of the SAME events is flat, so a Kuiper on it measures concentration and not the
    occupancy measure."""
    rng = np.random.default_rng(0)
    grid = np.sin(np.linspace(0.0, 400.0 * np.pi, 200001))
    dwell = P.empirical_dwell(grid)
    events = np.sin(rng.uniform(0.0, 2 * np.pi, 20000))          # uniform PHASE
    # raw levels: a third of them in the lowest quarter of range -- the artifact
    raw_low = float(np.mean(events <= -0.5))
    assert raw_low == pytest.approx(1.0 / 3.0, abs=0.02)
    assert raw_low > 0.30                       # would read as "concentration"
    pit = P.dwell_pit(events, dwell)
    assert float(np.mean(pit <= 0.25)) == pytest.approx(0.25, abs=0.02)
    assert float(np.mean(pit <= 0.50)) == pytest.approx(0.50, abs=0.02)


def test_dwell_pit_refuses_a_missing_dwell():
    with pytest.raises(P.DwellDensityMissing):
        P.dwell_pit(np.array([0.1]), None)


def test_level_columns_report_their_lowest_quarter_fraction():
    _c, audit = _cols(families=("solid_tide",))
    lv = [r for r in audit["columns"] if r["property"] == "tide_level"][0]
    assert lv["dwell_time_corrected"] is True
    assert 0.25 < lv["dwell_lowest_quarter_fraction"] < 0.55
    assert lv["arcsine_reference_lowest_quarter"] == pytest.approx(1.0 / 3.0)


# ----------------------------------------------------------- the sub-daily path --
def test_two_events_on_the_same_calendar_day_get_different_phases():
    """§S2.2(1): the count path CANNOT do this -- a daily bin has one hour phase."""
    cols, _a = _cols()
    same_day = np.array([10.05, 10.30, 10.80])
    for c in cols:
        if not c.subdaily:
            continue
        v = np.asarray(c.evaluate(same_day), dtype=np.float64)
        assert np.unique(v[np.isfinite(v)]).size > 1, c.name


def test_day_binned_columns_say_so_per_column():
    _c, audit = _cols()
    by = {r["property"]: r["subdaily"] for r in audit["columns"]}
    assert by["day_of_week"] is False and by["month_of_year"] is False
    assert by["tide_phase"] is True and by["utc_hour_phase"] is True


# ----------------------------------------------------- the calendar arithmetic --
def test_day_of_week_matches_the_calendar_over_a_thirty_year_span():
    cols = P.human_clock_columns(T0, 0.0)
    dow = [c for c in cols if c.name == "day_of_week"][0]
    for k in range(0, 11000, 53):
        want = (T0 + _dt.timedelta(days=k)).strftime("%A")
        got = P.DOW_NAMES[int(dow.evaluate(np.array([k + 0.5]))[0])]
        assert got == want, (k, want, got)


def test_month_of_year_matches_the_calendar():
    cols = P.human_clock_columns(T0, 0.0)
    mon = [c for c in cols if c.name == "month_of_year"][0]
    for k in range(0, 4000, 17):
        want = (T0 + _dt.timedelta(days=k)).month
        assert int(mon.evaluate(np.array([k + 0.5]))[0]) + 1 == want, k


def test_day_of_week_is_computed_by_the_same_code_for_nulls_and_observations():
    """One code path, so a one-day offset cannot enter through a datetime round-trip.

    Checked on the CALENDAR FUNCTIONS' own bodies, not on the enclosing builder: the
    builder's signature carries a `_dt.datetime` annotation for `t0`, which is a type
    and not a round-trip, and a text search over the whole source would fail on it.
    """
    import ast
    import inspect
    tree = ast.parse(inspect.getsource(P.human_clock_columns).lstrip())
    inner = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    assert {"dow_fn", "month_fn"} <= set(inner)
    for name in ("dow_fn", "month_fn"):
        code = chr(10).join(ast.unparse(x) for x in ast.walk(inner[name])
                            if isinstance(x, (ast.Call, ast.Attribute)))
        for forbidden in ("datetime", "strftime", "weekday", "isoweekday",
                          "fromtimestamp", "timedelta"):
            assert forbidden not in code, (name, forbidden)


def test_the_calendar_columns_are_vectorised_and_order_independent():
    """A per-element calendar path would make a null array and an observed array
    two different code paths; a vectorised one cannot."""
    cols = P.human_clock_columns(T0, 0.0)
    t = np.array([0.4, 900.9, 17.2, 4000.6, 123.1])
    for c in cols:
        one = np.array([float(c.evaluate(np.array([x]))[0]) for x in t])
        many = np.asarray(c.evaluate(t), dtype=np.float64)
        assert np.allclose(one, many), c.name
        o = np.argsort(t)
        assert np.allclose(np.asarray(c.evaluate(t[o]))[np.argsort(o)], many), c.name


# ------------------------------------------------------------- statistic values --
def test_phase_columns_return_angles_and_level_columns_return_pit_angles():
    cols, _a = _cols()
    t = np.linspace(1.0, 399.0, 300)
    for c in cols:
        v = np.asarray(c.statistic_values(t), dtype=np.float64)
        v = v[np.isfinite(v)]
        if c.ptype in ("phase", "level"):
            assert v.min() >= 0.0 and v.max() < 2 * np.pi + 1e-9, c.name
        else:
            assert set(np.unique(v)) <= set(map(float, c.categories)), c.name


def test_the_level_statistic_path_is_the_dwell_corrected_one():
    """A level's statistic values are UNIFORM under a uniform-phase draw, not arcsine."""
    cols, _a = _cols(families=("solid_tide",))
    lvl = [c for c in cols if c.name == "tide_level"][0]
    rng = np.random.default_rng(3)
    t = np.sort(rng.uniform(1.0, 399.0, 5000))
    ang = np.asarray(lvl.statistic_values(t)) / (2 * np.pi)
    assert float(np.mean(ang <= 0.25)) == pytest.approx(0.25, abs=0.03)


# ------------------------------------------------- absence rather than substitution --
def test_an_unbuilt_family_is_absent_and_says_why():
    _c, audit = P.build_property_matrix(T0, 0.0, 60.0, 38.0, 142.0,
                                        families=("season", "clocks"),
                                        sample_minutes=120.0)
    assert "NOT BUILT" in audit["family_notes"]["clocks"]
    assert "§P7-5(2)" in audit["family_notes"]["clocks"]


def test_an_unknown_family_raises_rather_than_being_ignored():
    with pytest.raises(ValueError):
        P.build_property_matrix(T0, 0.0, 60.0, 38.0, 142.0, families=("aurora",))


# ------------------------------------------------------------ redundancy collapse --
def test_redundancy_collapse_keeps_the_earlier_family_representative():
    """§S3.5: chosen by a PRE-DECLARED rule, never by which one scored better."""
    cols, _a = _cols()
    t = np.linspace(1.0, 399.0, 2000)
    kept, rep = P.redundancy_collapse(cols, t, threshold=0.90)
    names = {c.name for c in kept}
    # solar_annual_phase (family 'ephemeris') and day_of_year_phase (family 'season')
    # are the SAME quantity; ephemeris is earlier in FAMILY_ORDER, so it survives.
    assert ("solar_annual_phase" in names) and ("day_of_year_phase" not in names)
    assert any(p["dropped"] == "day_of_year_phase" for p in rep["collapsed"])
    assert rep["n_after"] < rep["n_before"]
    assert "PRE-DECLARED" in rep["rule"]


def test_redundancy_collapse_never_touches_categoricals():
    cols, _a = _cols()
    kept, _r = P.redundancy_collapse(cols, np.linspace(1.0, 399.0, 500))
    for c in cols:
        if c.ptype == "categorical":
            assert c in kept


# --------------------------------------------------------------- the declaration --
def test_the_declaration_digest_covers_the_columns_and_not_their_values():
    cols, a1 = _cols()
    _c2, a2 = _cols()
    assert P.property_declaration_digest(a1) == P.property_declaration_digest(a2)
    a3 = dict(a1)
    a3["columns"] = a1["columns"][:-1]
    assert P.property_declaration_digest(a3) != P.property_declaration_digest(a1)


def test_the_declared_magnitude_strata_are_the_searcher_s_four():
    assert P.DECLARED_MAG_STRATA == (4.5, 5.0, 5.5, 6.0)
    assert "NEW SEED" in P.MAG_STRATA_RULE
    assert P.HUMAN_SCHEDULE_MC_FLOOR == 6.0


def test_the_site_rule_is_declared_on_every_site_dependent_column():
    cols, audit = _cols()
    assert "REGION'S DECLARED CENTROID" in audit["site"]["rule"]
    for c in cols:
        if c.name in ("tide_phase", "tide_level", "local_solar_hour_phase"):
            assert any("centroid" in n.lower() for n in c.notes), c.name


# ------------------------------------------------------- family G: the marks --
def _mark_events(n=300, seed=0):
    rng = np.random.default_rng(seed)
    return dict(day_float=np.sort(rng.uniform(0.0, 400.0, n)),
                mag=rng.uniform(4.5, 7.0, n), depth=rng.uniform(0.0, 200.0, n),
                lat=rng.uniform(30.0, 40.0, n), lon=rng.uniform(130.0, 140.0, n))


def test_mark_columns_carry_sp2s_clustering_and_endogenous_classes():
    """The two classes whose null layer is a MEASUREMENT, not a construction."""
    cols, audit = P.mark_columns(_mark_events())
    by = {c.name: c for c in cols}
    assert by["depth_km"].pclass == "catalogue-endogenous"
    assert by["dt_prior_days"].pclass == "clustering-derived"
    assert by["dist_nearest_prior_km"].pclass == "clustering-derived"
    assert by["cluster_member"].pclass == "clustering-derived"
    assert set(audit["built_marks"]) >= {"depth", "dt_prior_days", "cluster_member"}


def test_no_mark_column_can_promote_until_a_mark_null_is_attached():
    """SP-2's disposition, not an omission: they RANK, they do not PROMOTE."""
    cols, _a = P.mark_columns(_mark_events())
    for c in cols:
        ok, why = c.may_promote()
        assert not ok, c.name
        assert "mark_null" in why
    cols[0].attach("mark_null")
    assert cols[0].may_promote()[0]


def test_mark_columns_take_EVENT_INDICES_not_times_and_say_so():
    """This family is not a function of time, and its signature reflects that."""
    ev = _mark_events()
    cols, _a = P.mark_columns(ev)
    depth = [c for c in cols if c.name == "depth_km"][0]
    got = np.asarray(depth.evaluate(np.array([0, 5, 17])), dtype=np.float64)
    assert np.allclose(got, np.asarray(ev["depth"])[[0, 5, 17]])
    assert any("EVENT INDICES" in n for n in depth.notes)
    assert depth.subdaily is False


def test_the_tranche_b_VIF_mark_warning_travels_on_the_two_marks_that_earned_it():
    """§S2.1 G: dt_prior and cluster_member were measured ABOVE the 4.575 pooled
    fallback, and the warning is on the column rather than in someone's memory."""
    cols, _a = P.mark_columns(_mark_events())
    by = {c.name: c for c in cols}
    for name in ("dt_prior_days", "cluster_member"):
        assert any("4.575" in n for n in by[name].notes), name
    assert not any("4.575" in n for n in by["depth_km"].notes)


def test_a_mark_dwell_is_labelled_as_NOT_the_arcsine_construction():
    """A catalogue mark has no time grid to dwell on; conflating the two would read
    §P7-23(C)'s argument into a place it does not apply."""
    cols, _a = P.mark_columns(_mark_events())
    depth = [c for c in cols if c.name == "depth_km"][0]
    assert depth.dwell is not None
    assert any("not a time-occupancy measure" in n for n in depth.notes)


def test_mark_columns_statistic_values_are_dwell_corrected_and_uniform():
    cols, _a = P.mark_columns(_mark_events(2000, seed=3))
    depth = [c for c in cols if c.name == "depth_km"][0]
    ang = np.asarray(depth.statistic_values(np.arange(2000))) / (2 * np.pi)
    assert float(np.mean(ang <= 0.25)) == pytest.approx(0.25, abs=0.03)
    assert float(np.mean(ang <= 0.75)) == pytest.approx(0.75, abs=0.03)


def test_the_marks_family_refuses_to_build_without_an_event_table():
    with pytest.raises(ValueError):
        P.build_property_matrix(T0, 0.0, 60.0, 38.0, 142.0, families=("marks",),
                                sample_minutes=120.0)
