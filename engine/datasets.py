"""Load replication/data/comcat_world/*.csv into one deduped catalog.

Column order varies by source file -> always select by NAME.
pandas 3.0: parse `time` with utc=True.
"""

from __future__ import annotations

import glob
import os

import numpy as np
import pandas as pd

DEFAULT_DATA_DIR = os.path.join("data", "comcat_world")

# The comcat_world boxes were downloaded with a stated M>=4.5 floor. Some ComCat
# exports carry sub-floor spurious rows; we filter to the stated floor and report.
CATALOG_MAG_FLOOR = 4.5

NEEDED = ["time", "latitude", "longitude", "depth", "mag", "id", "type"]


class LoadReport(dict):
    """Counted invariants of a catalog load; printed by the CLI."""

    def lines(self):
        out = ["catalog load:"]
        for f, n in self["per_file"]:
            out.append(f"  {f:<24s} rows={n}")
        out.append(f"  rows_read_total      = {self['rows_read']}")
        out.append(f"  dropped_non_earthquake = {self['dropped_type']}")
        out.append(f"  dropped_below_floor(M<{CATALOG_MAG_FLOOR}) = {self['dropped_floor']}")
        out.append(f"  dropped_bad_geometry = {self['dropped_bad']}")
        out.append(f"  duplicate_ids_removed = {self['dup_removed']}")
        out.append(f"  n_events             = {self['n_events']}")
        out.append(f"  time_span            = {self['t_min']} .. {self['t_max']}")
        out.append(f"  span_days            = {self['span_days']}")
        return out


def load_catalog(data_dir: str = DEFAULT_DATA_DIR, mag_floor: float = CATALOG_MAG_FLOOR):
    """Return (DataFrame, LoadReport). DataFrame columns: t, lat, lon, depth, mag, id, src.

    Sorted by time ascending, deduped by ComCat event id (boxes overlap).
    """
    files = sorted(glob.glob(os.path.join(data_dir, "*.csv")))
    if not files:
        raise FileNotFoundError(f"no CSVs under {data_dir!r}")

    frames = []
    per_file = []
    for path in files:
        df = pd.read_csv(path)
        missing = [c for c in NEEDED if c not in df.columns]
        if missing:
            raise ValueError(f"{path}: missing columns {missing}")
        df = df[NEEDED].copy()          # select by NAME; order varies by source
        df["src"] = os.path.basename(path)
        per_file.append((os.path.basename(path), len(df)))
        frames.append(df)

    cat = pd.concat(frames, ignore_index=True)
    rows_read = len(cat)
    if rows_read == 0:
        raise ValueError("bulk transform produced 0 rows (concat)")

    n0 = len(cat)
    cat = cat[cat["type"].astype("string") == "earthquake"]
    dropped_type = n0 - len(cat)

    n0 = len(cat)
    cat = cat[cat["mag"].astype("float64") >= mag_floor]
    dropped_floor = n0 - len(cat)

    # pandas 3.0: format='mixed' tolerates the two timestamp spellings in these files.
    t = pd.to_datetime(cat["time"], utc=True, format="mixed")
    cat = cat.assign(t=t)

    n0 = len(cat)
    ok = (
        cat["t"].notna()
        & cat["latitude"].between(-90, 90)
        & cat["longitude"].between(-180, 180)
        & cat["mag"].notna()
    )
    cat = cat[ok]
    dropped_bad = n0 - len(cat)

    n_before = len(cat)
    cat = cat.drop_duplicates(subset="id", keep="first")
    dup_removed = n_before - len(cat)

    cat = (
        cat.rename(columns={"latitude": "lat", "longitude": "lon"})
        .loc[:, ["t", "lat", "lon", "depth", "mag", "id", "src"]]
        .sort_values("t")
        .reset_index(drop=True)
    )

    n = len(cat)
    if n == 0:
        raise ValueError("bulk transform produced 0 rows (catalog load)")

    # --- Invariants, asserted loudly (SPEC "Invariants") ---
    if mag_floor == CATALOG_MAG_FLOOR:
        # size invariant applies to the catalogue at its stated floor
        assert 50_000 < n < 150_000, (
            f"catalog size {n} outside expected 50k..150k for comcat_world")
    assert dup_removed > 0, (
        "dedup by event id removed nothing; overlapping boxes should share events "
        "-- loader or data is wrong"
    )
    assert float(cat["mag"].min()) >= mag_floor

    rep = LoadReport(
        per_file=per_file,
        rows_read=rows_read,
        dropped_type=int(dropped_type),
        dropped_floor=int(dropped_floor),
        dropped_bad=int(dropped_bad),
        dup_removed=int(dup_removed),
        n_events=int(n),
        t_min=str(cat["t"].iloc[0]),
        t_max=str(cat["t"].iloc[-1]),
        span_days=int((cat["t"].iloc[-1] - cat["t"].iloc[0]).days),
        files=[f for f, _ in per_file],
        mag_floor=float(mag_floor),
    )
    return cat, rep


def catalog_arrays(cat: pd.DataFrame):
    """Plain numpy view of a catalog: (t_days_float, lat, lon, mag)."""
    t0 = cat["t"].iloc[0].normalize()
    days = (cat["t"] - t0).dt.total_seconds().to_numpy(dtype="float64") / 86400.0
    return (
        days,
        cat["lat"].to_numpy(dtype="float64"),
        cat["lon"].to_numpy(dtype="float64"),
        cat["mag"].to_numpy(dtype="float64"),
        t0,
    )
