"""EQ-23 predictability engine (v1). See SPEC.md."""

__version__ = "1.1.0"

BASELINE_CAVEAT = (
    "baseline=climatology-v1 (clustering NOT absorbed; sniff-grade only)"
)

ETAS_CAVEAT = (
    "baseline=etas-v1 (isotropic kernel; anisotropy/mechanism NOT absorbed)"
)

# Canonical baseline names, and the aliases the CLI accepts.
BASELINE_ALIASES = {
    "climatology": "climatology-v1",
    "climatology-v1": "climatology-v1",
    "etas": "etas",
    "etas-v1": "etas",
}

_CAVEATS = {"climatology-v1": BASELINE_CAVEAT, "etas": ETAS_CAVEAT}


def canonical_baseline(name: str) -> str:
    try:
        return BASELINE_ALIASES[name]
    except KeyError:
        raise KeyError(
            f"unknown baseline {name!r}; have {sorted(set(BASELINE_ALIASES))}"
        ) from None


def baseline_caveat(name: str) -> str:
    """The caveat line that MUST be printed on every report, per baseline.

    climatology-v1 keeps the v1 line (clustering not absorbed). etas-v1 absorbs
    Omori-Utsu clustering through an isotropic space-time kernel and says so --
    and says what it still does NOT absorb.
    """
    return _CAVEATS[canonical_baseline(name)]
