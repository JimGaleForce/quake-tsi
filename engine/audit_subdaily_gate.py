"""D-5 runner -- point `observer`'s sub-daily gate at a real session and read it out.

§K92-1 D-5: *"Clear the F9-10 sub-daily gate: F7-01/02/03 observer controls.
`observer.assert_subdaily_gate` hard-refuses until this exists. Nothing sub-daily runs
before it."* The controls were built in Tranche A; the debt this discharges is the
EVALUATION -- what the gate says TODAY, against a session that exists on disk, with the
measured control values quoted rather than characterised.

Reads only `engine/out/mine/<session>/checkpoint.json`. Draws no surrogates, writes
nothing into any existing session directory, spends no holdout hash, logs no
EXPLORE_COUNT line. Nothing here is evidence.
"""

from __future__ import annotations

import json
import os

from . import observer

MINE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "mine")

# The Tranche B session §P7-22's D-5 line points at.
DEFAULT_SESSION = "session_20260813T092628"


def load_checkpoint(session=DEFAULT_SESSION, root=MINE_DIR):
    path = os.path.join(root, session, "checkpoint.json")
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def run(session=DEFAULT_SESSION, root=MINE_DIR):
    return observer.session_subdaily_gate(load_checkpoint(session, root), session)


def _print(rec):
    print("D-5 SUB-DAILY GATE -- session %s (config hash %s, %d scored tests)"
          % (rec["session"], rec["config_hash"], rec["n_scored_tests"]))
    print()
    p = rec["stage_1_presence"]
    print("  STAGE 1 (GATES) -- presence of the required F7 controls")
    for n in p["required"]:
        print("    %-28s %s" % (n, "PRESENT" if n in p["present"] else "MISSING"))
    print("    satisfied = %s" % p["satisfied"])
    print()
    print("  STAGE 2 (REPORTED, NOT GATING) -- what the controls measured")
    print("    %-28s %9s %10s %10s %10s  %s"
          % ("feature", "chi2", "p_shift", "p_boot", "p_raw", "at floor"))
    for r in rec["stage_2_reading"]["rows"]:
        def f(x, w=10, d=5):
            return ("%*.*f" % (w, d, x)) if isinstance(x, (int, float)) else "%*s" % (w, "-")
        print("    %-28s %9.3f %s %s %s  %s"
              % (r["feature"], r["chi2_score"] or 0.0, f(r["p_circular_shift"]),
                 f(r["p_block_bootstrap"]), f(r["p_raw"]),
                 ",".join(r["at_resolution_floor"]) or "-"))
    print("    %s" % rec["stage_2_reading"]["reading"])
    print()
    print("  VERDICT: %s" % rec["verdict"])
    print("  basis  : %s" % rec["verdict_basis"])
    if rec["count_path_degeneracy"]:
        print("  why    : %s" % rec["count_path_degeneracy"])
    print()
    print("  A PASS WOULD NOT LICENSE: %s" % rec["what_a_pass_would_not_license"])


if __name__ == "__main__":          # pragma: no cover - operator entry point
    import sys
    sess = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_SESSION
    rec = run(sess)
    _print(rec)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out",
                       "audit_subdaily_gate.json")
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2, default=float)
    print("\nwrote", out)
