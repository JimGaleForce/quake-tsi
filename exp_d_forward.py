"""EXP-D - forward favorable-phase windows, per OVERNIGHT_PREDICTION_PROTOCOL.md.

The protocol permits a forward statement ONLY for bins passing EXP-A's frozen success rule
(pooled test p < 0.05 AND clean control). EXP-A failed its success rule (pooled S = -0.100,
p = 0.942; the single train-selected bin anti-aligned out of sample), so the required output
is the null statement below. This script exists so the pipeline's answer is explicit and
reproducible rather than an absence.
"""
import json
from pathlib import Path

HERE = Path(__file__).parent
a = json.loads((HERE / "results_exp_a.json").read_text(encoding="utf-8"))

passed = a["test"]["p"] < 0.05 and a["control"]["p"] > 0.05
out = {
    "protocol": "OVERNIGHT_PREDICTION_PROTOCOL.md#EXP-D",
    "exp_a_success_rule_passed": bool(passed),
    "exp_a_test_p": a["test"]["p"],
    "exp_a_pooled_S": a["test"]["pooled_S"],
    "forward_windows": [],
    "statement": (
        "No validated basis for a forward statement: EXP-A's out-of-sample phase-skill test "
        "failed its pre-registered success rule (train-selected bin anti-aligned in the test "
        "window; anti-leak control clean). Under the frozen protocol, no favorable-phase "
        "windows are issued." if not passed else "UNEXPECTED: rule passed; implement windows."
    ),
}
(HERE / "results_exp_d.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(out["statement"])
