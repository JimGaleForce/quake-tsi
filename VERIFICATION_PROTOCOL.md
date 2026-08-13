---
name: verified-claims
description: Use when any outgoing document (letter, paper, report, post, claim set) must be ACCURATE at expert-adversarial level in ONE build cycle plus ONE verification pass, instead of converging through many review rounds. Works for any subject or industry. Invoke before writing the first sentence, not after a draft exists.
---

# Verified Claims: expert-level accuracy in one cycle

Distilled 2026-08-13 from an empirical base: a two-letter scientific reply that took 15 blind
adversarial passes to converge, with every finding classified afterward. The method below reaches
the same end state in one build cycle + one mechanical harness run + one judgment pass, because it
removes the four causes that made 15 passes necessary. This is not a style guide; it is a
compiler discipline for factual documents.

## Why many passes never converge (the disease, measured)

1. Prose-first writing creates an unbounded claim surface; every fix round writes new unverified
   sentences (fixes are births). Measured: at least 4 of 15 rounds fixed errors introduced by
   earlier fixes.
2. Reviewers sample an open check-space instead of enumerating a closed one; each pass finds new
   things because "all checks" was never a finite set.
3. Reviewer findings get trusted and applied without independent verification. Reviewers inject
   errors at a real rate (measured: 3 confidently-asserted false findings in 15 passes, each
   costing a full round to detect, each caught only when a later pass COMPUTED instead of arguing).
4. Patch-editing (string replacement on prose) creates structural defects: duplicated sentences,
   corrupted text from regex boundaries, silent no-ops where an edit missed its target and the
   process reported success anyway.

## The method (five phases; only Phase 4 uses model judgment)

### Phase 0 - Fix the ground before building on it
Enumerate every source the document will rest on: counterpart statements (their emails, their
code, their paper), your own artifacts (files, data, published repo), external documents. Then
scan YOUR OWN sources for internal contradictions and staleness FIRST, and fix them before one
sentence of the document exists. A document cannot be righter than what it points at; mid-stream
source repair caused most of the churn in the empirical base (seven stale passages, one
inconsistent inventory, one wrong artifact field).

### Phase 1 - Compile the claim ledger (before any prose)
Every factual statement the document will make becomes a typed row:
claim | type | source pointer (file:field or doc:page, never "the paper") | check method | status.

Types and their non-negotiable rules:
- NUMBER: carries its convention in the row (units, zero-to-peak vs peak-to-peak, which
  denominator, which estimand, which window/region/population). A number without its convention
  is the single largest error class measured (it caused the retraction, the 6.6-vs-10.2
  conflation, the wrong-denominator percent, and a wrong-branch derivation "confirmed on four
  data points"). The check: a script extracts the number from the prose and recomputes or greps
  it from the pointed source.
- ATTRIBUTION ("X said/did/showed"): the row contains the VERBATIM quote or code lines. The claim
  may not exceed the quote. If the source does not state a property (a convention, a component, a
  mechanism), the claim says "which we read as..., a characterization that is ours" - ownership
  of readings is a load-bearing hedge, not politeness.
- OWN-CODE / OWN-DATA facts: verified by executing the code or parsing the data during Phase 1,
  never from memory. (Measured failure: "our port follows your construction" - false about our
  own code; "no annual result" - contradicted by our own repo.)
- MECHANISM ("this happens because..."): requires a derivation or an execution in the row.
  A plausible mechanism without one is rewritten as a QUESTION. (Measured: a proposed bias
  mechanism that cancels algebraically; a claimed library behavior that testing disproved.)
- POINTER: the preferred type. Point at the source table/figure instead of restating it.
  Every restatement is a claim; every claim is risk; the empirical base converged exactly as
  restatements were shed. If you cannot verify it now, point at it or cut it.

### Phase 2 - Render from the ledger, regenerate whole
Prose is written ONLY from ledger rows plus connective tissue. Every revision REGENERATES the
whole document from the ledger; never patch-edit prose with string replacement (the corruption /
duplication / silent-no-op class). If the document splits by audience (public letter + private
note), the ledger rows carry an audience tag, and the split is checked mechanically (Phase 3),
because audience leaks were among the most damaging findings.

### Phase 3 - Mechanical harness (scripted; minutes; blocks on failure)
Runs 100% of typed checks by construction:
- every number in the prose located, matched to its row, recomputed/grepped from its source;
- every attribution matched against its verbatim quote;
- every pointer resolved (file exists, is public if claimed public, says what the prose says);
- cross-reference integrity (no dangling "below"/"above"/"attached"), no duplicated sentences;
- audience isolation for split documents (no private-tagged content or its implications in the
  public part);
- style/encoding rules (this base: ASCII-only, no em-dashes - substitute your own).
CRITICAL: the harness's exit codes must GATE the pipeline. Two incidents in the base shipped
wrong states because a pipe or a compound command swallowed a failing check while the process
reported success. A verification whose result cannot stop the pipeline is decoration.

### Phase 4 - ONE adversarial judgment pass (the only model pass)
A single blind hostile reviewer, seeing only the final document + the sources (never the drafting
history), scoped to the untyped residue: tone, implicature, audience reading,
completeness-of-answer, framing, what-a-hostile-expert-infers. It carries the standing checklist
(below) so it does not rediscover known classes. Its binding rule - the one that ends the
injected-error cycle:

  A factual objection must carry its own derivation, execution, or verbatim citation.
  A finding without one is a QUESTION for the ledger, not a finding.

Reviewer findings that pass their own compute-or-quote bar become ledger rows (verified like any
row), the document regenerates, the harness re-runs. That is the loop's single iteration.

### Phase 5 - Stop rule (declared before starting, never moved after a result)
Ship when: harness 100% green AND the single judgment pass returns zero findings surviving its
own compute-or-quote rule. Do NOT run additional judgment passes "to be safe": each pass samples
new regions of an open space forever and injects errors at the measured rate. Confidence comes
from closed coverage of the typed surface plus one disciplined read of the residue - not from
exhaustion.

## The standing checklist (grows; every new finding class joins it permanently)
- Conventions travel with numbers (units, peak convention, denominator, estimand, window,
  population, alpha/power for statistics).
- A back-solved or fitted parameter is resolved by its PHYSICAL definition, never by fit alone;
  where two branches fit equally, the fit is evidence for neither, and a datum that fits only
  under a rationalization is a disconfirmation.
- Claims about counterparts limited to their verbatim words and public artifacts; their private
  material never described to third parties; split-audience content checked mechanically.
- Praise/grading of senior counterparts, priority claims, and "their field needed this" framings:
  flag-and-justify or cut.
- Level observations on oscillatory data need dwell-time nulls (what does pure chance look like
  on this waveform?) before "clustering" is claimed.
- Small-n strong-concentration is decisive; few events is not weak effect - do the binomial
  before dismissing.
- Both-sided disclosure: when quoting a family of measurements, quote the unflattering members
  or point at the full table.
- Every document that invites readers into a repository gets a repository self-consistency scan
  (stale text contradicting the document is a finding against the document).
- The fix for a violated invariant conditions on the invariant, not the example.

## Economics (measured on the base case)
15 judgment passes at ~120k tokens each, plus ~8 fix rounds: roughly 2M tokens, a full day, and
three injected-error detours. The method: one build worker (Phases 0-2), one harness run
(scripted, near-free), one judgment pass, one regeneration - roughly 5-7x cheaper and a day
faster, with confidence that is CALIBRATED (closed coverage + disciplined residue) rather than
exhausted into being.

## How to invoke on any subject
Tell the persona: "Produce X under the verified-claims discipline." The persona must (a) fix its
own sources first, (b) compile the typed ledger before writing, (c) render from the ledger and
regenerate whole on every change, (d) run the mechanical harness with gating exit codes,
(e) commission exactly one blind compute-or-quote adversarial read, (f) ship on the declared stop
rule. Domain changes the sources and the checklist entries; the discipline does not change.
