// D-1 NODE BRIDGE (§P7-23(C)) -- evaluates the FROZEN scalar with the observing
// application's own code, unmodified, per K092_SCALAR_PROVENANCE.md §5(3):
//
//   "score D-12/D-13 with the app's code verbatim (the Node import used in §3 above
//    runs astro.ts unmodified and is reproducible in one line), and report the
//    sitetide.py classification alongside as a robustness check."
//
// This file contains NO tidal arithmetic of its own. It imports
// earth-tides-globe/src/utils/astro.ts and calls solidTideDisplacementCm /
// solidTideRateCmPerHour. Any change to the app changes these numbers, which is the
// point: the frozen scalar IS the app's scalar.
//
// Usage:  node exp_k092_d1_bridge.mjs <in.json> <out.json>
// in.json : {"astro_url":..., "half_window_days":..,"step_minutes":..,
//            "events":[{"id","t_ms","lat","lon"},...]}
// out.json: per event -> {level_cm, rate_cm_per_h, grid:{t_ms0, step_ms, n, level[]}}
//           plus the self-test reproducing the frozen 1991 readout.

import { readFileSync, writeFileSync } from 'node:fs';

const [, , inPath, outPath] = process.argv;
if (!inPath || !outPath) { console.error('usage: node exp_k092_d1_bridge.mjs <in.json> <out.json>'); process.exit(2); }

const cfg = JSON.parse(readFileSync(inPath, 'utf8'));
const astro = await import(cfg.astro_url);
const { solidTideDisplacementCm, solidTideRateCmPerHour } = astro;
if (typeof solidTideDisplacementCm !== 'function' || typeof solidTideRateCmPerHour !== 'function') {
  console.error('astro.ts did not export the two frozen functions'); process.exit(3);
}

// --- provenance self-test: the frozen 1991-05-30 Sand Point readout, exactly. ---
// K092_SCALAR_PROVENANCE.md §3. If this does not match to the digit, the app's source
// has moved since the provenance was closed and D-1 must NOT be scored against it.
const T91 = Date.UTC(1991, 4, 30, 13, 17, 41);
const selftest = {
  displacement_cm: solidTideDisplacementCm(T91, 54.57, -161.61),
  rate_cm_per_h: solidTideRateCmPerHour(T91, 54.57, -161.61),
  expected_displacement_cm: -13.206888677138735,
  expected_rate_cm_per_h: -0.5757627956078135,
};
selftest.exact_match =
  selftest.displacement_cm === selftest.expected_displacement_cm &&
  selftest.rate_cm_per_h === selftest.expected_rate_cm_per_h;
if (!selftest.exact_match) {
  writeFileSync(outPath, JSON.stringify({ selftest, aborted: true }, null, 2));
  console.error('PROVENANCE SELF-TEST FAILED -- app source has moved; refusing to score.');
  process.exit(4);
}

const halfDays = cfg.half_window_days;
const stepMs = Math.round(cfg.step_minutes * 60000);
const nHalf = Math.round((halfDays * 86400000) / stepMs);
const n = 2 * nHalf + 1;

const out = [];
for (const ev of cfg.events) {
  const t0 = ev.t_ms - nHalf * stepMs;
  const level = new Array(n);
  for (let i = 0; i < n; i++) level[i] = solidTideDisplacementCm(t0 + i * stepMs, ev.lat, ev.lon);
  out.push({
    id: ev.id,
    level_cm: solidTideDisplacementCm(ev.t_ms, ev.lat, ev.lon),
    rate_cm_per_h: solidTideRateCmPerHour(ev.t_ms, ev.lat, ev.lon),
    grid: { t_ms0: t0, step_ms: stepMs, n, level },
  });
}

writeFileSync(outPath, JSON.stringify({
  selftest,
  scalar: 'earth-tides-globe/src/utils/astro.ts::solidTideDisplacementCm (cm, UP POSITIVE)',
  rate: 'solidTideRateCmPerHour (+-10 min central difference, cm/h)',
  node_version: process.version,
  half_window_days: halfDays, step_minutes: cfg.step_minutes,
  events: out,
}));
console.error(`bridge OK: ${out.length} events, grid n=${n} each, selftest exact.`);
