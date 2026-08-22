// D-1 WAVEFORM-MATCHED NULL BRIDGE.
//
// Same contract as exp_k092_d1_bridge.mjs -- imports the observing application's
// astro.ts unmodified and contains NO tidal arithmetic and NO analysis of its own --
// but emits a LONG grid per site as raw little-endian Float64 so Python can run the
// IDENTICAL local-cycle analysis over it. Emitting the grid rather than a summary is
// deliberate: computing the band fractions in JavaScript would create a second
// implementation of the analysis, and then the null and the observation would be
// compared across two code paths instead of two data sets.
//
// Usage: node exp_k092_d1_null_bridge.mjs <in.json> <out.bin> <out.json>

import { readFileSync, writeFileSync } from 'node:fs';

const [, , inPath, binPath, metaPath] = process.argv;
if (!inPath || !binPath || !metaPath) {
  console.error('usage: node exp_k092_d1_null_bridge.mjs <in.json> <out.bin> <out.json>');
  process.exit(2);
}

const cfg = JSON.parse(readFileSync(inPath, 'utf8'));
const astro = await import(cfg.astro_url);
const { solidTideDisplacementCm } = astro;
if (typeof solidTideDisplacementCm !== 'function') {
  console.error('astro.ts did not export solidTideDisplacementCm'); process.exit(3);
}

// Provenance self-test, identical to the observation bridge's.
const T91 = Date.UTC(1991, 4, 30, 13, 17, 41);
const disp = solidTideDisplacementCm(T91, 54.57, -161.61);
if (disp !== -13.206888677138735) {
  console.error('PROVENANCE SELF-TEST FAILED -- app source has moved; refusing to run.');
  process.exit(4);
}

const stepMs = Math.round(cfg.step_minutes * 60000);
const nHalf = Math.round((cfg.half_span_days * 86400000) / stepMs);
const n = 2 * nHalf + 1;

const buf = Buffer.alloc(cfg.sites.length * n * 8);
let off = 0;
for (const s of cfg.sites) {
  const t0 = s.t_ms - nHalf * stepMs;
  for (let i = 0; i < n; i++) {
    buf.writeDoubleLE(solidTideDisplacementCm(t0 + i * stepMs, s.lat, s.lon), off);
    off += 8;
  }
}
writeFileSync(binPath, buf);
writeFileSync(metaPath, JSON.stringify({
  selftest_displacement_cm: disp,
  selftest_exact: true,
  scalar: 'earth-tides-globe/src/utils/astro.ts::solidTideDisplacementCm (cm, UP POSITIVE)',
  node_version: process.version,
  n_sites: cfg.sites.length,
  n_per_site: n,
  step_minutes: cfg.step_minutes,
  half_span_days: cfg.half_span_days,
  order: cfg.sites.map((s) => s.id),
  dtype: 'float64-le',
}, null, 2));
console.error(`null bridge OK: ${cfg.sites.length} sites x ${n} samples`);
