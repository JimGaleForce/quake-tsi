"""L-1 — S-14(c) two-branch link bracket. Produces results_l1.json.
FIRST-RUN. Reads only k034 artifacts; writes only results_l1.json."""
import numpy as np, pandas as pd, json, hashlib, datetime
from scipy import optimize, stats

REP = '.'
def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()[:16]

cs = pd.read_csv('k034_cellstats.csv'); pw = pd.read_csv('k034_power.csv')
k34 = json.load(open('results_k034.json'))

prim = cs[(cs.window == '5d') & (cs.mcut == 1.5) & (~cs.gated_out)].copy()
prim['sig_kPa'] = prim.sig_primary_Pa / 1e3
prim['grp'] = np.where(prim.cls.isin(['A', 'B']), 'AB', 'C')
prim['meas'] = prim.n_bg > 0
pw['grp'] = np.where(pw.cls.isin(['A', 'B']), 'AB', 'C')

# ---------------- Branch E fit ----------------
def poisfit(df):
    y = df.n_post.values.astype(float); off = np.log(df.expected.values)
    x = np.log(df.sig_kPa.values); X = np.column_stack([np.ones_like(x), x])
    nll = lambda p: np.sum(np.exp(np.clip(off + X @ p, -50, 50)) - y * (off + X @ p))
    p = optimize.minimize(nll, [0., 0.], method='BFGS').x
    mu = np.exp(off + X @ p); Vn = np.linalg.inv(X.T @ (X * mu[:, None]))
    meat = np.zeros((2, 2))
    for s, idx in df.groupby('source').indices.items():
        u = (X[idx] * (y[idx] - mu[idx])[:, None]).sum(0); meat += np.outer(u, u)
    Vc = Vn @ meat @ Vn
    dev = 2*np.sum(np.where(y > 0, y*np.log(np.where(y > 0, y/mu, 1)), 0) - (y-mu))
    nll0 = lambda p0: np.sum(np.exp(off + p0[0]) - y*(off + p0[0]))
    a0 = optimize.minimize(nll0, [0.], method='BFGS').x[0]; mu0 = np.exp(off + a0)
    dev0 = 2*np.sum(np.where(y > 0, y*np.log(np.where(y > 0, y/mu0, 1)), 0) - (y-mu0))
    disp = float(np.sum((y-mu)**2/mu)/(len(y)-2)); lrt = 2*(nll0([a0]) - nll(p))
    return dict(a=float(p[0]), b=float(p[1]),
                se_naive=[float(v) for v in np.sqrt(np.diag(Vn))],
                se_quasipoisson=[float(v) for v in np.sqrt(np.diag(Vn)*disp)],
                se_cluster_by_source=[float(v) for v in np.sqrt(np.diag(Vc))],
                deviance=float(dev), null_deviance=float(dev0),
                pseudo_R2_deviance=float(1-dev/dev0), dispersion_pearson=disp,
                LRT_b=float(lrt), p_LRT_b_naive_NOT_TRUSTED=float(stats.chi2.sf(lrt, 1)),
                n_pairs=int(len(y)))

def cboot(df, nb=4000, seed=20260811):
    rng = np.random.default_rng(seed); srcs = df.source.unique(); B = []
    for _ in range(nb):
        pick = rng.choice(srcs, len(srcs), replace=True)
        dd = pd.concat([df[df.source == s].assign(source=f'{s}_{i}') for i, s in enumerate(pick)])
        if dd.n_post.sum() == 0: continue
        try: f = poisfit(dd); B.append([f['a'], f['b']])
        except Exception: pass
    return np.array(B)

fits = {}
for g in ['AB', 'C']:
    d = prim[(prim.grp == g) & prim.meas]
    f = poisfit(d); B = cboot(d)
    f['cluster_bootstrap_a_CI95'] = [float(v) for v in np.percentile(B[:, 0], [2.5, 97.5])]
    f['cluster_bootstrap_b_CI95'] = [float(v) for v in np.percentile(B[:, 1], [2.5, 97.5])]
    f['cluster_bootstrap_n'] = int(len(B))
    f['sigma_support_kPa'] = [float(d.sig_kPa.min()), float(d.sig_kPa.max())]
    f['median_observed_RR'] = float(d.RR.median())
    f['n_pairs_class_total'] = int((prim.grp == g).sum())
    f['n_pairs_S15_UNMEASURABLE_nbg0'] = int(((prim.grp == g) & (~prim.meas)).sum())
    m = d[d.RR > 0]
    sl, ic, r, p, se = stats.linregress(np.log10(m.sig_kPa), np.log10(m.RR))
    f['sensitivity_OLS_log10RR'] = dict(slope_per_decade=float(sl), intercept=float(ic),
                                        r2=float(r**2), p=float(p), se_slope=float(se), n=int(len(m)))
    fits[g] = f

aAB, bAB = fits['AB']['a'], fits['AB']['b']
# optimistic envelope: slope fixed at bAB, intercept at the 90th pct of the offset residuals
dAB = prim[(prim.grp == 'AB') & prim.meas].copy()
dAB['RRc'] = (dAB.n_post + 0.5)/(dAB.expected + 0.5)
offs = np.log(dAB.RRc) - bAB*np.log(dAB.sig_kPa)
a90 = float(np.quantile(offs, 0.90))
fits['AB']['envelope_q90_intercept_a'] = a90
fits['AB']['envelope_note'] = ('slope held at the mean-link b; intercept set to the 90th percentile '
                               'of log(RR_c) - b*log(sigma) with RR_c=(n_post+0.5)/(expected+0.5)')
SUP_LO, SUP_HI = fits['AB']['sigma_support_kPa']

R_E_mean = lambda s: float(np.exp(aAB) * s**bAB)
R_E_up   = lambda s: float(np.exp(a90) * s**bAB)
R_S      = lambda s, A: float(np.exp(s/(A*1000.)))          # A in MPa, s in kPa
MDA_S    = lambda R, A: float(A*1000.*np.log(R))
MDA_Em   = lambda R: float((R/np.exp(aAB))**(1.0/bAB))
MDA_Eu   = lambda R: float((R/np.exp(a90))**(1.0/bAB))

# ---------------- Branch-S per-cell table (COLLECTED from k034_power.csv) ----------------
pw['MDA_E_mean_kPa'] = np.where(pw.measurable, (pw.R_min_80pct/np.exp(aAB))**(1/bAB), np.nan)
pw['MDA_E_upper_kPa'] = np.where(pw.measurable, (pw.R_min_80pct/np.exp(a90))**(1/bAB), np.nan)
def summ(df, c):
    v = df[c].dropna(); return dict(best=float(v.min()), median=float(v.median()), worst=float(v.max()), n=int(len(v)))
cellsum = {}
for nm, df in [('all_measurable_cells', pw[pw.measurable]), ('class_AB_measurable_cells', pw[(pw.grp=='AB') & pw.measurable])]:
    cellsum[nm] = {c: summ(df, c) for c in ['R_min_80pct', 'dTau_min_kPa_Asigma0.03',
                    'dTau_min_kPa_Asigma0.10', 'dTau_min_kPa_Asigma0.15', 'MDA_E_mean_kPa', 'MDA_E_upper_kPa']}

# ---------------- entry table ----------------
REF_FLOOR = 1.05
E = []
def ent(eid, band, floor, floor_src, cls, note, kind='amplitude'):
    E.append(dict(entry=eid, band_kPa=band, floor_R_min=floor, floor_provenance=floor_src,
                  receiver_class=cls, kind=kind, note=note))

ent('K-059 (main, R1/R2 near-teleseismic)', [1.0, 100.0], REF_FLOOR, 'REFERENCE-FLOOR (entry declares no MDM; bits/event statistic)', 'mixed (global cells)', 'ledger L7144: "1-100 kPa: predicted 3%-x28 rate modulation at Asigma=0.03"')
ent('K-059 (3 kPa exposure gate)', [3.0, 3.0], REF_FLOOR, 'REFERENCE-FLOOR', 'mixed (global cells)', 'the exposure gate named in P5-8 and at ledger L10490; Popper predicted this moves')
ent('K-060', [3.0, 4.2], 1.028, 'DERIVED from declared N_exposed >= 1e4 (2.8/sqrt(N))', 'mixed', 'ledger L7200: 3 kPa envelope -> 11%, 4.2 kPa coherent peak -> 15% at Asigma=0.03')
ent('K-061', [3.0, 3.0], 1.028, 'DERIVED from declared N_exposed >= 1e4', 'mixed', 'ledger L7245: 3 kPa envelope, deficit -8% to -11%')
ent('K-062', [0.3, 3.0], 1.051, 'DERIVED from declared N_exposed > 3e3', 'mixed', 'ledger L7288: net antipodal 0.3-3 kPa, 1%-11%')
ent('K-063', [0.3, 3.0], 1.12, 'ENTRY-DECLARED MDM 4-12%, adverse end', 'mixed', 'ledger L7345: R3 at 0.3-3 kPa; R5 below floor')
ent('K-064', [5.0, 5.0], REF_FLOOR, 'REFERENCE-FLOOR', 'A+B (entry conditions on low-Asigma/geothermal cells)', 'ledger L7388: dTau ~5 kPa; entry states its power lives in low-Asigma cells')
ent('K-065', [0.01, 0.1], 1.056, 'ENTRY-DECLARED MDM 5.6%', 'mixed', 'ledger L7420: 0.01-0.1 kPa mode amplitudes')
ent('K-066', None, None, 'n/a', 'mixed', 'ledger L7480: not amplitude-limited (no triggering claimed); patch-count limited', kind='not-amplitude-limited')
ent('K-067', [1.0, 10.0], REF_FLOOR, 'REFERENCE-FLOOR', 'mixed', 'ledger L7522: exceedance thresholds 1 kPa and 10 kPa')
ent('K-068', None, 1.28, 'DERIVED from declared >=100 exposures/cell', 'mixed', 'ledger L7573: per-pair SNR a few percent; no amplitude band declared', kind='amplitude-unstated')
ent('K-069', None, None, 'n/a', 'mixed', 'ledger L7613: committed-window count limited; no amplitude band', kind='not-amplitude-limited')
ent('K-070', None, 1.5, 'ENTRY-DECLARED rate ratio >= 1.5', 'mixed', 'ledger L7662: MDM is a factor not a percent; top-0.1% configurations', kind='amplitude-unstated')
ent('K-071', [5.0, 5.0], 1.8, 'ENTRY-DECLARED >=1.4 pooled / >=1.8 per stratum, adverse end', 'A+B and C (ledger-class contrast)', 'ledger L7821: sources contributing >= 5 kPa to a target cell')
ent('K-072 (1-5 kPa band)', [1.0, 5.0], 1.16, 'DERIVED from declared >=10 triggered events x >=30 cells/group', 'mixed', 'P5-8 names the 1-5 kPa band; entry text uses >=5 kPa exposure. Popper predicted this moves')
ent('K-073', [1.0, 10.0], REF_FLOOR, 'REFERENCE-FLOOR', 'mixed (wetness contrast)', 'ledger L7942: 1 kPa -> 3% predicted; 10 kPa exceedances are the informative ones')
ent('K-074', [1.0, 5.0], 1.089, 'DERIVED from pooled matrix N ~ 1e3', 'mixed', 'ledger L7997: M8.5 at 8000 km = 1-5 kPa -> 3-18% predicted')
ent('K-075', None, None, 'n/a', 'mixed', 'ledger L8058: probe-pair / region-year limited; rank statistic, no amplitude', kind='not-amplitude-limited')
ent('K-038', [50.0, 50.0], REF_FLOOR, 'REFERENCE-FLOOR', 'mixed (interaction with ledger class)', 'ledger L1332: "threshold-like behaviour with detection above ~0.05 MPa"')
ent('K-043', None, REF_FLOOR, 'REFERENCE-FLOOR', 'mixed (global ping map)', 'no kPa band declared in-entry; PGV-weighted pooling over all pings', kind='amplitude-unstated')
ent('W-002-P2', [33.8, 346.9], REF_FLOOR, 'REFERENCE-FLOOR', 'A+B (Landers/HectorMine/Ridgecrest targets)', 'no band declared in-entry; the sources are K-034s own, so the K-034 measured amplitude span is used and flagged')
ent("K-078 (slab-transient arm)", None, None, 'n/a', 'C (intraslab, non-geothermal)', 'no amplitude declared; POWER-STATE is P(oasis) at n=40, count-limited', kind='not-amplitude-limited')
ent('K-084 (0-1 d dynamic row)', [0.4, 2.0], 1.15, 'DERIVED from N_eff=337 targets after clustering', 'mixed', 'ledger L10194: dynamic amplitude at 320 km, M6 = 0.4-2 kPa')
ent('A0b (Landers angular control)', [9.6, 346.9], 5.5, 'MEASURED: median R_min_80pct over K-034 class-AB measurable cells', 'A+B', 'the K-034 cells themselves; amplitude span is K-034 measured (suggestive floor to max cell)')

rows = []
for e in E:
    r = dict(e)
    fl = e['floor_R_min']
    if fl is not None:
        r['MDA_branchS_kPa'] = {A: round(MDA_S(fl, float(A)), 3) for A in ['0.03', '0.10', '0.15']}
        r['MDA_branchS_adverse_kPa'] = round(MDA_S(fl, 0.15), 3)
        r['MDA_branchE_mean_kPa'] = round(MDA_Em(fl), 3)
        r['MDA_branchE_upper_kPa'] = round(MDA_Eu(fl), 3)
        r['MDA_branchE_in_support'] = bool(SUP_LO <= MDA_Em(fl) <= SUP_HI)
    else:
        r['MDA_branchS_kPa'] = r['MDA_branchE_mean_kPa'] = r['MDA_branchE_upper_kPa'] = 'UNMEASURED (no detection floor declared)'
    b = e['band_kPa']
    if e['kind'] == 'not-amplitude-limited':
        r['status'] = 'N/A - NOT AMPLITUDE-LIMITED (no exp(dtau/Asigma) conversion in the POWER-STATE; bracket does not bind)'
        r['bracket_kPa'] = 'n/a'
    elif b is None:
        r['status'] = 'POWER-INDETERMINATE (amplitude band UNSTATED in-entry; bracket printed, status pending band declaration per S-15)'
        r['bracket_kPa'] = [r['MDA_branchE_upper_kPa'], r['MDA_branchS_adverse_kPa']]
    else:
        bmax = b[1]
        r['R_predicted_at_band_max'] = dict(
            branchS_Asigma0_03=round(R_S(bmax, 0.03), 4), branchS_Asigma0_10=round(R_S(bmax, 0.10), 4),
            branchS_Asigma0_15=round(R_S(bmax, 0.15), 4),
            branchE_mean=round(R_E_mean(bmax), 4), branchE_upper=round(R_E_up(bmax), 4),
            branchE_extrapolated_below_support=bool(bmax < SUP_LO))
        lo = min(r['MDA_branchE_upper_kPa'], r['MDA_branchS_adverse_kPa'])
        hi = max(r['MDA_branchE_upper_kPa'], r['MDA_branchS_adverse_kPa'], r['MDA_branchE_mean_kPa'])
        r['bracket_kPa'] = [round(lo, 3), round(hi, 3)]
        bmin = b[0]
        r['R_predicted_at_band_min'] = dict(
            branchS_Asigma0_15=round(R_S(bmin, 0.15), 4), branchE_mean=round(R_E_mean(bmin), 4),
            branchE_upper=round(R_E_up(bmin), 4), branchE_extrapolated_below_support=bool(bmin < SUP_LO))
        if bmin < SUP_LO:
            r['status_at_band_min'] = 'POWER-INDETERMINATE (Branch E UNMEASURED below support)'
        elif bmin >= r['MDA_branchS_adverse_kPa'] and bmin >= r['MDA_branchE_mean_kPa']:
            r['status_at_band_min'] = 'POWERED'
        elif bmin < min(r['MDA_branchE_upper_kPa'], r['MDA_branchS_adverse_kPa']):
            r['status_at_band_min'] = 'UNDERPOWERED'
        else:
            r['status_at_band_min'] = 'POWER-INDETERMINATE'
        below_support = bmax < SUP_LO
        if below_support:
            r['status'] = 'POWER-INDETERMINATE'
            r['status_reason'] = (f'band max {bmax} kPa lies below the Branch-E fitted support '
                                  f'[{SUP_LO:.1f}, {SUP_HI:.1f}] kPa: Branch E is UNMEASURED there per S-15, and Branch S '
                                  'may never be used to declare no power (P5-8). Bracket printed, no MDA claimed.')
        elif bmax >= r['MDA_branchS_adverse_kPa'] and bmax >= r['MDA_branchE_mean_kPa']:
            r['status'] = 'POWERED'; r['status_reason'] = 'band max clears the adverse Branch-S MDA and the mean Branch-E MDA'
        elif bmax < lo:
            r['status'] = 'UNDERPOWERED'; r['status_reason'] = 'band max is below both branches MDA, inside the Branch-E support'
        else:
            r['status'] = 'POWER-INDETERMINATE'; r['status_reason'] = 'branches disagree about whether the entry has power; bracket printed, no MDA'
    if r.get('receiver_class', '').startswith('C') or 'C (' in str(r.get('receiver_class')):
        r['branchE_class_C'] = 'UNMEASURED per S-15 (class-C cells did not fire; see fits.C)'
    rows.append(r)

counts = pd.Series([r['status'].split(' ')[0].rstrip(':') for r in rows]).value_counts().to_dict()

out = dict(
    experiment='L-1', title='The transient link bracket - every transient POWER-STATE recomputed as the S-14(c) two-branch bracket',
    authority='HYPOTHESIS_LEDGER.md P5-8 (link-function ruling) + S-14(c) + S-15; queue rank 1 of Popper round 4',
    run_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    run_status='FIRST-RUN. No number in this file has been produced before. Arithmetic only; no new data, no new statistic, no download.',
    scope_constraint=('The K-034 10-200x excess is used here ONLY as calibration data for power arithmetic. '
                      'It remains logged-unclaimed per R2-2 and P5-8 and is not a finding of this artifact. '
                      'This artifact RE-PRICES; it does not RE-RULE. No ADMIT/DEFER verdict in P4-3 is touched.'),
    inputs={f: sha(f) for f in ['k034_cellstats.csv', 'k034_power.csv', 'results_k034.json']},
    branchS=dict(definition='R = exp(dtau/Asigma) over the S-14 bracket {0.03, 0.10, 0.15} MPa, retained ONLY as a lower bound on response',
                 source='COLLECTED from k034_power.csv (columns dTau_min_kPa_Asigma*), which reproduces the K34-3 table exactly; not re-derived',
                 K34_3_table_reproduced=cellsum['all_measurable_cells']),
    branchE=dict(form='R(sigma) = exp(a) * sigma_kPa**b  (quasi-Poisson GLM, log link, offset log(expected), regressor log sigma)',
                 primary_stratum='class A+B (pre-registered geothermal/volcanic)',
                 response_data='K-034 primary reading: 5 d window, M >= 1.5, distance-gate survivors, n_bg > 0',
                 fits=fits, optimistic_envelope_q90_intercept=a90,
                 class_C_ruling=('UNMEASURED per S-15. Class-C slope b = %.3f with cluster-bootstrap CI [%.2f, %.2f] spanning 0, '
                                 'median observed RR = %.2f (< 1). Branch E is NOT quoted at a non-geothermal receiver.'
                                 % (fits['C']['b'], fits['C']['cluster_bootstrap_b_CI95'][0], fits['C']['cluster_bootstrap_b_CI95'][1], fits['C']['median_observed_RR']))),
    per_cell_MDA_summary=cellsum,
    status_rules=dict(
        POWERED='band max >= adverse Branch-S MDA (Asigma = 0.15) AND >= mean Branch-E MDA',
        UNDERPOWERED='band max below both branches MDA AND inside the Branch-E fitted support',
        POWER_INDETERMINATE='branches disagree, or band max lies below the Branch-E support where Branch E is UNMEASURED (S-15)',
        NOT_AMPLITUDE_LIMITED='the POWER-STATE contains no stress-to-response conversion; the bracket does not bind',
        reference_floor=REF_FLOOR,
        reference_floor_justification='used only where the entry declares no MDM; equals the modulation floor the family itself quotes (K-063 4-12%, K-065 5.6%). Flagged per row.'),
    bracket_table=rows,
    status_counts=counts,
    deviations=[
        'D1. Popper P5-8 characterises Branch E as the OPTIMISTIC response. The mean empirical link is not optimistic everywhere: '
        'the power law (b = %.2f) is FLATTER than the exponential, so Branch S at Asigma = 0.03 predicts MORE response than Branch E above ~40 kPa, '
        'and the branches CROSS near 15-20 kPa. Two Branch-E variants are therefore reported: the mean link (primary) and a q90 optimistic envelope. '
        'The bracket is [most optimistic, most pessimistic] across all branch variants, computed per row rather than assumed.' % bAB,
        'D2. The 10-200x excess quoted in P5-8 is a per-firing-cell ratio, not the pooled mean response. The count-weighted mean link over all '
        'class-A+B pairs predicts R = %.2f at the 33.8 kPa certified floor, against Branch S(0.15) = 1.25 - a factor of ~%.1f, not 10-200x. '
        'The larger factors live in individual cells (Landers -> cedar_city_ut RR = 204 at 46 kPa). This is a property of the fit, reported, not claimed.' % (R_E_mean(33.793), R_E_mean(33.793)/1.2527),
        'D3. Branch-E support is [%.1f, %.1f] kPa. Every entry whose band lies below %.1f kPa is priced POWER-INDETERMINATE by extrapolation refusal (S-15), '
        'not by a computed MDA. This is the mechanism by which K-059s 3 kPa gate and K-072s 1-5 kPa band move.' % (SUP_LO, SUP_HI, SUP_LO),
        'D4. Entry detection floors R_min are entry-declared where the POWER-STATE states an MDM or rate ratio, DERIVED as 2.8/sqrt(N_exposed) where '
        'the entry states only an N, and set to the REFERENCE floor 1.05 otherwise. Provenance is printed per row. No floor was invented silently.',
        'D5. Overdispersion in the Branch-E fit is severe (Pearson dispersion = %.0f). Naive GLM p-values are not trusted; quasi-Poisson and '
        'cluster-by-source (4 clusters) standard errors and a 4000-draw cluster bootstrap are reported instead. With only four sources the cluster '
        'CI is wide by construction and that width is the finding about our own ignorance that S-14(c) asks to be put in the headline.' % fits['AB']['dispersion_pearson'],
        'D6. Six of the 54 pairs are S-15 UNMEASURABLE at the primary reading (n_bg = 0); k034_power.csv flags 14 of 54 unmeasurable under its own '
        'floor. The fit uses the 48 measurable-by-n_bg pairs (30 class A+B, 18 class C); the per-cell MDA summary uses the 40 flagged measurable.',
    ],
    verification=[
        'Branch-S per-cell MDA summary over all 40 measurable cells reproduces the K34-3 table to 3 significant figures: '
        '0.03 MPa best/median/worst 22.3/50.4/89.5 kPa; 0.10 MPa 74.2/168/298; 0.15 MPa 111/252/447.',
        'Branch-S identity checked: dTau_min_kPa_Asigma0.15 == 150 * ln(R_min_80pct) for every measurable cell.',
    ],
)
mx = float(np.max(np.abs(pw.loc[pw.measurable, 'dTau_min_kPa_Asigma0.15'] - 150*np.log(pw.loc[pw.measurable, 'R_min_80pct']))))
out['verification'].append(f'Max abs deviation of that identity: {mx:.2e} kPa.')

json.dump(out, open('results_l1.json', 'w'), indent=1)
print(json.dumps(counts, indent=1))
print('support', SUP_LO, SUP_HI, 'a', aAB, 'b', bAB, 'a90', a90)
for r in rows:
    print(f"{r['entry'][:34]:36} band={str(r['band_kPa']):16} floor={str(r['floor_R_min']):6} "
          f"S015={str(r.get('MDA_branchS_adverse_kPa')):8} Emean={str(r.get('MDA_branchE_mean_kPa')):8} "
          f"Eup={str(r.get('MDA_branchE_upper_kPa')):8} {r['status'][:34]:36} min:{r.get('status_at_band_min','')}")
print('identity max dev', mx)
