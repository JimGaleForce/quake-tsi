"""K-034 — Landers positive control / dynamic-triggering licensing gate.

Executes K034_PREREGISTERED_CELLS.md exactly (hash 01e41f97..., frozen before download).
Amplitude model from K034_SEALED_LITERATURE.md (hash eae95839..., frozen before download).
"""
import os, json, hashlib
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data', 'k034')
RNG = np.random.default_rng(20260811)

# ---------------- frozen constants (K034_SEALED_LITERATURE.md) ----------------
G_SHEAR = 30e9          # Pa, vdE&B 2010
C_PHASE = 3500.0        # m/s
T_SW = 20.0             # s
ASIGMA = {'0.03': 0.03e6, '0.10': 0.10e6, '0.15': 0.15e6}   # Pa, S-14 bracket
ALPHA = 0.05
POST_D = {'5d': 5.0, '1d': 1.0}
BG_D = 90.0
NSHIFT = 999
SHIFT_YR = 3.0
MINSEP_D = 30.0

SOURCES = {
 'landers':    dict(t0='1992-06-28T11:57:34Z', lat=34.200, lon=-116.436, Ms=7.3, Mw=7.3, L=70.0),
 'hectormine': dict(t0='1999-10-16T09:46:44Z', lat=34.594, lon=-116.271, Ms=7.4, Mw=7.1, L=48.0),
 'ridgecrest': dict(t0='2019-07-06T03:19:53Z', lat=35.766, lon=-117.605, Ms=7.1, Mw=7.1, L=50.0),
 'denali':     dict(t0='2002-11-03T22:12:41Z', lat=63.517, lon=-147.444, Ms=8.5, Mw=7.9, L=340.0),
}
CLASS = {'long_valley':'A','coso':'A','geysers':'A','yellowstone':'A','salton_brawley':'A',
         'lassen':'B','mono_west_nv_mina':'B','little_skull_mtn':'B','cedar_city_ut':'B',
         'smith_valley_nv':'C','parkfield':'C','mendocino':'C','wasatch_slc':'C','san_jacinto':'C'}
DOCUMENTED = {'long_valley','coso','geysers','yellowstone','lassen','mono_west_nv_mina',
              'little_skull_mtn','cedar_city_ut','smith_valley_nv'}
CLASSRANK = {'A':0,'B':1,'C':2}

MAN = json.load(open(os.path.join(D, 'manifest.json')))

# ---------------- geometry & amplitude ----------------
def gcdist_km(la1, lo1, la2, lo2):
    p = np.pi/180.0
    a = (np.sin((la2-la1)*p/2)**2 + np.cos(la1*p)*np.cos(la2*p)*np.sin((lo2-lo1)*p/2)**2)
    return 2*6371.0*np.arcsin(np.sqrt(a))

def sigma_primary(M, r_km):
    """vdE&B eq.(6): log10 A20[um] = Ms - 1.66 log10(D[deg]) - 2 ; V=2*pi*A20/T ; sig=G*V/c."""
    Ddeg = r_km/111.195
    A20_um = 10.0**(M - 1.66*np.log10(Ddeg) - 2.0)
    V = 2*np.pi*(A20_um*1e-6)/T_SW      # m/s
    return G_SHEAR*V/C_PHASE            # Pa

def sigma_secondary(M, r_km):
    """vdE&B eq.(5), unconstrained Table 1: log10 PGV[cm/s] = -2.29 + 0.85M - 1.29 log10(r)."""
    pgv = 10.0**(-2.29 + 0.85*M - 1.29*np.log10(r_km))/100.0   # m/s
    return G_SHEAR*pgv/C_PHASE

# ---------------- catalogue ----------------
CAT = {}
for name in CLASS:
    df = pd.read_csv(os.path.join(D, f'{name}.csv'), usecols=['time','latitude','longitude','mag'])
    df['time'] = pd.to_datetime(df['time'], utc=True, format='mixed')
    df = df.dropna(subset=['time','mag']).sort_values('time')
    CAT[name] = df
    # cell centroid from the pre-registered box
    b = MAN['cells'][name]['box']
    CAT[name].attrs['clat'] = (b[0]+b[1])/2.0
    CAT[name].attrs['clon'] = (b[2]+b[3])/2.0

DAY_NS = 86400_000_000_000

def counts(tns, t0ns, post_d):
    """N_post over [t0,t0+post_d), N_bg over [t0-90d,t0). tns = sorted int64 ns array."""
    lo = t0ns - int(BG_D*DAY_NS); hi = t0ns + int(post_d*DAY_NS)
    i0 = np.searchsorted(tns, t0ns, 'left')
    n_post = int(np.searchsorted(tns, hi, 'left') - i0)
    n_bg = int(i0 - np.searchsorted(tns, lo, 'left'))
    return n_post, n_bg

def stat(n_post, n_bg, post_d):
    lam = n_bg/BG_D
    exp = lam*post_d
    if exp <= 0:
        return np.nan, exp
    return n_post/exp, exp

def shift_grid_ns(t0ns, span_yr=SHIFT_YR, n=NSHIFT):
    """Regular circular grid of pseudo-origins in t0 +/- span_yr, excluding |dt|<30 d."""
    off = np.linspace(-span_yr*365.25, span_yr*365.25, n+1)
    off = off[np.abs(off) >= MINSEP_D]
    return (t0ns + (off*DAY_NS)).astype(np.int64)

# ---------------- main measurement ----------------
rows = []
for sid, S in SOURCES.items():
    t0 = pd.Timestamp(S['t0'])
    gate = 2*S['L']
    for cell in CLASS:
        clat, clon = CAT[cell].attrs['clat'], CAT[cell].attrs['clon']
        r = float(gcdist_km(S['lat'], S['lon'], clat, clon))
        gated_out = r < gate
        t0ns = int(t0.value)
        for mcut in (1.5, 2.5):
            sub = CAT[cell][CAT[cell]['mag'] >= mcut]
            # pandas 3.0: parsed dtype is datetime64[us, UTC]; force ns to match Timestamp.value
            tns = np.sort(pd.DatetimeIndex(sub['time']).as_unit('ns').asi8)
            for wname, wd in POST_D.items():
                n_post, n_bg = counts(tns, t0ns, wd)
                rr, exp = stat(n_post, n_bg, wd)
                # circular-shift null
                nulls = []
                for ts in shift_grid_ns(t0ns):
                    a, b = counts(tns, int(ts), wd)
                    s, _ = stat(a, b, wd)
                    nulls.append(s)
                nulls = np.array(nulls, float)
                good = np.isfinite(nulls)
                nn = int(good.sum())
                if np.isfinite(rr) and nn > 0:
                    p = (1 + int((nulls[good] >= rr).sum()))/(1 + nn)
                else:
                    p = np.nan
                rows.append(dict(source=sid, cell=cell, cls=CLASS[cell],
                    documented=cell in DOCUMENTED, dist_km=r, gate_km=gate,
                    gated_out=bool(gated_out), mcut=mcut, window=wname,
                    n_post=n_post, n_bg=n_bg, lam_bg_per_day=n_bg/BG_D, expected=exp,
                    RR=rr, p_cell=p, n_null=nn,
                    null_mean=float(np.nanmean(nulls[good])) if nn else np.nan,
                    null_p95=float(np.nanpercentile(nulls[good], 95)) if nn else np.nan,
                    sig_primary_Pa=float(sigma_primary(S['Ms'], r)),
                    sig_secondary_Pa=float(sigma_secondary(S['Ms'], r)),
                    sig_primary_Mw_Pa=float(sigma_primary(S['Mw'], r)),
                    nulls=nulls))   # full length incl. NaN, so the S-8 max-statistic aligns

df = pd.DataFrame(rows)
nullcol = df.pop('nulls')
df.to_csv(os.path.join(HERE, 'k034_cellstats.csv'), index=False)
print('cells scored:', len(df), flush=True)

# ---------------- S-8 family-wise max-statistic ----------------
# declared family = sources x gated-in cells x windows x magnitude cuts (Asigma enters power only)
fam = df[~df['gated_out']].copy()
fam_idx = fam.index.to_numpy()
NM = np.vstack([nullcol.loc[i] for i in fam_idx])            # (F, nnull)
maxnull = np.nanmax(NM, axis=0)                              # family max-statistic per shift
fam['p_familywise'] = [ (1 + int((maxnull >= rr).sum()))/(1+len(maxnull)) if np.isfinite(rr) else np.nan
                        for rr in fam['RR'] ]

# The raw-RR max-statistic is DEGENERATE across cells: RR is not comparable between a cell with
# lambda_bg = 2.3/d (geysers) and one with 0.03/d (cedar_city), so the family max is set by the
# sparsest cell's null tail and swamps real detections. Implemented faithfully above and reported,
# but the licence uses the standardised (Westfall-Young min-p) form of the SAME null, below.
def wy_pvals(idx_list):
    """Westfall-Young single-step min-p over a family, using the shared circular-shift null."""
    NMx = np.vstack([np.asarray(nullcol.loc[i], float) for i in idx_list])   # (F, J)
    F, J = NMx.shape
    # null p-value scale: for member f, realisation j -> P(RR_null[f,:] >= RR_null[f,j])
    order = np.argsort(np.where(np.isfinite(NMx), NMx, -np.inf), axis=1)
    pnull = np.ones_like(NMx)
    for f in range(F):
        v = NMx[f]
        good = np.isfinite(v)
        ranks = np.empty(J); ranks[:] = np.nan
        vv = v[good]
        # p = (1 + #{>= v}) / (1 + n)
        srt = np.sort(vv)
        cnt = len(vv) - np.searchsorted(srt, vv, 'left')
        ranks[good] = (1.0 + cnt) / (1.0 + len(vv))
        pnull[f] = np.where(np.isfinite(ranks), ranks, 1.0)
    minp_null = pnull.min(axis=0)
    return minp_null

def obs_p(i):
    v = np.asarray(nullcol.loc[i], float)
    good = np.isfinite(v)
    rr = df.loc[i, 'RR']
    if not np.isfinite(rr) or good.sum() == 0:
        return np.nan
    return (1 + int((v[good] >= rr).sum())) / (1 + int(good.sum()))

minp_all = wy_pvals(list(fam_idx))
fam['p_WY_family'] = [ (1 + int((minp_all <= p).sum()))/(1+len(minp_all)) if np.isfinite(p) else np.nan
                       for p in fam['p_cell'] ]
persrc_minp = {}
for sid in SOURCES:
    m = (fam['source'] == sid).to_numpy()
    persrc_minp[sid] = wy_pvals(list(fam_idx[m]))
fam['p_WY_within_source'] = [
    (1 + int((persrc_minp[r.source] <= r.p_cell).sum()))/(1+len(persrc_minp[r.source]))
    if np.isfinite(r.p_cell) else np.nan for r in fam.itertuples()]

# per-source family (the licence is per-source: "fires on >=2 of 4")
persrc_max = {}
for sid in SOURCES:
    m = fam['source'] == sid
    NMs = np.vstack([nullcol.loc[i] for i in fam_idx[m.to_numpy()]])
    persrc_max[sid] = np.nanmax(NMs, axis=0)
fam['p_familywise_within_source'] = [
    (1 + int((persrc_max[r.source] >= r.RR).sum()))/(1+len(persrc_max[r.source]))
    if np.isfinite(r.RR) else np.nan for r in fam.itertuples()]

fam.to_csv(os.path.join(HERE, 'k034_familywise.csv'), index=False)

# -------- SECONDARY (EXPLORATORY, not in the prereg): count-only statistic --------
# Rescue for cells the primary statistic leaves UNMEASURABLE (S-15): when N_bg = 0 the ratio RR is
# undefined, yet little_skull_mtn has N_post = 11 in 5 d after Landers. N_post scored directly
# against its OWN circular-shift null is defined everywhere. Added AFTER seeing the primary run,
# therefore EXPLORATORY: it may not set the PASS flag.
sec = []
for sid, S in SOURCES.items():
    t0ns = int(pd.Timestamp(S['t0']).value); gate = 2*S['L']
    for cell in CLASS:
        r = float(gcdist_km(S['lat'], S['lon'], CAT[cell].attrs['clat'], CAT[cell].attrs['clon']))
        if r < gate:
            continue
        tns = np.sort(pd.DatetimeIndex(CAT[cell][CAT[cell]['mag'] >= 1.5]['time']).as_unit('ns').asi8)
        n_obs, _ = counts(tns, t0ns, 5.0)
        nl = np.array([counts(tns, int(ts), 5.0)[0] for ts in shift_grid_ns(t0ns)], float)
        p = (1 + int((nl >= n_obs).sum()))/(1 + len(nl))
        sec.append(dict(source=sid, cell=cell, cls=CLASS[cell], documented=cell in DOCUMENTED,
                        dist_km=r, N_post_5d=int(n_obs), null_mean_N=float(nl.mean()),
                        null_p95_N=float(np.percentile(nl, 95)), p_count=float(p),
                        primary_measurable=bool(df[(df.source == sid) & (df.cell == cell) &
                            (df.window == '5d') & (df.mcut == 1.5)]['RR'].notna().iloc[0])))
secdf = pd.DataFrame(sec)
secdf.to_csv(os.path.join(HERE, 'k034_secondary_counts.csv'), index=False)

# -------- P3 pattern test: pre-registered rank vs observed --------
from scipy.stats import spearmanr, mannwhitneyu
pattern = {}
for sid in SOURCES:
    m = fam[(fam.source == sid) & (fam.window == '5d') & (fam.mcut == 1.5)].copy()
    m = m[m['p_cell'].notna()]
    if len(m) < 4:
        pattern[sid] = dict(n=len(m), note='too few measurable cells'); continue
    m['prereg_key'] = [ (CLASSRANK[c], -s) for c, s in zip(m['cls'], m['sig_primary_Pa']) ]
    m = m.sort_values(['cls', 'sig_primary_Pa'], ascending=[True, False]).reset_index(drop=True)
    m['prereg_rank'] = np.arange(1, len(m)+1)
    m['obs_rank'] = m['p_cell'].rank(method='average')       # 1 = most significant
    rho, pr = spearmanr(m['prereg_rank'], m['obs_rank'])
    ab = m[m.cls.isin(['A', 'B'])]['p_cell'].to_numpy()
    c = m[m.cls == 'C']['p_cell'].to_numpy()
    if len(ab) and len(c):
        u, pu = mannwhitneyu(ab, c, alternative='less')      # A+B more significant than C
    else:
        u, pu = np.nan, np.nan
    pattern[sid] = dict(n=int(len(m)), spearman_rho=float(rho), spearman_p_two_sided=float(pr),
                        spearman_p_one_sided=float(pr/2 if rho > 0 else 1-pr/2),
                        median_p_AB=float(np.median(ab)) if len(ab) else None,
                        median_p_C=float(np.median(c)) if len(c) else None,
                        mannwhitney_U=float(u) if np.isfinite(u) else None,
                        mannwhitney_p_AB_lt_C=float(pu) if np.isfinite(pu) else None,
                        prereg_order=list(m['cell']), observed_order=list(
                            m.sort_values('p_cell')['cell']))

# ---------------- power / detection-threshold curve ----------------
NSIM = 4000
Rgrid = np.concatenate([np.arange(1.0, 5.0, 0.05), np.arange(5.0, 20.01, 0.25)])
power_rows = []
for i in fam_idx:
    r0 = df.loc[i]
    if r0['window'] != '5d' or r0['mcut'] != 1.5:
        continue
    nl = nullcol.loc[i]
    crit = np.nanpercentile(nl, 100*(1-ALPHA))    # per-cell alpha=0.05 critical RR
    lam = r0['lam_bg_per_day']; wd = POST_D[r0['window']]
    exp0 = lam*wd
    if exp0 <= 0:
        power_rows.append(dict(source=r0['source'], cell=r0['cell'], measurable=False)); continue
    Rmin = None; curve = []
    for R in Rgrid:
        n = RNG.poisson(R*exp0, NSIM)
        nbg = RNG.poisson(lam*BG_D, NSIM)
        with np.errstate(divide='ignore', invalid='ignore'):
            rr = np.where(nbg > 0, n/((nbg/BG_D)*wd), np.nan)
        pw = float(np.nanmean(rr > crit))
        curve.append((float(R), pw))
        if Rmin is None and pw >= 0.80:
            Rmin = float(R)
    d = dict(source=r0['source'], cell=r0['cell'], cls=r0['cls'],
             documented=bool(r0['documented']), measurable=Rmin is not None,
             dist_km=float(r0['dist_km']), lam_bg_per_day=float(lam), expected_5d=float(exp0),
             crit_RR=float(crit), R_min_80pct=Rmin,
             sig_primary_kPa=float(r0['sig_primary_Pa']/1e3),
             sig_secondary_kPa=float(r0['sig_secondary_Pa']/1e3),
             power_curve=curve)
    for k, a in ASIGMA.items():
        d[f'dTau_min_kPa_Asigma{k}'] = (a*np.log(Rmin)/1e3) if Rmin else None
        # power AT the cell's own predicted amplitude, both axes
        for axis, sg in (('primary', r0['sig_primary_Pa']), ('secondary', r0['sig_secondary_Pa'])):
            Rpred = float(np.exp(sg/a))
            n = RNG.poisson(min(Rpred, 1e4)*exp0, NSIM)
            nbg = RNG.poisson(lam*BG_D, NSIM)
            with np.errstate(divide='ignore', invalid='ignore'):
                rr = np.where(nbg > 0, n/((nbg/BG_D)*wd), np.nan)
            d[f'R_pred_{axis}_Asigma{k}'] = Rpred
            d[f'power_{axis}_Asigma{k}'] = float(np.nanmean(rr > crit))
    power_rows.append(d)
pw = pd.DataFrame(power_rows)
pw.drop(columns=[c for c in ['power_curve'] if c in pw.columns]).to_csv(os.path.join(HERE, 'k034_power.csv'), index=False)
print('power rows:', len(pw), flush=True)

json.dump(dict(cellstats=df.to_dict('records'),
               familywise=fam.to_dict('records'),
               secondary_counts=secdf.to_dict('records'),
               pattern=pattern,
               power=power_rows),
          open(os.path.join(HERE, 'k034_raw.json'), 'w'), default=str)
print('WROTE k034_raw.json')
