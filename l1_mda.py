"""L-1 step 2: Branch-E link in natural-log form, MDA per cell, bracket."""
import numpy as np, pandas as pd
from scipy import optimize, stats

cs = pd.read_csv('k034_cellstats.csv'); pw = pd.read_csv('k034_power.csv')
prim = cs[(cs.window == '5d') & (cs.mcut == 1.5) & (~cs.gated_out)].copy()
prim['sig_kPa'] = prim.sig_primary_Pa / 1e3
prim['grp'] = np.where(prim.cls.isin(['A', 'B']), 'AB', 'C')
prim['meas'] = prim.n_bg > 0
pw['grp'] = np.where(pw.cls.isin(['A', 'B']), 'AB', 'C')
print('power measurable by grp:\n', pw.groupby(['grp', 'measurable']).size())

def poisfit(df):
    y = df.n_post.values.astype(float); off = np.log(df.expected.values)
    x = np.log(df.sig_kPa.values)            # NATURAL log -> R = e^a * sig^b
    X = np.column_stack([np.ones_like(x), x])
    nll = lambda p: np.sum(np.exp(np.clip(off + X @ p, -50, 50)) - y * (off + X @ p))
    p = optimize.minimize(nll, [0., 0.], method='BFGS').x
    mu = np.exp(off + X @ p)
    Vn = np.linalg.inv(X.T @ (X * mu[:, None]))
    meat = np.zeros((2, 2))
    for s, idx in df.groupby('source').indices.items():
        u = (X[idx] * (y[idx] - mu[idx])[:, None]).sum(0); meat += np.outer(u, u)
    Vc = Vn @ meat @ Vn
    dev = 2*np.sum(np.where(y > 0, y*np.log(np.where(y > 0, y/mu, 1)), 0) - (y-mu))
    nll0 = lambda p0: np.sum(np.exp(off + p0[0]) - y*(off + p0[0]))
    a0 = optimize.minimize(nll0, [0.], method='BFGS').x[0]
    mu0 = np.exp(off + a0)
    dev0 = 2*np.sum(np.where(y > 0, y*np.log(np.where(y > 0, y/mu0, 1)), 0) - (y-mu0))
    lrt = 2*(nll0([a0]) - nll(p))
    disp = np.sum((y-mu)**2/mu)/(len(y)-2)
    return dict(a=float(p[0]), b=float(p[1]),
                se_naive=[float(v) for v in np.sqrt(np.diag(Vn))],
                se_cluster_by_source=[float(v) for v in np.sqrt(np.diag(Vc))],
                se_quasi=[float(v) for v in np.sqrt(np.diag(Vn)*disp)],
                deviance=float(dev), null_deviance=float(dev0),
                pseudo_R2_dev=float(1-dev/dev0), dispersion_pearson=float(disp),
                LRT_b=float(lrt), p_LRT_b_naive=float(stats.chi2.sf(lrt, 1)),
                n_pairs=int(len(y)), a_null_only=float(a0))

def cboot(df, nb=4000, seed=20260811):
    rng = np.random.default_rng(seed); srcs = df.source.unique(); B = []
    for _ in range(nb):
        pick = rng.choice(srcs, len(srcs), replace=True)
        dd = pd.concat([df[df.source == s].assign(source=f'{s}_{i}') for i, s in enumerate(pick)])
        if dd.n_post.sum() == 0: continue
        try:
            f = poisfit(dd); B.append([f['a'], f['b']])
        except Exception: pass
    return np.array(B)

res = {}
for g in ['AB', 'C']:
    d = prim[(prim.grp == g) & prim.meas]
    f = poisfit(d); B = cboot(d)
    f['boot_a_CI95'] = [float(v) for v in np.percentile(B[:, 0], [2.5, 97.5])]
    f['boot_b_CI95'] = [float(v) for v in np.percentile(B[:, 1], [2.5, 97.5])]
    f['boot_n'] = int(len(B))
    f['sigma_support_kPa'] = [float(d.sig_kPa.min()), float(d.sig_kPa.max())]
    f['median_RR'] = float(d.RR.median())
    m = d[(d.RR > 0)]
    sl, ic, r, p, se = stats.linregress(np.log10(m.sig_kPa), np.log10(m.RR))
    f['ols_logRR_slope_per_decade'] = float(sl); f['ols_logRR_p'] = float(p)
    f['ols_logRR_r2'] = float(r**2); f['ols_n'] = int(len(m))
    res[g] = f
    print(f'\n--- {g} ---')
    for k, v in f.items(): print(' ', k, v)

# ---- MDA per cell ----
a, b = res['AB']['a'], res['AB']['b']
pw['MDA_E_kPa'] = np.where(pw.measurable, (pw.R_min_80pct/np.exp(a))**(1.0/b), np.nan)
for asig in ['0.03', '0.10', '0.15']:
    pw[f'MDA_S_{asig}'] = np.where(pw.measurable, pw[f'dTau_min_kPa_Asigma{asig}'], np.nan)

def summ(df, col):
    v = df[col].dropna()
    return dict(best=float(v.min()), median=float(v.median()), worst=float(v.max()), n=int(len(v)))

ab = pw[(pw.grp == 'AB') & pw.measurable]; allm = pw[pw.measurable]
print('\nA+B measurable cells:', len(ab), ' all measurable:', len(allm))
for name, df in [('AB', ab), ('ALL', allm)]:
    print(name, 'R_min_80', summ(df, 'R_min_80pct'))
    for c in ['MDA_S_0.03', 'MDA_S_0.10', 'MDA_S_0.15', 'MDA_E_kPa']:
        print('  ', c, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in summ(df, c).items()})
pw.to_csv(r'C:\Users\jimgale\AppData\Local\Temp\claude\D--CODE-git-quake\d13596dc-3656-4a0e-86a3-168a941d435b\scratchpad\l1_mda_percell.csv', index=False)
import json
json.dump(res, open(r'C:\Users\jimgale\AppData\Local\Temp\claude\D--CODE-git-quake\d13596dc-3656-4a0e-86a3-168a941d435b\scratchpad\l1_fits.json', 'w'), indent=1)
