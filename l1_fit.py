"""L-1: S-14(c) two-branch link bracket. Branch-E fit exploration."""
import numpy as np, pandas as pd
from scipy import optimize, stats

cs = pd.read_csv('k034_cellstats.csv')
pw = pd.read_csv('k034_power.csv')
print('cellstats', cs.shape, 'power', pw.shape)
print('gated_out counts:', cs.gated_out.value_counts().to_dict())

prim = cs[(cs.window == '5d') & (cs.mcut == 1.5) & (~cs.gated_out)].copy()
print('primary rows (5d, M>=1.5, not gated):', len(prim))
prim['sig_kPa'] = prim.sig_primary_Pa / 1e3
prim['grp'] = np.where(prim.cls.isin(['A', 'B']), 'AB', 'C')
prim['measurable'] = prim.n_bg > 0
print(prim.groupby(['grp', 'measurable']).size())

def poisfit(df):
    """n_post ~ Pois(expected * exp(a + b*log10 sig))"""
    y = df.n_post.values.astype(float)
    off = np.log(df.expected.values)
    x = np.log10(df.sig_kPa.values)
    X = np.column_stack([np.ones_like(x), x])
    def nll(p):
        eta = off + X @ p
        mu = np.exp(np.clip(eta, -50, 50))
        return np.sum(mu - y * eta)
    r = optimize.minimize(nll, [0.0, 0.0], method='BFGS')
    p = r.x
    eta = off + X @ p
    mu = np.exp(eta)
    # naive information
    I = X.T @ (X * mu[:, None])
    Vn = np.linalg.inv(I)
    # cluster-robust by source
    meat = np.zeros((2, 2))
    for s, idx in df.groupby('source').indices.items():
        u = (X[idx] * (y[idx] - mu[idx])[:, None]).sum(axis=0)
        meat += np.outer(u, u)
    Vc = Vn @ meat @ Vn
    dev = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y / mu, 1)), 0) - (y - mu))
    # null model (b=0)
    def nll0(p0):
        eta = off + p0[0]
        mu0 = np.exp(np.clip(eta, -50, 50))
        return np.sum(mu0 - y * eta)
    r0 = optimize.minimize(nll0, [0.0], method='BFGS')
    a0 = r0.x[0]
    mu0 = np.exp(off + a0)
    dev0 = 2 * np.sum(np.where(y > 0, y * np.log(np.where(y > 0, y / mu0, 1)), 0) - (y - mu0))
    lrt = 2 * (nll0([a0]) - nll(p))
    pear = np.sum((y - mu) ** 2 / mu)
    return dict(a=p[0], b=p[1], se_naive=np.sqrt(np.diag(Vn)).tolist(),
                se_cluster=np.sqrt(np.diag(Vc)).tolist(),
                deviance=dev, null_deviance=dev0, df_resid=len(y) - 2,
                dispersion_pearson=pear / (len(y) - 2),
                lrt_b=lrt, p_lrt_b=float(stats.chi2.sf(lrt, 1)),
                n=len(y), a0=a0)

def olsfit(df):
    m = df[(df.n_bg > 0) & (df.RR > 0)]
    x = np.log10(m.sig_kPa.values); y = np.log10(m.RR.values)
    b, a, r, p, se = stats.linregress(x, y)
    return dict(a=a, b=b, r2=r**2, p=p, se_b=se, n=len(m))

out = {}
for g in ['AB', 'C']:
    d = prim[(prim.grp == g) & prim.measurable]
    print('\n=== group', g, 'n=', len(d))
    f = poisfit(d); print('poisson', {k: (np.round(v, 4) if isinstance(v, float) else v) for k, v in f.items()})
    o = olsfit(d); print('ols-logRR', {k: np.round(v, 4) for k, v in o.items()})
    out[g] = (f, o)
    # raw description
    print('median RR', d.RR.median(), 'sig range kPa', d.sig_kPa.min().round(1), d.sig_kPa.max().round(1))

# cluster bootstrap over sources
def boot(d, nb=4000, seed=20260811):
    rng = np.random.default_rng(seed)
    srcs = d.source.unique()
    B = []
    for _ in range(nb):
        pick = rng.choice(srcs, len(srcs), replace=True)
        dd = pd.concat([d[d.source == s].assign(source=s + f'_{i}') for i, s in enumerate(pick)])
        try:
            f = poisfit(dd); B.append([f['a'], f['b']])
        except Exception:
            pass
    B = np.array(B)
    return B

for g in ['AB', 'C']:
    d = prim[(prim.grp == g) & prim.measurable]
    B = boot(d, nb=1000)
    print(g, 'boot a CI', np.percentile(B[:, 0], [2.5, 50, 97.5]).round(3),
          'b CI', np.percentile(B[:, 1], [2.5, 50, 97.5]).round(3), 'nb', len(B))
