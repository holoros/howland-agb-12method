#!/usr/bin/env python3
"""Close the open analysis items of the Howland AGB projection manuscript from public data on Cardinal.

Inputs: ORNL DAAC 2490 version 2 maps, DigiSylva 2.0 plots (plots.tab), and the AcadianGY 12.3.6
climate-site-index trajectories from the October 23 check (csi_trajectories.csv). The study area
polygon is not public, so the domain is the convex hull of the Howland plots buffered by 30 m (the
buffer whose cell count is closest to the paper's 1,967), with 0 and 60 m as sensitivity.

Produces, for the map-initialized methods that can be rebuilt from the paper's own curve parameters
(static, constant increment, area-based curve refit on the domain and with the paper's parameters,
FVS-ACD and FVS-NE curves with the Table 3 parameters, and trajectory imputation from the CSI 12 m
AcadianGY trajectories): scores at 2021 and 2025 with block bootstrap intervals and block counts,
paired contrasts at 90, 270 and 540 m, offset-adjusted rescoring with the offset re-estimated in each
replicate and at fixed 3.5, 4.3 and 5.1, change efficiency, plot-level scoring at the 2026 plots,
domain and registration sensitivity of the offset, three-cornered hat for every triplet, and the
FVS-ACD curve increment at CSI 10, 12 and 14 m. Units Mg C ha-1.
"""
import json, os, sys, hashlib
import numpy as np, pandas as pd, rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import Affine
from scipy.optimize import least_squares
from shapely.geometry import MultiPoint
from shapely import vectorized

D = sys.argv[1]; CSI = sys.argv[2]; OUT = sys.argv[3]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261026)
YEARS = [2012, 2015, 2017, 2021, 2023, 2025]
PAPER = {"area_based_90m": (111.4, 0.0405, 10.0), "fvs_ne": (180.0, 0.0224, 10.0),
         "fvs_acd": (365.3, 0.0023, 0.86), "fvs_acd_ingrowth": (1000.0, 0.0005, 0.89)}
CARBON = 0.5  # [ASSUMPTION] dry to carbon for the AcadianGY trajectories; paper used species fractions 0.48 to 0.52
md5 = {f: hashlib.md5(open(os.path.join(D, f), 'rb').read()).hexdigest() for f in sorted(os.listdir(D)) if f.endswith(('.tif', '.tab'))}
md5['csi_trajectories.csv'] = hashlib.md5(open(CSI, 'rb').read()).hexdigest()
json.dump(md5, open(f"{OUT}/inputs_md5.json", 'w'), indent=1)

def cr(a1, A, k, pp, dt):
    x = np.clip(a1 / A, 1e-9, None) ** (1 / pp)
    return A * (1 - (1 - x) * np.exp(-k * dt)) ** pp

def load_maps(shift):
    with rasterio.open(f"{D}/Howland_AGB_2017.tif") as r:
        ref_t, ref_shape, ref_crs = r.transform, r.shape, r.crs
    maps = {}
    for y in YEARS:
        with rasterio.open(f"{D}/Howland_AGB_{y}.tif") as r:
            a = r.read(1).astype('float64'); a[a == r.nodata] = np.nan
            t = r.transform; crs = r.crs or ref_crs
            if y == 2025 and shift:
                t = Affine(t.a, t.b, t.c + 5.0, t.d, t.e, t.f - 5.0)
            dst = np.full(ref_shape, np.nan)
            reproject(a, dst, src_transform=t, src_crs=crs, dst_transform=ref_t, dst_crs=ref_crs,
                      resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
            maps[y] = dst * 10.0
    return maps, ref_t, ref_shape

def to_cells(a, h3, w3):
    b = a[: h3 * 3, : w3 * 3].reshape(h3, 3, w3, 3)
    complete = np.isfinite(b).all(axis=(1, 3))
    return np.where(complete, np.nanmean(b, axis=(1, 3)), np.nan)

# plots
p = pd.read_csv(f"{D}/plots.tab", sep=None, engine='python')
t = pd.read_csv(f"{D}/trees.tab", sep=None, engine='python')
hl = p[p.Plot.str.startswith('HL')].copy()
t['plot_id'] = t.Plot.str.rsplit('_', n=1).str[0]
hl['meas_year'] = hl.Plot.map(t[t.plot_id.str.startswith('HL')].groupby('plot_id').Year.max())
for y in (2017, 2021, 2023, 2025):
    hl[f'AGB_{y}_J'] = hl[f'AGB_{y}_J'] * 10.0
hull = MultiPoint(list(zip(hl.Easting, hl.Northing))).convex_hull

# trajectories (CSI 12) for trajectory imputation and CSI sensitivity
tr = pd.read_csv(CSI)
tr['agbC'] = tr.agb_J * CARBON
res = {'carbon_fraction_assumed': CARBON, 'csi_plots': int(tr.PLOT.nunique())}

def curve_from_traj(csi, dt):
    s = tr[tr.CSI == csi].sort_values(['PLOT', 'year'])
    pairs = []
    for pl, g in s.groupby('PLOT'):
        v = g.agbC.to_numpy(); yrs = g.year.to_numpy()
        for i in range(len(v) - dt):
            pairs.append((v[i], v[i + dt]))
    pr = np.array(pairs)
    fit = least_squares(lambda th: cr(pr[:, 0], *th, dt) - pr[:, 1], x0=[300, 0.003, 0.9], bounds=([75, 1e-4, 0.3], [1000, 1, 10]))
    A, k, pp = fit.x
    return {'A': A, 'k': k, 'p': pp, 'pairs': len(pr), 'increment_at_82.3': float((cr(82.3, A, k, pp, 8) - 82.3) / 8)}
res['fvs_acd_curve_by_csi'] = {str(c): curve_from_traj(c, 8) for c in sorted(tr.CSI.unique())}

def traj_impute(s17, dt, csi=12, nn=4):
    """IDW of the 4 simulated plot states closest in AGB; increments over dt years from the CSI trajectories."""
    s = tr[tr.CSI == csi].sort_values(['PLOT', 'year'])
    states, incs = [], []
    for pl, g in s.groupby('PLOT'):
        v = g.agbC.to_numpy()
        for i in range(len(v) - dt):
            states.append(v[i]); incs.append(v[i + dt] - v[i])
    states, incs = np.array(states), np.array(incs)
    out = np.empty_like(s17)
    for j, a in enumerate(s17):
        d = np.abs(states - a); idx = np.argpartition(d, nn)[:nn]
        w = 1 / np.maximum(d[idx], 1e-3)
        out[j] = a + (w * incs[idx]).sum() / w.sum()
    return out

def score(pred, obs):
    e = pred - obs
    return {'bias': float(e.mean()), 'rmse_pct': float(100 * np.sqrt((e ** 2).mean()) / obs.mean()),
            'r2': float(1 - (e ** 2).sum() / ((obs - obs.mean()) ** 2).sum())}

def run(buf, shift, main=False):
    maps, ref_t, ref_shape = load_maps(shift)
    H, W = ref_shape; h3, w3 = H // 3, W // 3
    C = {y: to_cells(maps[y], h3, w3) for y in YEARS}
    rows, cols = np.meshgrid(np.arange(h3), np.arange(w3), indexing='ij')
    xs = ref_t.c + (cols * 3 + 1.5) * ref_t.a; ys = ref_t.f + (rows * 3 + 1.5) * ref_t.e
    g = hull.buffer(buf)
    dom = vectorized.contains(g, xs, ys) & np.isfinite(C[2017]) & np.isfinite(C[2025]) & np.isfinite(C[2012])
    all6 = dom & np.all([np.isfinite(C[y]) for y in YEARS], axis=0)
    c21 = dom & np.isfinite(C[2021])
    V = {y: C[y][all6] for y in YEARS}
    r = {'buffer_m': buf, 'shift_2025': shift, 'area_ha': g.area / 1e4, 'cells_dom': int(dom.sum()), 'cells_all6': int(all6.sum()),
         'means_all6': {str(y): float(V[y].mean()) for y in YEARS}}
    yrs = np.array([2012, 2017, 2021.0]); Y = np.vstack([V[2012], V[2017], V[2021]]).T
    xb = yrs.mean(); b = ((yrs - xb) * (Y - Y.mean(1, keepdims=True))).sum(1) / ((yrs - xb) ** 2).sum(); a0 = Y.mean(1) - b * xb
    r['trend_rate'] = float(b.mean())
    r['offset_2025_trend'] = float(np.mean(V[2025] - (a0 + b * 2025)))
    r['offset_2025_rate_2012_2017'] = float(np.mean(V[2025] - (V[2017] + (V[2017] - V[2012]) / 5 * 8)))
    r['offset_2025_rate_2017_2021'] = float(np.mean(V[2025] - (V[2017] + (V[2021] - V[2017]) / 4 * 8)))
    r['offsets_3dep'] = {str(y): float(np.mean(V[y] - (a0 + b * y))) for y in (2015, 2023)}
    r['mean_change_2017_2025'] = float(np.mean(V[2025] - V[2017]))
    # offset interval, 270 m blocks
    bid = (np.floor(xs[all6] / 270).astype(int) * 100000 + np.floor(ys[all6] / 270).astype(int))
    ub = np.unique(bid); groups = [np.flatnonzero(bid == u) for u in ub]
    offs = []
    dev = V[2025] - (a0 + b * 2025)
    for _ in range(1000):
        idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))]); offs.append(dev[idx].mean())
    r['offset_2025_ci95_270m'] = [float(np.percentile(offs, 2.5)), float(np.percentile(offs, 97.5))]
    for y0, y1 in [(2017, 2021), (2017, 2025), (2023, 2025)]:
        m = dom & np.isfinite(C[y0]) & np.isfinite(C[y1])
        r[f'loss_gt20pct_{y0}_{y1}'] = float(np.mean((C[y1][m] - C[y0][m]) / C[y0][m] < -0.2))
    if not main:
        return r, None
    # ---------- methods ----------
    pix = np.zeros(ref_shape, bool); pix[: h3 * 3, : w3 * 3] = np.repeat(np.repeat(dom, 3, 0), 3, 1)
    inc = float(np.nanmean((maps[2017] - maps[2012])[pix]) / 5); r['constant_increment'] = inc
    h9, w9 = h3 // 3, w3 // 3
    def blocks(a, mask):
        aa = np.where(mask, a, np.nan)[: h9 * 3, : w9 * 3].reshape(h9, 3, w9, 3)
        full = np.isfinite(aa).all(axis=(1, 3)); return np.where(full, np.nanmean(aa, axis=(1, 3)), np.nan)
    B12, B17 = blocks(C[2012], dom), blocks(C[2017], dom); ok = np.isfinite(B12) & np.isfinite(B17)
    fit = least_squares(lambda th: cr(B12[ok], *th, 5.0) - B17[ok], x0=[110, 0.04, 3], bounds=([75, 1e-4, 0.3], [1000, 1, 10]))
    Af, kf, pf = fit.x; r['area_based_refit'] = {'A': Af, 'k': kf, 'p': pf, 'pairs': int(ok.sum())}
    def methods(s17, dt):
        return {'static': s17, 'constinc': s17 + inc * dt, '2_refit': cr(s17, Af, kf, pf, dt),
                '2_paper': cr(s17, *PAPER['area_based_90m'], dt), '3_fvs_acd': cr(s17, *PAPER['fvs_acd'], dt),
                '3b_fvs_ne': cr(s17, *PAPER['fvs_ne'], dt), '3_ingrowth': cr(s17, *PAPER['fvs_acd_ingrowth'], dt),
                '9_traj_imp': traj_impute(s17, dt)}
    tables = []; preds = {}
    for yv, dt, m in [(2021, 4, c21), (2025, 8, dom)]:
        s17, obs = C[2017][m], C[yv][m]; P = methods(s17, dt); preds[yv] = (P, obs, m)
        varch = np.var(obs - s17)
        for kk, v in P.items():
            sc = score(v, obs); sc['E'] = float(1 - ((v - obs) ** 2).mean() / varch)
            tables.append({'year': yv, 'offset_removed': 0.0, 'method': kk, 'n_cells': int(m.sum()), **sc})
    # offset-adjusted at 2025: domain estimate and fixed values; bootstrap with offset re-estimated per replicate (270 m blocks for offset, 90 m for scores)
    P, obs, m = preds[2025]; s_est = r['offset_2025_trend']
    for s in [s_est, 3.5, 4.3, 5.1]:
        varch = np.var(obs - s - P['static'])
        for kk, v in P.items():
            sc = score(v, obs - s); sc['E'] = float(1 - ((v - (obs - s)) ** 2).mean() / varch)
            tables.append({'year': 2025, 'offset_removed': float(s), 'method': kk, 'n_cells': int(m.sum()), **sc})
    pd.DataFrame(tables).to_csv(f"{OUT}/method_scores.csv", index=False)
    # block bootstrap: RMSE% intervals, paired contrasts, offset-adjusted with re-estimated offset
    bb = []; sc_ci = []
    for size in (90, 270, 540):
        for yv in (2021, 2025):
            P, obs, m = preds[yv]
            bid = (np.floor(xs[m] / size).astype(int) * 100000 + np.floor(ys[m] / size).astype(int))
            ub = np.unique(bid); groups = [np.flatnonzero(bid == u) for u in ub]
            # map all6 deviation onto m for per-replicate offset: use cells in both
            keys = list(P.keys()); rm_all = {k: [] for k in keys}; rm_adj = {k: [] for k in keys}; bias_all = {k: [] for k in keys}
            diffs = {(a, b_): [] for a in keys for b_ in ('static', 'constinc') if a != b_}
            if yv == 2025:
                # offset per replicate from the 270 m trend deviation restricted to resampled cells needs all6; approximate with cell-level deviation where available
                dev_full = np.full(m.shape, np.nan)
                tmp = np.full(all6.sum(), np.nan); tmp[:] = V[2025] - (a0 + b * 2025); dev_full[all6] = tmp
                dev_m = dev_full[m]
            for _ in range(1000):
                idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
                o = obs[idx]; om = o.mean()
                rm = {k: 100 * np.sqrt(((v[idx] - o) ** 2).mean()) / om for k, v in P.items()}
                for k in keys:
                    rm_all[k].append(rm[k]); bias_all[k].append((P[k][idx] - o).mean())
                for key in diffs: diffs[key].append(rm[key[0]] - rm[key[1]])
                if yv == 2025:
                    s_rep = np.nanmean(dev_m[idx]); oa = o - s_rep
                    for k in keys: rm_adj[k].append(100 * np.sqrt(((P[k][idx] - oa) ** 2).mean()) / oa.mean())
            for k in keys:
                sc_ci.append({'block_m': size, 'n_blocks': len(ub), 'year': yv, 'method': k,
                              'rmse_pct_lo95': np.percentile(rm_all[k], 2.5), 'rmse_pct_hi95': np.percentile(rm_all[k], 97.5),
                              'bias_lo95': np.percentile(bias_all[k], 2.5), 'bias_hi95': np.percentile(bias_all[k], 97.5),
                              'rmse_adj_lo95': np.percentile(rm_adj[k], 2.5) if yv == 2025 else np.nan,
                              'rmse_adj_hi95': np.percentile(rm_adj[k], 97.5) if yv == 2025 else np.nan,
                              'rmse_adj_mean': np.mean(rm_adj[k]) if yv == 2025 else np.nan})
            for key, v in diffs.items():
                v = np.array(v)
                bb.append({'block_m': size, 'n_blocks': len(ub), 'year': yv, 'contrast': f'{key[0]}-{key[1]}',
                           'mean': v.mean(), 'lo95': np.percentile(v, 2.5), 'hi95': np.percentile(v, 97.5)})
    pd.DataFrame(bb).to_csv(f"{OUT}/contrasts.csv", index=False)
    pd.DataFrame(sc_ci).to_csv(f"{OUT}/score_intervals.csv", index=False)
    # three-cornered hat, all triplets
    import itertools
    tch = []
    for tri in itertools.combinations(YEARS, 3):
        a_, b_, c_ = (V[y] - V[y].mean() for y in tri)
        vab, vac, vbc = np.var(a_ - b_), np.var(a_ - c_), np.var(b_ - c_)
        est = {tri[0]: (vab + vac - vbc) / 2, tri[1]: (vab + vbc - vac) / 2, tri[2]: (vac + vbc - vab) / 2}
        tch.append({'triplet': '_'.join(map(str, tri)), **{f'sd_{y}': (float(np.sqrt(v)) if v > 0 else np.nan) for y, v in est.items()}})
    pd.DataFrame(tch).to_csv(f"{OUT}/three_cornered_hat_all.csv", index=False)
    # plots: every method at the plot cell, 2021 (back-cast) and 2025
    inv = ~ref_t; pr = []
    for _, r0 in hl.iterrows():
        col, row = inv * (r0.Easting, r0.Northing); col, row = int(col), int(row)
        rec = {'plot': r0.Plot, 'meas_year': r0.meas_year}
        for y in YEARS:
            w = maps[y][max(row - 1, 0): row + 2, max(col - 1, 0): col + 2]
            rec[f'map_{y}'] = float(np.nanmean(w)) if np.isfinite(w).any() else np.nan
        for y in (2017, 2021, 2023, 2025): rec[f'plot_{y}'] = r0[f'AGB_{y}_J']
        pr.append(rec)
    pr = pd.DataFrame(pr); okp = pr.map_2017.notna()
    for yv, dt in [(2021, 4), (2025, 8)]:
        P = methods(pr.map_2017.to_numpy(), dt)
        for k, v in P.items(): pr[f'pred_{k}_{yv}'] = v
    pr.to_csv(f"{OUT}/plot_predictions.csv", index=False)
    ps = []
    for grp, sel in [('first_measured_2026', pr.meas_year == 2026), ('measured_2024_2025', pr.meas_year.isin([2024, 2025])), ('all', pr.meas_year.notna())]:
        for yv in (2021, 2025):
            s = pr[sel & okp & pr[f'plot_{yv}'].notna() & pr[f'map_{yv}'].notna()]
            if len(s) < 3: continue
            obs = s[f'plot_{yv}'].to_numpy()
            for nm in [f'map_{yv}'] + [f'pred_{k}_{yv}' for k in ['static', 'constinc', '2_refit', '2_paper', '3_fvs_acd', '3b_fvs_ne', '9_traj_imp']]:
                pred = s[nm].to_numpy(); e = pred - obs
                bs = [e[rng.integers(0, len(e), len(e))].mean() for _ in range(2000)]
                ps.append({'plots': grp, 'year': yv, 'n': len(s), 'prediction': nm, 'plot_mean': obs.mean(), 'bias': e.mean(),
                           'bias_lo95': np.percentile(bs, 2.5), 'bias_hi95': np.percentile(bs, 97.5),
                           'rmse_pct': 100 * np.sqrt((e ** 2).mean()) / obs.mean(), 'r': float(np.corrcoef(pred, obs)[0, 1])})
    pd.DataFrame(ps).to_csv(f"{OUT}/plot_scores.csv", index=False)
    # plot vs map increments at the same plots
    s = pr[okp & pr.map_2021.notna() & pr.plot_2017.notna() & pr.plot_2021.notna()]
    r['plot_vs_map_increment_2017_2021'] = {'n': int(len(s)), 'plot': float((s.plot_2021 - s.plot_2017).mean()), 'map': float((s.map_2021 - s.map_2017).mean())}
    s = pr[okp & pr.map_2025.notna() & pr.plot_2017.notna() & pr.plot_2025.notna()]
    r['plot_vs_map_increment_2017_2025'] = {'n': int(len(s)), 'plot': float((s.plot_2025 - s.plot_2017).mean()), 'map': float((s.map_2025 - s.map_2017).mean())}
    r['plots_by_meas_year'] = {str(k): int(v) for k, v in hl.meas_year.value_counts().sort_index().items()}
    return r, None

summary = {'main': run(30, True, main=True)[0], 'sensitivity': []}
for buf in (0, 30, 60):
    for shift in (True, False):
        if buf == 30 and shift: continue
        summary['sensitivity'].append(run(buf, shift)[0])
summary['inputs'] = res
json.dump(summary, open(f"{OUT}/summary.json", 'w'), indent=1, default=float)
print(json.dumps(summary, indent=1, default=float))
for f in ['method_scores.csv', 'contrasts.csv', 'score_intervals.csv', 'plot_scores.csv', 'three_cornered_hat_all.csv']:
    print('\n==', f); print(pd.read_csv(f"{OUT}/{f}").round(2).to_string())
