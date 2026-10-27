#!/usr/bin/env python3
"""Independent reproduction of the map-only numbers in the October 22 Howland manuscript,
from the public ORNL DAAC 2490 version 2 maps and DigiSylva version 2.0 plots.

Covers what does not need the authors' method predictions: map means, G-LiHT rates and trend,
project offsets, static and constant-increment references, the area-based curve (Method 2),
deciles of observed change binned on 2017, 2015 and 2012, loss screen, map pair agreement,
three-cornered hat, block counts and paired contrasts, and plot-level scoring of the
map-initialized predictions at plots first measured in 2026 (independent of the 2024 and 2025
calibration plots named in the DAAC guide).

The 227.6 ha study area polygon is not public, so the domain is approximated by the convex
hull of the Howland plots, with buffers chosen to bracket the paper's 1,967 cells. Every
output states the domain used. Values in Mg C ha-1 (map kg C m-2 x 10; plot kg C m-2 x 10).
"""
import json
import sys

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import Affine
from scipy.optimize import least_squares
from shapely.geometry import MultiPoint, Point
from shapely import vectorized

D = sys.argv[1] if len(sys.argv) > 1 else "."
OUT = sys.argv[2] if len(sys.argv) > 2 else "repro_out"
import os
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20261023)
YEARS = [2012, 2015, 2017, 2021, 2023, 2025]
res = {}

# ---------------- maps onto the 2017 grid ----------------
with rasterio.open(f"{D}/Howland_AGB_2017.tif") as r:
    ref_t, ref_shape, ref_crs = r.transform, r.shape, r.crs
maps = {}
for y in YEARS:
    with rasterio.open(f"{D}/Howland_AGB_{y}.tif") as r:
        a = r.read(1).astype("float64")
        a[a == r.nodata] = np.nan
        t = r.transform
        crs = r.crs or ref_crs
        if y == 2025:  # registration shift stated in the paper: 5 m east and 5 m south
            t = Affine(t.a, t.b, t.c + 5.0, t.d, t.e, t.f - 5.0)
        dst = np.full(ref_shape, np.nan)
        reproject(a, dst, src_transform=t, src_crs=crs, dst_transform=ref_t, dst_crs=ref_crs,
                  resampling=Resampling.bilinear, src_nodata=np.nan, dst_nodata=np.nan)
        maps[y] = dst * 10.0  # kg C m-2 -> Mg C ha-1

# ---------------- 30 m cells ----------------
H, W = ref_shape
h3, w3 = H // 3, W // 3
def cells(a):
    b = a[: h3 * 3, : w3 * 3].reshape(h3, 3, w3, 3)
    complete = np.isfinite(b).all(axis=(1, 3))
    m = np.where(complete, np.nanmean(b, axis=(1, 3)), np.nan)
    return m
C = {y: cells(maps[y]) for y in YEARS}
rows, cols = np.meshgrid(np.arange(h3), np.arange(w3), indexing="ij")
xs = ref_t.c + (cols * 3 + 1.5) * ref_t.a
ys = ref_t.f + (rows * 3 + 1.5) * ref_t.e

# ---------------- plots ----------------
p = pd.read_csv(f"{D}/plots.tab", sep=None, engine="python")
t = pd.read_csv(f"{D}/trees.tab", sep=None, engine="python")
hl = p[p.Plot.str.startswith("HL")].copy()
t["plot_id"] = t.Plot.str.rsplit("_", n=1).str[0]
yr = t[t.plot_id.str.startswith("HL")].groupby("plot_id").Year.max()
hl["meas_year"] = hl.Plot.map(yr)
res["hl_plots_total"] = int(len(hl))
res["hl_plots_by_meas_year"] = {str(k): int(v) for k, v in hl.meas_year.value_counts().sort_index().items()}
for c in ["AGB_2017_J", "AGB_2021_J", "AGB_2025_J", "AGB_2025_Y"]:
    if c in hl:
        hl[c] = hl[c] * 10.0

# ---------------- domain from plot hull ----------------
hull = MultiPoint(list(zip(hl.Easting, hl.Northing))).convex_hull
res["hull_area_ha"] = hull.area / 1e4
base_ok = np.isfinite(C[2017]) & np.isfinite(C[2025]) & np.isfinite(C[2012])
best = None
for buf in [-60, -30, 0, 30, 60, 90, 120]:
    g = hull.buffer(buf)
    inside = vectorized.contains(g, xs, ys)
    n = int((inside & base_ok).sum())
    res.setdefault("domain_candidates", []).append({"buffer_m": buf, "area_ha": g.area / 1e4, "cells": n})
    if best is None or abs(n - 1967) < abs(best[1] - 1967):
        best = (buf, n, inside)
buf, n_dom, inside = best
res["domain_used"] = {"buffer_m": buf, "cells_2017_2025_2012": n_dom, "paper_cells_2025": 1967}
dom = inside & base_ok
all6 = dom & np.all([np.isfinite(C[y]) for y in YEARS], axis=0)
c21 = dom & np.isfinite(C[2021])
res["cells_all_six"] = int(all6.sum()); res["cells_complete_2021"] = int(c21.sum())

V = {y: C[y][all6] for y in YEARS}
res["map_means_all6"] = {y: float(np.mean(V[y])) for y in YEARS}
res["rate_2012_2017"] = float(np.mean(V[2017] - V[2012]) / 5)
res["rate_2017_2021"] = float(np.mean(V[2021] - V[2017]) / 4)
yrs = np.array([2012, 2017, 2021.0])
Y = np.vstack([V[2012], V[2017], V[2021]]).T
xb = yrs.mean()
b = ((yrs - xb) * (Y - Y.mean(1, keepdims=True))).sum(1) / ((yrs - xb) ** 2).sum()
a0 = Y.mean(1) - b * xb
res["trend_rate"] = float(np.mean(b))
res["offsets_vs_trend"] = {y: float(np.mean(V[y] - (a0 + b * y))) for y in (2015, 2023, 2025)}
res["mean_change_2017_2025"] = float(np.mean(V[2025] - V[2017]))
for y0, y1 in [(2017, 2021), (2017, 2025), (2023, 2025)]:
    m = dom & np.isfinite(C[y0]) & np.isfinite(C[y1])
    res[f"loss_gt20pct_{y0}_{y1}"] = float(np.mean((C[y1][m] - C[y0][m]) / C[y0][m] < -0.2))

# ---------------- references and Method 2 ----------------
# constant increment on the 10 m pixels of the domain cells
pix_dom = np.zeros(ref_shape, bool)
pix_dom[: h3 * 3, : w3 * 3] = np.repeat(np.repeat(dom, 3, 0), 3, 1)
d10 = (maps[2017] - maps[2012])[pix_dom]
inc = float(np.nanmean(d10) / 5)
res["constant_increment_pixels"] = inc

# area-based curve on 90 m block means of 2012 and 2017 (3x3 cells)
h9, w9 = h3 // 3, w3 // 3
def blocks(a, mask):
    aa = np.where(mask, a, np.nan)[: h9 * 3, : w9 * 3].reshape(h9, 3, w9, 3)
    full = np.isfinite(aa).all(axis=(1, 3))
    return np.where(full, np.nanmean(aa, axis=(1, 3)), np.nan)
B12, B17 = blocks(C[2012], dom), blocks(C[2017], dom)
ok = np.isfinite(B12) & np.isfinite(B17)
res["area_based_pairs_90m"] = int(ok.sum())
def cr(a1, A, k, pp, dt):
    x = np.clip(a1 / A, 1e-9, None) ** (1 / pp)
    return A * (1 - (1 - x) * np.exp(-k * dt)) ** pp
fit = least_squares(lambda th: cr(B12[ok], *th, 5.0) - B17[ok], x0=[110, 0.04, 3], bounds=([75, 1e-4, 0.3], [1000, 1, 10]))
A, k, pp = fit.x
res["area_based_params"] = {"A": A, "k": k, "p": pp}

def score(pred, obs):
    e = pred - obs
    return {"bias": float(e.mean()), "rmse_pct": float(100 * np.sqrt((e ** 2).mean()) / obs.mean()),
            "r2": float(1 - (e ** 2).sum() / ((obs - obs.mean()) ** 2).sum())}
dom25 = dom & np.isfinite(C[2025])
preds = {}
for yv, dt, m in [(2021, 4, c21), (2025, 8, dom25)]:
    s17, obs = C[2017][m], C[yv][m]
    P = {"static": s17, "constinc": s17 + inc * dt, "2": cr(s17, A, k, pp, dt)}
    preds[yv] = (P, obs, m)
    res[f"scores_{yv}"] = {kk: score(v, obs) for kk, v in P.items()}
    res[f"cells_scored_{yv}"] = int(m.sum())
# offset-removed 2025
P, obs, m = preds[2025]
off = res["offsets_vs_trend"][2025]
res["scores_2025_offset_removed"] = {kk: score(v, obs - off) for kk, v in P.items()}

# ---------------- paired contrasts, block bootstrap ----------------
bb = []
for size in (90, 270, 450, 540):
    for yv in (2021, 2025):
        P, obs, m = preds[yv]
        bid = (np.floor(xs[m] / size).astype(int) * 100000 + np.floor(ys[m] / size).astype(int))
        ub = np.unique(bid); groups = [np.flatnonzero(bid == u) for u in ub]
        diffs = {("2", "static"): [], ("2", "constinc"): [], ("constinc", "static"): []}
        for _ in range(1000):
            idx = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
            o = obs[idx]; rm = {kk: 100 * np.sqrt(((v[idx] - o) ** 2).mean()) / o.mean() for kk, v in P.items()}
            for key in diffs:
                diffs[key].append(rm[key[0]] - rm[key[1]])
        for key, v in diffs.items():
            v = np.array(v)
            bb.append({"block_m": size, "n_blocks": len(ub), "year": yv, "contrast": f"{key[0]}-{key[1]}",
                       "mean": v.mean(), "lo95": np.percentile(v, 2.5), "hi95": np.percentile(v, 97.5)})
pd.DataFrame(bb).to_csv(f"{OUT}/contrasts.csv", index=False)

# ---------------- deciles binned on 2017, 2015, 2012 ----------------
dec = []
for base in (2017, 2015, 2012):
    m = all6
    q = pd.qcut(C[base][m], 10, labels=False)
    ch = C[2025][m] - C[2017][m]
    for d in range(10):
        dec.append({"binned_on": base, "decile": d + 1, "mean_base": float(C[base][m][q == d].mean()),
                    "obs_change_2017_2025": float(ch[q == d].mean())})
pd.DataFrame(dec).to_csv(f"{OUT}/deciles.csv", index=False)

# ---------------- map pair agreement and three-cornered hat ----------------
pairs = {}
for y0, y1 in [(2017, 2021), (2021, 2023), (2017, 2023), (2023, 2025)]:
    dlt = V[y1] - V[y0]
    pairs[f"{y0}_{y1}"] = float(100 * np.std(dlt) / np.mean(V[y1]))
res["map_pair_sd_pct"] = pairs
tch = {}
for tri in [(2012, 2015, 2017), (2015, 2017, 2021), (2017, 2021, 2023), (2021, 2023, 2025)]:
    a_, b_, c_ = (V[y] - V[y].mean() for y in tri)
    vab, vac, vbc = np.var(a_ - b_), np.var(a_ - c_), np.var(b_ - c_)
    est = {tri[0]: (vab + vac - vbc) / 2, tri[1]: (vab + vbc - vac) / 2, tri[2]: (vac + vbc - vab) / 2}
    tch["_".join(map(str, tri))] = {str(y): float(np.sqrt(v)) if v > 0 else float("nan") for y, v in est.items()}
res["three_cornered_hat_sd"] = tch

# ---------------- plots: map agreement and independent scoring ----------------
inv = ~ref_t
pr = []
for _, r0 in hl.iterrows():
    col, row = inv * (r0.Easting, r0.Northing)
    col, row = int(col), int(row)
    rec = {"plot": r0.Plot, "meas_year": r0.meas_year}
    for y in YEARS:
        w = maps[y][max(row - 1, 0): row + 2, max(col - 1, 0): col + 2]
        rec[f"map_{y}"] = float(np.nanmean(w)) if np.isfinite(w).any() else np.nan
    for c in ["AGB_2021_J", "AGB_2025_J", "AGB_2025_Y"]:
        rec[c] = r0.get(c, np.nan)
    pr.append(rec)
pr = pd.DataFrame(pr)
pr["pred_static_2025"] = pr.map_2017
pr["pred_constinc_2025"] = pr.map_2017 + inc * 8
pr["pred_2_2025"] = cr(pr.map_2017.to_numpy(), A, k, pp, 8)
pr.to_csv(f"{OUT}/plot_extract.csv", index=False)
ps = []
for grp, sel in [("first_measured_2026", pr.meas_year == 2026), ("measured_2024_2025", pr.meas_year.isin([2024, 2025])), ("all", pr.meas_year.notna())]:
    s = pr[sel & pr.AGB_2025_J.notna() & pr.map_2017.notna() & pr.map_2025.notna()]
    if len(s) < 3:
        continue
    obs = s.AGB_2025_J.to_numpy()
    for nm in ["map_2025", "pred_static_2025", "pred_constinc_2025", "pred_2_2025"]:
        pred = s[nm].to_numpy(); e = pred - obs
        bs = [e[rng.integers(0, len(e), len(e))].mean() for _ in range(2000)]
        ps.append({"plots": grp, "n": len(s), "prediction": nm, "plot_mean": obs.mean(), "bias": e.mean(),
                   "bias_lo95": np.percentile(bs, 2.5), "bias_hi95": np.percentile(bs, 97.5),
                   "rmse_pct": 100 * np.sqrt((e ** 2).mean()) / obs.mean(), "r": float(np.corrcoef(pred, obs)[0, 1])})
pd.DataFrame(ps).to_csv(f"{OUT}/plot_scores.csv", index=False)

json.dump(res, open(f"{OUT}/summary.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
print(pd.DataFrame(ps).round(2).to_string())
print(pd.DataFrame(bb).round(2).to_string())
print(pd.DataFrame(dec).round(1).to_string())
