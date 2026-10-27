#!/usr/bin/env python3
"""Howland AGB projection: reanalysis checks requested in the October 22, 2026 red team.

Runs on the per-cell and per-plot tables of the current (version 2 map, 70 plot) analysis.
Inputs are plain CSV files so the script is independent of the original pipeline.

Checks
  1. Offset-adjusted rescoring of 2025 (offset estimated per cell from the G-LiHT trend, or fixed).
  2. Block bootstrap at several block sizes with block counts, paired RMSE contrasts against the
     static and constant-increment references, Holm adjustment across the contrast family, and
     Monte Carlo error of the interval endpoints. Optional refit of the area-based curve
     (Method 2) and the constant increment inside every replicate.
  3. Residual by decile binned on an earlier map (2012 or 2015) as well as 2017, with intervals,
     to separate the decile gradient from regression to the mean.
  4. Rescoring without disturbance-flagged cells and, if supplied, without an independent mask.
  5. Plot-level scoring of every method at plots that did not calibrate the maps.

Usage
  python howland_reanalysis_checks.py --cells cells.csv --plots plots.csv --out results/ \
      [--offset 4.8] [--blocks 90 270 450 540] [--reps 2000] [--refit] [--seed 20261023]

cells.csv columns (one row per scored 30 m cell; AGB in Mg C ha-1):
  cell_id, x, y                         projected coordinates of the cell center (m)
  map_2012, map_2015, map_2017, map_2021, map_2023, map_2025
  complete_2021                         1 if the cell has a complete 2021 map, else 0
  disturbed                             optional, 1 if flagged by the >20% loss screen
  indep_mask                            optional, 1 if disturbed by independent data (Landsat, imagery)
  pred_<method>_2021, pred_<method>_2025  one pair per method, for example pred_2_2021
  The static and constant-increment references are built by the script; do not supply them.
plots.csv columns:
  plot_id, cell_id, agb_2021, agb_2025, calibrated   calibrated = 1 if the plot calibrated the maps
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

RNG_DEFAULT = 20261023


def cr_pred(agb1, A, k, p, dt):
    """Algebraic difference Chapman-Richards (Equation 1)."""
    x = np.clip(agb1 / A, 1e-9, None) ** (1.0 / p)
    return A * (1.0 - (1.0 - x) * np.exp(-k * dt)) ** p


def cr_fit(agb1, agb2, dt, start=(110.0, 0.04, 3.0)):
    from scipy.optimize import least_squares
    def res(th):
        return cr_pred(agb1, th[0], th[1], th[2], dt) - agb2
    r = least_squares(res, x0=np.array(start), bounds=([75, 1e-4, 0.3], [1000, 1.0, 10.0]))
    return r.x


def score(pred, obs):
    e = pred - obs
    return {
        "bias": float(np.mean(e)),
        "rmse": float(np.sqrt(np.mean(e ** 2))),
        "rmse_pct": float(100 * np.sqrt(np.mean(e ** 2)) / np.mean(obs)),
    }


def methods_of(df):
    m = sorted({c.split("_")[1] for c in df.columns if c.startswith("pred_")})
    return m


def build_refs(df, ci_rate=None):
    rate = ci_rate if ci_rate is not None else float(np.mean((df.map_2017 - df.map_2012) / 5.0))
    out = df.copy()
    for yr, dt in ((2021, 4), (2025, 8)):
        out[f"pred_static_{yr}"] = out.map_2017
        out[f"pred_constinc_{yr}"] = out.map_2017 + rate * dt
    return out, rate


def offset_from_trend(df):
    """Per-cell linear trend through the three G-LiHT years, evaluated at 2025."""
    yrs = np.array([2012, 2017, 2021], float)
    Y = df[["map_2012", "map_2017", "map_2021"]].to_numpy(float)
    xb = yrs.mean()
    b = ((yrs - xb) * (Y - Y.mean(1, keepdims=True))).sum(1) / ((yrs - xb) ** 2).sum()
    a = Y.mean(1) - b * xb
    trend25 = a + b * 2025
    return df.map_2025.to_numpy(float) - trend25


def table_scores(df, methods, label):
    rows = []
    for yr in (2021, 2025):
        sub = df if yr == 2025 else df[df.complete_2021 == 1]
        obs = sub[f"obs_{yr}"].to_numpy(float)
        for m in methods:
            s = score(sub[f"pred_{m}_{yr}"].to_numpy(float), obs)
            rows.append({"set": label, "year": yr, "method": m, "n_cells": len(sub), **s})
    return pd.DataFrame(rows)


def block_ids(df, size):
    return (np.floor(df.x / size).astype(int).astype(str) + "_" + np.floor(df.y / size).astype(int).astype(str)).to_numpy()


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty_like(p)
    run = 0.0
    for i, idx in enumerate(o):
        run = max(run, (len(p) - i) * p[idx])
        adj[idx] = min(1.0, run)
    return adj


def bootstrap(df, methods, size, reps, rng, refit=False):
    ids = block_ids(df, size)
    ub = np.unique(ids)
    groups = [np.flatnonzero(ids == b) for b in ub]
    refs = ["static", "constinc"]
    stats = {}
    for r in range(reps):
        pick = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[i] for i in pick])
        sub = df.iloc[idx]
        if refit:
            rate = float(np.mean((sub.map_2017 - sub.map_2012) / 5.0))
            A, k, p = cr_fit(sub.map_2012.to_numpy(float), sub.map_2017.to_numpy(float), 5.0)
        for yr, dt in ((2021, 4), (2025, 8)):
            s2 = sub if yr == 2025 else sub[sub.complete_2021 == 1]
            obs = s2[f"obs_{yr}"].to_numpy(float)
            rm = {}
            for m in methods + refs:
                pred = s2[f"pred_{m}_{yr}"].to_numpy(float)
                if refit and m == "constinc":
                    pred = s2.map_2017.to_numpy(float) + rate * dt
                if refit and m == "2":
                    pred = cr_pred(s2.map_2017.to_numpy(float), A, k, p, dt)
                rm[m] = 100 * np.sqrt(np.mean((pred - obs) ** 2)) / np.mean(obs)
            for m in methods:
                for ref in refs:
                    stats.setdefault((yr, m, ref), []).append(rm[m] - rm[ref])
    rows = []
    for (yr, m, ref), v in stats.items():
        v = np.asarray(v)
        lo, hi = np.percentile(v, [2.5, 97.5])
        # Monte Carlo error of the endpoints from 10 batches
        bat = np.array_split(rng.permutation(v), 10)
        se_lo = np.std([np.percentile(b, 2.5) for b in bat]) / np.sqrt(10)
        se_hi = np.std([np.percentile(b, 97.5) for b in bat]) / np.sqrt(10)
        p2 = 2 * min(np.mean(v <= 0), np.mean(v >= 0))
        rows.append({"block_m": size, "n_blocks": len(ub), "year": yr, "method": m, "reference": ref,
                     "mean_diff_pp": v.mean(), "lo95": lo, "hi95": hi, "mc_se_lo": se_lo,
                     "mc_se_hi": se_hi, "p_boot": max(p2, 1.0 / len(v)), "refit": refit})
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out.p_boot.to_numpy())
    return out


def deciles(df, methods, base, reps, rng):
    rows = []
    q = pd.qcut(df[base], 10, labels=False, duplicates="drop")
    for yr in (2021, 2025):
        sub = df if yr == 2025 else df[df.complete_2021 == 1]
        qq = q.loc[sub.index]
        for m in methods + ["static", "constinc"]:
            res = sub[f"pred_{m}_{yr}"] - sub[f"obs_{yr}"]
            for d in sorted(qq.unique()):
                v = res[qq == d].to_numpy(float)
                bs = [rng.choice(v, len(v)).mean() for _ in range(min(reps, 500))]
                rows.append({"binned_on": base, "year": yr, "method": m, "decile": int(d) + 1,
                             "mean_base": float(sub[base][qq == d].mean()), "mean_resid": v.mean(),
                             "lo95": np.percentile(bs, 2.5), "hi95": np.percentile(bs, 97.5)})
        ch = (sub[f"obs_{yr}"] - sub.map_2017)
        for d in sorted(qq.unique()):
            rows.append({"binned_on": base, "year": yr, "method": "observed_change", "decile": int(d) + 1,
                         "mean_base": float(sub[base][qq == d].mean()), "mean_resid": float(ch[qq == d].mean()),
                         "lo95": np.nan, "hi95": np.nan})
    return pd.DataFrame(rows)


def plot_scores(plots, cells, methods, reps, rng):
    m = plots.merge(cells, on="cell_id", how="left", suffixes=("", "_cell"))
    rows = []
    for subset, sel in (("independent", m.calibrated == 0), ("calibration", m.calibrated == 1), ("all", m.calibrated >= 0)):
        sub = m[sel]
        if len(sub) == 0:
            continue
        for yr in (2021, 2025):
            obs = sub[f"agb_{yr}"].to_numpy(float)
            for meth in methods + ["static", "constinc"]:
                pred = sub[f"pred_{meth}_{yr}"].to_numpy(float)
                ok = np.isfinite(pred) & np.isfinite(obs)
                if ok.sum() < 3:
                    continue
                s = score(pred[ok], obs[ok])
                bs = []
                for _ in range(min(reps, 1000)):
                    i = rng.integers(0, ok.sum(), ok.sum())
                    bs.append(np.mean(pred[ok][i] - obs[ok][i]))
                rows.append({"plots": subset, "year": yr, "method": meth, "n": int(ok.sum()), **s,
                             "bias_lo95": np.percentile(bs, 2.5), "bias_hi95": np.percentile(bs, 97.5)})
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True)
    ap.add_argument("--plots")
    ap.add_argument("--out", required=True)
    ap.add_argument("--offset", type=float, default=None, help="fixed 2025 offset; default per-cell trend")
    ap.add_argument("--blocks", type=float, nargs="+", default=[90, 270, 450, 540])
    ap.add_argument("--reps", type=int, default=2000)
    ap.add_argument("--refit", action="store_true")
    ap.add_argument("--seed", type=int, default=RNG_DEFAULT)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    log = open(os.path.join(a.out, "error_log.txt"), "a")
    try:
        cells = pd.read_csv(a.cells)
    except Exception as e:  # noqa: BLE001
        log.write(f"read cells: {e}\n"); sys.exit(1)
    methods = methods_of(cells)
    cells, rate = build_refs(cells)
    cells["obs_2021"] = cells.map_2021
    cells["obs_2025"] = cells.map_2025
    allm = methods + ["static", "constinc"]

    # 1. offset-adjusted rescoring
    off_cell = offset_from_trend(cells)
    off = a.offset if a.offset is not None else float(np.mean(off_cell))
    base = table_scores(cells, allm, "as_published")
    adj = cells.copy(); adj["obs_2025"] = adj.map_2025 - off
    t1 = pd.concat([base, table_scores(adj, allm, f"offset_removed_{off:.2f}")])
    t1.to_csv(os.path.join(a.out, "01_offset_rescoring.csv"), index=False)

    # 2. block bootstrap contrasts
    bb = pd.concat([bootstrap(cells, methods, b, a.reps, rng, refit=False) for b in a.blocks])
    if a.refit:
        bb = pd.concat([bb] + [bootstrap(cells, methods, b, a.reps, rng, refit=True) for b in a.blocks])
    bb.to_csv(os.path.join(a.out, "02_bootstrap_contrasts.csv"), index=False)

    # 3. deciles on 2017, 2015, 2012
    dd = pd.concat([deciles(cells, methods, b, a.reps, rng) for b in ("map_2017", "map_2015", "map_2012")])
    dd.to_csv(os.path.join(a.out, "03_deciles_by_base_map.csv"), index=False)

    # 4. disturbance sensitivity
    parts = [table_scores(cells, allm, "all_cells")]
    for col in ("disturbed", "indep_mask"):
        if col in cells:
            parts.append(table_scores(cells[cells[col] != 1], allm, f"without_{col}"))
    pd.concat(parts).to_csv(os.path.join(a.out, "04_disturbance_sensitivity.csv"), index=False)

    # 5. plot-level scoring
    if a.plots:
        plots = pd.read_csv(a.plots)
        plot_scores(plots, cells, methods, a.reps, rng).to_csv(os.path.join(a.out, "05_plot_level_scores.csv"), index=False)

    with open(os.path.join(a.out, "00_summary.md"), "w") as f:
        f.write(f"Seed {a.seed}; {len(cells)} cells; methods {', '.join(methods)}; constant increment {rate:.3f} Mg C ha-1 yr-1.\n")
        f.write(f"2025 offset used {off:.2f} Mg C ha-1 (per-cell trend mean {np.mean(off_cell):.2f}).\n")
        nb = bb.groupby("block_m").n_blocks.first().to_dict()
        f.write(f"Blocks per size: {nb}.\n")
        f.write("Tables 01 to 05 are written beside this file. Interpret p_holm across the whole contrast family.\n")
    log.close()


if __name__ == "__main__":
    main()
