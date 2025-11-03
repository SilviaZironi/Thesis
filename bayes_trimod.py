"""
bayes_trimod.py — Minimal Bayesian layer on top of TRIMOD, using runner.py
- box prior on (v01, v02, a1, a2)
- Gaussian likelihood on E3 with sigma (theory) + optional sigma_exp
- posterior = prior + likelihood in log space
- includes a simple 1D grid sweep utility (for quick plots)

Run examples:
    python bayes_trimod.py --exe ./run --workdir ./WORK --sweep1d v01 -120 -20 12
    python bayes_trimod.py --exe ./run --workdir ./WORK --theta -50 -30 1.5 2.0
"""
from __future__ import annotations
import argparse, math, csv
from dataclasses import dataclass
from typing import Iterable, Tuple, List, Dict, Optional

import numpy as np

def logsumexp(a):
    a = np.asarray(a, dtype=float)
    m = np.nanmax(a)
    if not np.isfinite(m):
        return -np.inf
    return m + np.log(np.nansum(np.exp(a - m)))

from runner import TrimodRunner

# --- Reference triton energy and uncertainties ---
E3_EXP   = -8.482   # MeV (experimental)
SIGMA_TH = 0.5      # MeV (theory/model uncertainty) — adjust with your advisor
SIGMA_EXP= 0.001    # MeV (experimental; tiny, often negligible)
SIGMA    = math.hypot(SIGMA_TH, SIGMA_EXP)

# --- Prior ranges (box) — adjust ranges to sensible physics ---
PRIOR_RANGES = {
    'v01': (-620.0, -530.0),
    'v02': (1200.0, 1600.0),
    'a1' : (1.40, 1.70),
    'a2' : (2.80, 3.40),
}

def logprior(theta: Iterable[float]) -> float:
    v01,v02,a1,a2 = theta
    if not(PRIOR_RANGES["v01"][0] <= v01 <= PRIOR_RANGES["v01"][1]): return -math.inf
    if not(PRIOR_RANGES["v02"][0] <= v02 <= PRIOR_RANGES["v02"][1]): return -math.inf
    if not(PRIOR_RANGES["a1"][0]  <= a1  <= PRIOR_RANGES["a1"][1]):  return -math.inf
    if not(PRIOR_RANGES["a2"][0]  <= a2  <= PRIOR_RANGES["a2"][1]):  return -math.inf
    # uniform box prior: constant (drop additive constant)
    return 0.0

def loglike(theta: Iterable[float], tr: TrimodRunner) -> float:
    # Write params, run TRIMOD, read E3
    E3 = tr.run_one(theta)
    r  = (E3 - E3_EXP) / SIGMA
    return -0.5 * (r*r) - math.log(SIGMA * math.sqrt(2.0*math.pi))

def logposterior(theta: Iterable[float], tr: TrimodRunner) -> float:
    lp = logprior(theta)
    if not np.isfinite(lp):
        return -math.inf
    try:
        return lp + loglike(theta, tr)
    except Exception:
        # non-converged run ⇒ posterior ~ 0
        return -math.inf

# --- CLI utilities ---

def sweep1d(param: str, pmin: float, pmax: float, n: int,
            theta_fixed: Tuple[float, float, float, float],
            tr: TrimodRunner, csv_out: str | None = None) -> List[Dict[str, float]]:
    idx = {"v01":0,"v02":1,"a1":2,"a2":3}[param]
    grid = np.linspace(pmin, pmax, n)
    rows: List[Dict[str,float]] = []
    for val in grid:
        theta = list(theta_fixed); theta[idx] = float(val)
        lp = logposterior(theta, tr)
        if np.isfinite(lp):
            try: E3 = tr.run_one(theta)
            except Exception: E3 = float("nan")
        else:
            E3 = float("nan")
        rows.append({"param": param, "x": float(val), "E3": E3, "logpost": float(lp)})

    logw = np.array([r["logpost"] for r in rows], dtype=float)
    mask = np.isfinite(logw)
    if mask.any():
        logw_norm = logw.copy()
        logw_norm[mask] -= logsumexp(logw[mask])
        w = np.exp(logw_norm)
    else:
        w = np.zeros_like(logw)

    for r, wi in zip(rows, w):
        r["post_norm"] = float(wi)

    if csv_out:
        with open(csv_out, "w", newline="") as f:
            wtr = csv.DictWriter(f, fieldnames=["param","x","E3","logpost","post_norm"])
            wtr.writeheader(); wtr.writerows(rows)

    return rows


def main():
    ap = argparse.ArgumentParser(description="Bayesian layer for TRIMOD (Fortran)")
    ap.add_argument("--exe", default="./run")
    ap.add_argument("--workdir", default="./WORK")
    ap.add_argument("--theta", nargs=4, type=float,
                    help="Theta=(v01 v02 a1 a2) to evaluate a single posterior.")
    ap.add_argument("--sweep1d", nargs=4, metavar=("param","min","max","n"),
                    help="Do a 1D sweep over 'param' with fixed others. Example: v01 -120 -20 12")
    ap.add_argument("--theta-fixed", nargs=4, type=float, default=[-50.0,-30.0,1.5,2.0],
                    help="Fixed theta for sweep1d (defaults: -50 -30 1.5 2.0)")
    ap.add_argument("--csv", default="sweep1d.csv")
    args = ap.parse_args()

    tr = TrimodRunner(exe=args.exe, workdir=args.workdir)

    if args.theta:
        th = tuple(args.theta)
        lp = logposterior(th, tr)
        print(f"logposterior({th}) = {lp}")
        return

    if args.sweep1d:
        param, pmin, pmax, n = args.sweep1d
        if param not in ("v01","v02","a1","a2"):
            raise SystemExit("param must be one of: v01 v02 a1 a2")
        rows = sweep1d(param, float(pmin), float(pmax), int(n),
                        theta_fixed=tuple(args.theta_fixed), tr=tr,
                        csv_out=args.csv)
        # print short summary
        xs = np.array([r["x"] for r in rows])
        e3 = np.array([r["E3"] for r in rows], dtype=float)
        post = np.array([r["post_norm"] for r in rows], dtype=float)
        print(f"Saved {len(rows)} rows to {args.csv}.")
        print("x[min,max] =", float(np.nanmin(xs)), float(np.nanmax(xs)))
        print("E3[finite] count =", int(np.isfinite(e3).sum()))
        print("posterior sum ~", float(post.sum()))
        return

    ap.print_help()

if __name__ == "__main__":
    main()
