from __future__ import annotations
import argparse, math, csv
from dataclasses import dataclass
from typing import Iterable, Tuple, List, Dict, Optional

import numpy as np

from runner import TrimodRunner

E3_EXP   = -8.482   # MeV (trizio reale) --> energia di legame del trizio, negativa per stato legato
SIGMA_TH = 0.5      # MeV (stima incertezza teorica: va poi sistemata)
SIGMA_EXP = 0.001    # MeV (incertezza sperimentale, trascurabile qui)
SIGMA    = math.hypot(SIGMA_TH, SIGMA_EXP)  # sigma totale in quadratura

PRIOR_RANGES = {             #prior uniforme (box) sui 4 parametri 
    "v01": (-200.0, -1.0), # MeV
    "v02": (-200.0, -1.0), # MeV
    "a1":  (0.2, 5.0), #fm
    "a2":  (0.2, 5.0), #fm
}

def logprior(theta: Iterable[float]) -> float: 
    """
    prior uniforme: 0 se dentro il box, -inf se fuori
    """
    v01,v02,a1,a2 = theta
    if not(PRIOR_RANGES["v01"][0] <= v01 <= PRIOR_RANGES["v01"][1]): return -math.inf
    if not(PRIOR_RANGES["v02"][0] <= v02 <= PRIOR_RANGES["v02"][1]): return -math.inf
    if not(PRIOR_RANGES["a1"][0]  <= a1  <= PRIOR_RANGES["a1"][1]):  return -math.inf
    if not(PRIOR_RANGES["a2"][0]  <= a2  <= PRIOR_RANGES["a2"][1]):  return -math.inf
    return 0.0

def loglike(theta: Iterable[float], tr: TrimodRunner) -> float:  #likelihood gaussiana su E3
    """
    calcola la log-likelihood gaussiana:
      L ∝ exp( - (E3(theta)-E3_exp)^2 / (2 sigma^2) )
     e restituisce la versione log con la costante di normalizzazione
    """
    E3 = tr.run_one(theta) #lancia trimod e legge E3
    r  = (E3 - E3_EXP) / SIGMA
    return -0.5 * (r*r) - math.log(SIGMA * math.sqrt(2.0*math.pi))

def logposterior(theta: Iterable[float], tr: TrimodRunner) -> float:
    """
    logposterior = logprior + loglike.
    Se TRIMOD non converge o lancia eccezioni, restituisce -inf (posterior ~ 0)
    """
    lp = logprior(theta)
    if not np.isfinite(lp):
        return -math.inf
    try:
        return lp + loglike(theta, tr)
    except Exception:
        # non-converged run ⇒ posterior ~ 0
        return -math.inf

def logsumexp(a):
    """
    log( sum_i exp(a_i) ) 
    """
    a = np.asarray(a, dtype=float)
    m = np.nanmax(a)
    if not np.isfinite(m):
        return -np.inf
    return m + np.log(np.nansum(np.exp(a - m)))
    
def sweep1d(param: str, pmin: float, pmax: float, n: int,
            theta_fixed: Tuple[float, float, float, float],
            tr: TrimodRunner, csv_out: str | None = None) -> List[Dict[str, float]]:
    """
    Valuta logposterior su una griglia 1D del parametro `param` ∈ {v01,v02,a1,a2}
    Gli altri 3 parametri restano fissi
    """
    idx = {"v01":0,"v02":1,"a1":2,"a2":3}[param] #indice del parametro da variare
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

    if args.theta:  #valutazione singola del posterior
        th = tuple(args.theta)
        lp = logposterior(th, tr)
        print(f"logposterior({th}) = {lp}")
        return

    if args.sweep1d:  #sweep 1D
        param, pmin, pmax, n = args.sweep1d
        if param not in ("v01","v02","a1","a2"):
            raise SystemExit("param must be one of: v01 v02 a1 a2")
        rows = sweep1d(param, float(pmin), float(pmax), int(n),
                        theta_fixed=tuple(args.theta_fixed), tr=tr,
                        csv_out=args.csv)
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






