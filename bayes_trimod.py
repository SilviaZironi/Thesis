from __future__ import annotations
import numpy as np
from runner import write_params, run_trimod

E3_exp   = -8.482   # MeV (trizio reale) --> energia di legame del trizio, negativa per stato legato
sigma_th = 0.5      # MeV (stima incertezza teorica: va poi sistemata)
sigma_exp= 0.001    # MeV (incertezza sperimentale, trascurabile qui)

def logprior(theta) -> float:
    """
    prior semplice su range (controlla che i parametri del potenziale siano ragionevoli)
    -> poi eventualmente si possono sistemare per renderla migliore
    """
    v01, v02, a1, a2 = theta
    if not (-1500 < v01 <   0): return -np.inf
    if not (    0 < v02 < 5000): return -np.inf
    if not (  0.5 <  a1 <  3.5): return -np.inf
    if not (  1.0 <  a2 <  6.0): return -np.inf
    return 0.0

def loglike(theta, exe="./trimod.x", workdir=".") -> float:
    """
    likelihood gaussiana: confronta il valore di E3 calcolata da trimod e quella sperimentale
    """
    write_params(theta, workdir)                 # 1) scrivi params.in per fortran
    val = run_trimod(exe=exe, workdir=workdir)
    E3 = float(val[0]) if isinstance(val, tuple) else float(val) # 2) lancia trimod e legge E3(theta)
    if not np.isfinite(E3):
        raise RuntimeError("E3 non finita (NaN/inf)")

    sigma = float(np.hypot(sigma_th, sigma_exp))
    r = (E3 - E3_exp) / sigma
    return -0.5 * (r * r) - np.log(sigma * np.sqrt(2.0 * np.pi))

def logposterior(theta, exe="./trimod.x", workdir=".") -> float:
    lp = logprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    try:
        return lp + loglike(theta, exe=exe, workdir=workdir)
    except Exception as e:
        # Se TRIMOD non converge o fallisce → punto “impossibile”
        print(f"[logposterior] eccezione: {type(e).__name__}: {e}")
        return -np.inf
