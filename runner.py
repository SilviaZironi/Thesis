from __future__ import annotations
from pathlib import Path
import re
import subprocess
import sys

# regex per estrarre la riga corrispodente all'energia di legame
_RE_E = re.compile(r"3-body binding energy\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*MeV")

def write_params(theta, workdir="."):
    """
    theta = (v01, v02, a1, a2)
    scrive params.in come UNA riga: v01 v02 a1 a2
    """
    v01, v02, a1, a2 = theta
    Path(workdir, "params.in").write_text(f"{v01} {v02} {a1} {a2}\n", encoding="utf-8")


def run_trimod(exe: str | Path = "./trimod.x", workdir: str | Path | None = None) -> float:
    """
    lancia l'eseguibile Fortran (trimod) e restituisce la 3-body binding energy in MeV.
    - exe: path dell'eseguibile (./trimod.x oppure ./run)
    - workdir: directory di lavoro dove si trova output.dat (default: cartella dell'eseguibile)
    """
    exe = Path(exe).resolve()
    if not exe.exists():
        raise FileNotFoundError(f"Eseguibile non trovato: {exe}")

    wd = Path(workdir) if workdir else exe.parent
    wd = wd.resolve()
    if not wd.exists():
        raise FileNotFoundError(f"Working directory non trovata: {wd}")

    # lancia il programma Fortran nella sua cartella
    subprocess.run([str(exe)], cwd=wd, check=True)

    # Legge il file di output
    out = (wd / "output.dat").read_text(encoding="utf-8")

    # Estrae l'energia
    m = _RE_E.search(out)
    if not m:
        raise RuntimeError("Non ho trovato la riga con '3-body binding energy' in output.dat")
    return float(m.group(1))



