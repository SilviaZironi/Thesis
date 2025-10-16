from __future__ import annotations
from pathlib import Path
import os, re, subprocess
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

# regex per estrarre la riga corrispodente all'energia di legame
_RE_E = re.compile(r"3-body\s+binding\s+energy\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*MeV")

@dataclass
class TrimodRunner: # WRAPPER per trimod: gestisce dove si trova l'eseguibile (EXE) e dove lavorare/salvare i file di I/O (WORKDIR)
    exe: str | os.PathLike = "./run" #path all'eseguibile fortran
    workdir: str | os.PathLike = "./WORK" #directory di lavoro

    def __post_init__(self):
     self.exe_path = Path(self.exe)  #path assoluto dell'eseguibile
     self.exe_path = self.exe_path.resolve()  #trasforma path in un percorso assoluto, così da poter cambiare cartella di
                                              #lavoro senza perdere l'eseguibile

    #workdir: converte workdir in un oggetto Path
     self.workdir = Path(self.workdir)
     self.workdir.mkdir(parents=True, exist_ok=True) #creao la cartella se non esiste e non da errore se già c'è

     if not self.exe_path.exists():  #solo per sicurezza, se non trova l'eseguibile dove gli è stato detto, prova con ./fortran/run, ma abbastanza inutile 
        alt = Path("./fortran/run").resolve()
        if alt.exists():
            self.exe_path = alt

     if not self.exe_path.exists(): #check: se ancora non lo trova si ferma e da errore
        raise FileNotFoundError(f"TRIMOD executable not found: {self.exe}")

     try:  #prova ad aggiungere il permesso di esecuzione (equivalente di "chmod +x run" in shell) al file run, ma se non riesce va comunque avanti
         self.exe_path.chmod(self.exe_path.stat().st_mode | 0o111)
     except Exception:
        pass


    def write_params(self, theta: Iterable[float]) -> Path: #scrive i parametri del potenziale in params.in come riga
        v01, v02, a1, a2 = map(float, tuple(theta))
        p = self.workdir / "params.in"
        p.write_text(f"{v01:.12g} {v02:.12g} {a1:.12g} {a2:.12g}\n", encoding="utf-8")
        return p


    def _clean_output(self) -> None: #per avere run puliti da precedenti output
        for name in ("output.dat", "stdout.txt", "stderr.txt"):
            f = self.workdir / name
            if f.exists():
                try:
                    f.unlink()
                except Exception:
                    pass

    @staticmethod
    def _parse_energy(text: str): #estrae l'energia dal testo di output.dat usando la regex sopra
        m = _RE_E.search(text)
        return float(m.group(1)) if m else None

    def run_once(self): #esegue trimod una volta in WORKDIR
        self._clean_output()   #pulizia (vedi prima)
        proc = subprocess.run([str(self.exe_path)], cwd=str(self.workdir),   #esegue l'eseguibile run
                              capture_output=True, text=True)
        (self.workdir / "stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
        (self.workdir / "stderr.txt").write_text(proc.stderr or "", encoding="utf-8")    #scrive quello che fortran ha stampato
        out_text = ""
        out_file = self.workdir / "output.dat"
        if out_file.exists():    #lettura del file se esiste
            try:
                out_text = out_file.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                out_text = out_file.read_bytes().decode("utf-8", errors="ignore")
        combined = (proc.stdout or "") + "\n" + (proc.stderr or "") + "\n" + out_text
        lowered = combined.lower()
        nonconv = any(k in lowered for k in ["too many iterations", "itermax", "iemax"])
        E3 = self._parse_energy(out_text)   #chiama la funzione sopra per estrarre l'energia
        if E3 is None or nonconv:
            return (None, combined)
        return (E3, combined)

    def run_one(self, theta: Iterable[float]) -> float: #scrive params.in, esegue trimod, ritorna E3
        self.write_params(theta)     #scrive params.in (v01 v02 a1 a2) nella WORKDIR
        E3, log = self.run_once()    #lancia l’eseguibile Fortran e legge output/log
        if E3 is None:               #se non trova E3 stampa le ultime 25 righe 
            last = "\n".join((log or "").splitlines()[-25:])
            raise RuntimeError("TRIMOD not converged or energy not found.\n--- tail ---\n" + last)
        return E3


#per usare lo script da terminale (quindi theta posso inserirlo da terminale)
if __name__ == "__main__":
    import argparse, sys
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./run")
    ap.add_argument("--workdir", default="./WORK")
    ap.add_argument("--theta", nargs=4, type=float)
    a = ap.parse_args()
    tr = TrimodRunner(exe=a.exe, workdir=a.workdir)
    th = tuple(a.theta) if a.theta else (-50.0, -30.0, 1.5, 2.0)  #theta default se non gli viene fornito da terminale
    try:
        print(f"{tr.run_one(th):.6f}")
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)


