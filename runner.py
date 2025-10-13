from __future__ import annotations
from pathlib import Path
import os, re, subprocess
from dataclasses import dataclass
from typing import Iterable, Optional, Tuple

# regex per estrarre la riga corrispodente all'energia di legame
_RE_E = re.compile(r"3-body\s+binding\s+energy\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s*MeV")

@dataclass
class TrimodRunner: #wrapper per trimod
    exe: str | os.PathLike = "./run" #path all'eseguibile fortran
    workdir: str | os.PathLike = "./WORK" #directory di lavoro

    def __post_init__(self):
     self.exe_path = Path(self.exe)  #path assoluto dell'eseguibile
     self.exe_path = self.exe_path.resolve()  

    #workdir
     self.workdir = Path(self.workdir)
     self.workdir.mkdir(parents=True, exist_ok=True)

     if not self.exe_path.exists():
        alt = Path("./fortran/run").resolve()
        if alt.exists():
            self.exe_path = alt

     if not self.exe_path.exists(): #check
        raise FileNotFoundError(f"TRIMOD executable not found: {self.exe}")

     try:
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
    def _parse_energy(text: str) -> Optional[float]: #estrae l'energia dal testo di output.dat usando la regex sopra
        m = _RE_E.search(text)
        return float(m.group(1)) if m else None

    def run_once(self) -> Tuple[Optional[float], str]: #esegue trimod una volta in WORKDIR
        self._clean_output()
        proc = subprocess.run([str(self.exe_path)], cwd=str(self.workdir),
                              capture_output=True, text=True)
        (self.workdir / "stdout.txt").write_text(proc.stdout or "", encoding="utf-8")
        (self.workdir / "stderr.txt").write_text(proc.stderr or "", encoding="utf-8")
        out_text = ""
        out_file = self.workdir / "output.dat"
        if out_file.exists():
            try:
                out_text = out_file.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                out_text = out_file.read_bytes().decode("utf-8", errors="ignore")
        combined = (proc.stdout or "") + "\n" + (proc.stderr or "") + "\n" + out_text
        lowered = combined.lower()
        nonconv = any(k in lowered for k in ["too many iterations", "itermax", "iemax"])
        E3 = self._parse_energy(out_text)
        if E3 is None or nonconv:
            return (None, combined)
        return (E3, combined)

    def run_one(self, theta: Iterable[float]) -> float: #scrive params.in, esegue trimod, ritorna E3
        self.write_params(theta)
        E3, log = self.run_once()
        if E3 is None:
            last = "\n".join((log or "").splitlines()[-25:])
            raise RuntimeError("TRIMOD not converged or energy not found.\n--- tail ---\n" + last)
        return E3

if __name__ == "__main__":
    import argparse, sys
    ap = argparse.ArgumentParser()
    ap.add_argument("--exe", default="./run")
    ap.add_argument("--workdir", default="./WORK")
    ap.add_argument("--theta", nargs=4, type=float)
    a = ap.parse_args()
    tr = TrimodRunner(exe=a.exe, workdir=a.workdir)
    th = tuple(a.theta) if a.theta else (-50.0, -30.0, 1.5, 2.0)
    try:
        print(f"{tr.run_one(th):.6f}")
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

