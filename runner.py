from __future__ import annotations
import os, re, subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Tuple

# accetta sia "3-body binding energy = ..." sia "emev = ..."
_RE_E = re.compile(
    r"(?:3-body\s+binding\s+energy|emev)\s*=\s*([-+]?[\d.]+(?:[eE][-+]?\d+)?)",
    re.IGNORECASE,
)

@dataclass
class TrimodRunner:
    exe: str | os.PathLike = "./run"   
    workdir: str | os.PathLike = "./WORK"

    def __post_init__(self):
        # path eseguibile
        self.exe_path = Path(self.exe).resolve()
        # workdir
        self.workdir = Path(self.workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)

        # fallback automatico 
        if not self.exe_path.exists():
            alt = Path("./fortran/run").resolve() 
            if alt.exists():
                self.exe_path = alt

        if not self.exe_path.exists():
            raise FileNotFoundError(f"TRIMOD executable not found: {self.exe}")

        try:
            self.exe_path.chmod(self.exe_path.stat().st_mode | 0o111)
        except Exception:
            pass


    def write_params(self, theta: Iterable[float]) -> Path:
        v01, v02, a1, a2 = map(float, tuple(theta))
        p = self.workdir / "params.in"
        p.write_text(f"{v01:.12g} {v02:.12g} {a1:.12g} {a2:.12g}\n", encoding="utf-8")
        return p

    def _clean_output(self) -> None:
        for name in ("output.dat", "stdout.txt", "stderr.txt", "energy.out"):
            f = self.workdir / name
            if f.exists():
                try: f.unlink()
                except Exception: pass

    @staticmethod
    def _parse_energy(text: str) -> Optional[float]:
        vals = [float(m.group(1)) for m in _RE_E.finditer(text)]
        return vals[-1] if vals else None

    # accetto theta e la passo via CLI
    def run_once(self, theta: Optional[Iterable[float]] = None) -> Tuple[Optional[float], str]:
        self._clean_output()

        args = []
        if theta is not None:
            v01, v02, a1, a2 = map(float, tuple(theta))
            args = [f"{v01:.12g}", f"{v02:.12g}", f"{a1:.12g}", f"{a2:.12g}"]

        proc = subprocess.run(
            [str(self.exe_path), *args],
            cwd=str(self.workdir),
            capture_output=True,
            text=True,
            check=False,
        )

        combined = (proc.stdout or "") + "\n" + (proc.stderr or "")

        # fallback: aggrega eventuale output.dat se esiste
        out_file = self.workdir / "output.dat"
        if out_file.exists():
            try:
                combined += "\n" + out_file.read_text(encoding="utf-8", errors="ignore")
            except Exception:
                combined += "\n" + out_file.read_bytes().decode("utf-8", errors="ignore")

        E3 = self._parse_energy(combined)
        low = combined.lower()
        nonconv = any(k in low for k in ("too many iterations", "itermax", "iemax", "not converg"))

        if E3 is None or nonconv:
            return (None, combined)
        return (E3, combined)

    def run_one(self, theta: Iterable[float]) -> float:
        #salva i parametri per traccia, ma NON li rileggo
        self.write_params(theta)
        E3, log = self.run_once(theta)   #PASSO THETA QUI
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
    th = tuple(a.theta) if a.theta else (-573.0, 1420.0, 1.55, 3.11)
    try:
        print(f"{tr.run_one(th):.6f}")
    except Exception as e:
        print(str(e), file=sys.stderr)
        sys.exit(2)

