import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")     #backend non-interattivo per salvare figure
from matplotlib import pyplot as plt

csv_file = sys.argv[1] if len(sys.argv) > 1 else "sweep_v01.csv"
df = pd.read_csv(csv_file)

# Plot 1: E3 vs parametro
plt.figure(figsize=(6,4))
plt.plot(df["x"], df["E3"], marker="o", color="tab:blue")
plt.title("3-body binding energy vs v01")
plt.xlabel("v01 [MeV]")
plt.ylabel("E3 [MeV]")
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot 2: Posterior normalizzato
plt.figure(figsize=(6,4))
plt.plot(df["x"], df["post_norm"], marker="o", color="tab:red")
plt.title("Posterior normalized vs v01")
plt.xlabel("v01 [MeV]")
plt.ylabel("Posterior weight (normalized)")
plt.grid(True)
plt.tight_layout()
plt.show()

import os, time
os.makedirs("figs", exist_ok=True)
ts = time.strftime("%Y%m%d-%H%M%S")

plt.tight_layout()
plt.savefig(f"figs/E3_vs_param_{ts}.png", dpi=200)

#per la seconda figura
plt.tight_layout()
plt.savefig(f"figs/posterior_vs_param_{ts}.png", dpi=200)
