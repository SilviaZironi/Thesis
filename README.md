# TRIMOD Bayesian Analysis
This project integrates the Fortran code **TRIMOD** with a Python-based Bayesian analysis framework,
allowing:
- automated execution of TRIMOD from Python;
- Bayesian inference on potential parameters;
- uncertainty quantification and stability studies.


---

## 📦 Project Structure

```
├── fortran/               # Original Fortran sources (main.f90, modules, makefile)
│   ├── m_trimod_sub.f90
│   ├── m_interpolation.f90
│   ├── m_nr_sub.f90
│   ├── m_nrutil.f90
│   ├── m_constants.f90
│   ├── m_nrtype.f90
│   ├── main.f90
│   └── makefile
│
├── runner.py              # Python wrapper: runs TRIMOD executable and parses output
├── bayes_trimod.py        # Bayesian layer: prior, likelihood, posterior, and 1D sweep
├── sweep.py               # Plotting utility for the sweep results (E3 & posterior)
├── TRIMOD_bayes_notebook.ipynb    # Jupyter notebook for interactive runs and plots
└── README.md              # This file
```

---

## ⚙️ Requirements

- **Python ≥ 3.8**
- **NumPy**, **Pandas**, **Matplotlib**

Install dependencies (once):
```bash
pip install numpy pandas matplotlib
```

Ensure your Fortran code is compiled:
```bash
cd fortran #path 
make
```

This will produce an executable named `run` (default name expected by the Python code).

---

## 🚀 1. Running TRIMOD from Python

### File: `runner.py`

`runner.py` provides a class `TrimodRunner` that automatically:
- Writes the potential parameters to `WORK/params.in`
- Executes the Fortran binary inside the working directory
- Reads `output.dat` and extracts the **3-body binding energy** value
- Detects non-convergence (e.g. "Too many iterations", "itermax", etc.)

#### Example usage (Python)
```python
from runner import TrimodRunner

tr = TrimodRunner(exe="./run", workdir="./WORK")
E3 = tr.run_one((-50.0, -30.0, 1.5, 2.0))
print("Triton binding energy:", E3, "MeV")
```

#### Example usage (command line)
```bash
python runner.py --exe ./run --workdir ./WORK --theta -50 -30 1.5 2.0
```

---

## 🧠 2. Bayesian Inference Layer

### File: `bayes_trimod.py`

This module defines the **Bayesian model** built on top of TRIMOD.

| Function | Description |
|-----------|--------------|
| `logprior(theta)` | Uniform prior over predefined parameter ranges |
| `loglike(theta, tr)` | Gaussian likelihood on \(E_3\) with experimental + model uncertainty |
| `logposterior(theta, tr)` | Returns \(\log P(\theta \mid \text{data})\) |
| `sweep1d(param, pmin, pmax, n, theta_fixed, tr, csv_out)` | Evaluates posterior along one parameter keeping the others fixed |

#### Example: single posterior evaluation
```bash
python bayes_trimod.py --exe ./run --workdir ./WORK --theta -50 -30 1.5 2.0
```

#### Example: 1D parameter sweep
```bash
python bayes_trimod.py --exe ./run --workdir ./WORK \
    --sweep1d v01 -120 -20 12 \
    --theta-fixed -50 -30 1.5 2.0 \
    --csv sweep_v01.csv
```

Output → `sweep_v01.csv` containing columns:
```
param, x, E3, logpost, post_norm
```

---

## 📊 3. Visualization of Results

### File: `sweep.py`

`sweep.py` reads the CSV file produced by `bayes_trimod.py` and plots:
- **E3 vs parameter**
- **Normalized posterior vs parameter**

Default input: `sweep_v01.csv`.  
Plots are saved automatically in `./figs/` with timestamps.

#### Example
```bash
python sweep.py
```

Produces:
```
figs/E3_vs_param_YYYYMMDD-HHMMSS.png
figs/posterior_vs_param_YYYYMMDD-HHMMSS.png
```

If running in Jupyter, you can visualize inline instead.

---

## 🧩 4. Typical Workflow

1. **Compile** the Fortran code (`make` inside `fortran/`)
2. **Run TRIMOD** via Python for a single parameter set:
   ```bash
   python runner.py --theta -50 -30 1.5 2.0
   ```
3. **Perform a 1D sweep** over one parameter:
   ```bash
   python bayes_trimod.py --sweep1d v01 -120 -20 12
   ```
4. **Plot** the results:
   ```bash
   python sweep.py
   ```
5. (Optional) Combine in Jupyter for analysis, plotting, or parameter fitting.

---

## ⚗️ 5. Notes on Integration

- All runs are executed inside `WORK/`, which can be safely cleaned between simulations.

- Non-convergent runs are handled gracefully:  
  they return `logposterior = -inf` and are automatically down-weighted.

---


## 🙌 Acknowledgements

This project was developed as part of the Master’s Thesis in **Nuclear and Subnuclear Physics**,  
focused on the Bayesian analysis of the triton ground state using TRIMOD.

Supervisor: **Prof. Paolo Finelli**  
Student: **Silvia Zironi**

