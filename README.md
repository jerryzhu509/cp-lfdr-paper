# A Local False Discovery Rate Framework for Compound p-values

Reproducible simulations and applications accompanying the manuscript.
The repository contains the continuous Gaussian-means study, the GSE9101
gene-expression analysis, and the Romanian admissions
regression-discontinuity analysis.

## Repository contents

- `notebooks/01_continuous_simulations.ipynb`: Gaussian-means simulation.
- `notebooks/02_gene_expression.ipynb`: GSE9101 application.
- `notebooks/03_regression_discontinuity.ipynb`: Romanian admissions
  RD application.
- `notebooks/03_regression_discontinuity_placebo.ipynb`: midpoint-placebo
  falsification analysis for the RD application.
- `src/cp_lfdr/`: shared, tested implementations.
- `tests/`: unit and regression tests.
- `data/README.md`: data provenance and local setup.
- `results/`: generated figures, summaries, and simulation pickles.

## Running the analyses

After placing the data as described in `data/README.md`, open any of the
four notebooks in Jupyter and choose **Run All Cells**. Each notebook
finds the local `src/cp_lfdr` package automatically, can be run
independently, and writes its own outputs under `results/`.

The first code cell in each notebook contains all editable run settings:

| Setting | Notebook | Purpose |
| --- | --- | --- |
| `draws` | Continuous | Defaults to 100 repetitions per setting; change to `10_000` for the full run. |
| `SAVE_FIGURES` | All | Saves PDF and 300-dpi PNG figures when `True`. |
| `SAVE_TEX` | All | Also saves PGF-backed `.tex` figures when `True`. |
| `SAVE_RESULTS` | Continuous | Saves the simulation pickle and summary CSV. |
| `DATA_DIR` | Gene/RD | Can be edited if data are stored elsewhere. |
| `BANDWIDTHS` | RD placebo | Selects the midpoint-placebo bandwidth grid. |

No environment variables or terminal commands are needed to select a
run mode.

## Paper figures and simulation output

Figures are saved to `results/figures/`. With `SAVE_TEX = True`, the
notebooks use Matplotlib's PGF backend to write `.tex` files suitable
for inclusion in LaTeX.

The continuous preview (`draws = 100`) writes:

```text
results/continuous_simulations_100_draws.pkl
results/continuous_simulations_100_draws_summary.csv
results/figures/t_test_0.0.pdf
results/figures/t_test_0.3.pdf
```

The continuous simulation uses \(m=2000\), \(n=5\), \(\pi_0=0.99\),
and \(\alpha=0.1\). The two figures present the homogeneous
\((\tau=0)\) and heterogeneous \((\tau=0.3)\) scenarios separately, with
SL/BH rows and bFDR/FDR/power columns. Both procedures use the same data.
The curves are t-test, density-compound, weighted, and e-weighted; the
data-driven construction is commented out. Pickle and CSV filenames
include the repetition count, so a 10,000-draw run preserves the preview
results while replacing the figures. Each pickle records its run settings.
The included paper figures and saved notebook plots use 10,000 repetitions
per setting (seed 2); the notebook defaults to 100 for a quick preview.
Set `draws = 10_000` to reproduce the paper figures.

The gene-expression notebook writes the theoretical inflation curve,
rejection curves, p-value histogram, and Q-Q tail plot using the paper
styles. The main RD notebook writes the local RD fits, rejection curves,
theoretical inflation curve, p-value histogram, Q-Q plot, and exact-pool
comparison. The RD placebo notebook constructs midpoint placebos between
adjacent real cutoffs and writes the Town 20787 illustration, rejection
curves, and p-value histograms across bandwidths. These three selected
placebo figure sets are retained in the repository as paper artifacts.

## Python environment

If the notebook kernel already contains the dependencies listed in
`pyproject.toml`, no installation step is needed. A new environment can
be prepared with:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev,notebooks]"
pytest
```

## Reproducibility

The continuous and RD analyses use seed `2`. For each RD quasi-experiment,
they sample 69 distinct control assignments uniformly without replacement
after excluding the observed assignment. Adding the observation gives a
common block size of 70. Each placebo design restarts the same seed-2
`RandomState` stream.

## Data and licensing

The processed GSE9101 matrix is retained for reproducibility. The much
larger RD dataset is not committed and must be obtained from the original
replication package. See `data/README.md` for provenance and setup.
