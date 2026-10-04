# Generated results

Notebook-generated tables and simulation pickles are written in this
directory. Figures are written under `results/figures/`.

Figures are saved automatically as publication PDF and 300-dpi PNG
files. Set `SAVE_TEX = True` in a notebook's first code cell to add a
PGF-backed `.tex` version. The curated repository retains only the PNG
versions of the paper figures.

The continuous notebook defaults to 100 repetitions per setting with
\(m=2000\), \(n=5\), \(\pi_0=0.99\), and seed 2, and writes:

```text
continuous_simulations_100_draws.pkl
continuous_simulations_100_draws_summary.csv
```

Set `draws = 10_000` in the first code cell for the full run; the filenames
change to `continuous_simulations_10000_draws.pkl` and
`continuous_simulations_10000_draws_summary.csv`. The pickle retains the
`results[tau][pi0]` layout, with SL results in `metrics`, BH results in
`bh_metrics`, and the run settings in `configuration` for each setting.

The two `t_test_*.png` figures show SL/BH rows and bFDR/FDR/power columns,
one figure per heterogeneity level. Rerunning replaces these figures.
The included figures use the saved 10,000-repetition run, not the preview.
Generated files are ignored by Git except the selected continuous results
and PNG figures, including the three placebo-analysis paper figures.
