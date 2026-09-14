# Generated results

Notebook-generated tables and simulation pickles are written in this
directory. Figures are written under `results/figures/`.

Figures are saved automatically as publication PDF and 300-dpi PNG
files. Set `SAVE_TEX = True` in a notebook's first code cell to add a
PGF-backed `.tex` version. The curated repository retains only the PNG
versions of the paper figures.

A 20,000-draw run of the continuous notebook writes:

```text
continuous_simulations.pkl
continuous_simulations_summary.csv
```

Preview or custom runs include their draw count in the filename. All
generated files are ignored by Git by default. The three placebo-analysis
figures used in the paper are retained as selected release artifacts.
