# grapsupOU-simulation-estimation-application

R code to simulate and estimate Lévy-driven graph supOU processes, accompanying the article:

> Mehta, S. and Veraart, A. E. D. (2026). *Statistical inference for Lévy-driven graph supOU processes: From short- to long-memory in high-dimensional time series.* arXiv:2502.08838 [stat.ME]. https://arxiv.org/abs/2502.08838

---

## Repository structure

```
.
├── GraphSupOU-Functions-Final.R   # Core function library (source this first)
├── Portugal-Sim.R                 # Monte Carlo data generation (Section 4)
├── Portugal-Sim-Analysis.R        # Simulation study analysis (Section 4)
├── EmpiricalStudy-Final.R         # Empirical study (Section 5)
├── adj_mat.csv                    # 24-node Portuguese subnetwork adjacency matrix
├── stationary_data.csv            # Deseasonalised wind capacity factors
└── Previous versions/             # Original code from v1.0.0 release
```

---

## File descriptions

### `GraphSupOU-Functions-Final.R`
The core function library. Must be sourced before running any other script. Contains:
- **Simulation** — accelerated exact simulation of graph supOU processes via eigendecomposition of the drift matrix (Section 4.1, equation (eq:fastsim))
- **Graph normalisation** — column normalisation of the adjacency matrix (Section 2.3)
- **Estimation** — empirical autocovariance matrix and max real eigenvalue sequence (Section 3.1, equation (13))
- **Objective functions** — loss functions and fitted eigenvalue curves for the Gamma supOU, graph OU, and sum-of-two-exponentials models (equations (11), (12))
- **Unconstrained reparametrisation** — bijective transforms removing the constraints α > 1 and c ∈ (−1, 1) for gradient-based optimisation
- **Auxiliary matrix functions** — drift matrix constructor K(c) (equation (5))
- **Error measures** — RMSE and median absolute error

### `Portugal-Sim.R`
Generates 1000 independent Monte Carlo realisations of the graph supOU process on the 24-node Portuguese subnetwork, using the long-memory Gamma(1.5, 1) specification with parameters motivated by the empirical study (α = 1.5, c = −0.8, ϱ = 5, Δ = 1, N = 1000). Runs in parallel using all available cores minus one. Output is saved to `simdata/POR_sim_data_i.txt` for i = 1, …, 1000.

### `Portugal-Sim-Analysis.R`
Reads the simulated datasets produced by `Portugal-Sim.R` and reproduces the figures of Section 4:
- Boxplots of α̂ and ĉ across N* (Figures 2a, 2b)
- RMSE and MAE vs N* (Figure 2c)
- Violin plots at N* = 40 (Figure 2d)
- Median estimated Lévy basis mean and variance (Figures 2e, 2f)
- Scatter/density plots of (ĉ, α̂) at selected lags (Figure 3)

### `EmpiricalStudy-Final.R`
Reproduces the empirical study of Section 5, applied to wind capacity factors at 24 nodes in Portugal from the RE-Europe dataset (Jensen & Pinson, 2017). The script covers data loading, two-pass STL deseasonalisation, parameter estimation, and all figures in Section 5 of the paper (Figures 4–7).

---

## Data

### Included
- **`adj_mat.csv`** — 24×24 adjacency matrix for the Portuguese subnetwork shown in Figure 4a of the paper.
- **`stationary_data.csv`** — deseasonalised and detrended wind capacity factor time series (N = 26,304 hourly observations, d = 24 nodes), ready for use in `EmpiricalStudy-Final.R`.

### Required download
The raw wind capacity factor data (`wind_signal_COSMO.csv`) must be downloaded separately from the RE-Europe dataset:

> Jensen, T. and Pinson, P. (2017). RE-Europe, a large-scale dataset for modeling a highly renewable European electricity system. *Scientific Data* 4, 170175. https://doi.org/10.1038/sdata.2017.175

Place `wind_signal_COSMO.csv` in the working directory before running `EmpiricalStudy-Final.R`.

---

## How to reproduce the results

### Simulation study (Section 4)
```r
source("GraphSupOU-Functions-Final.R")  # load functions
source("Portugal-Sim.R")                # generate 1000 MC datasets (~minutes)
source("Portugal-Sim-Analysis.R")       # produce all figures
```

### Empirical study (Section 5)
```r
source("GraphSupOU-Functions-Final.R")  # load functions
source("EmpiricalStudy-Final.R")        # produce all figures
```
Note: `stationary_data.csv` is provided so the deseasonalisation step (which requires the raw `wind_signal_COSMO.csv`) can be skipped if desired. The estimation and plotting sections read from `stationary_data.csv` directly.

---

## Dependencies

| Package | Role |
|---|---|
| `expm` | Matrix exponential and logarithm |
| `igraph` | Graph utilities |
| `MASS` | Matrix operations |
| `ggplot2` | Plotting |
| `reshape2` | Data reshaping |
| `tidyr` | Data reshaping |
| `dplyr` | Data manipulation |
| `patchwork` | Combining plots |
| `ggrastr` | Rasterising large plots |
| `forecast` | STL decomposition |
| `parallel`, `doParallel`, `foreach` | Parallelisation |

Minimum R version: 4.1.0.

---

## Previous versions

The folder `Previous versions/` contains the original code as released in v1.0.0 (February 2025), prior to the addition of documentation, comments, and the restructuring described above.

---

## License

MIT — see `LICENSE` for details.

---

## Citation

If you use this code, please cite:

```bibtex
@misc{mehta2025statisticalinferencelevydrivengraph,
  title     = {Statistical inference for Levy-driven graph supOU processes:
               From short- to long-memory in high-dimensional time series},
  author    = {Shreya Mehta and Almut E. D. Veraart},
  year      = {2025},
  eprint    = {2502.08838},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  url       = {https://arxiv.org/abs/2502.08838}
}
```
## Acknowledgments

Claude.ai (Sonnet 4.6) was used to structure the code and repositories, improve the documentation, make the code more efficient, improve the graphics and provide first drafts of the README files.