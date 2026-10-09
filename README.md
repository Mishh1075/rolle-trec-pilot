# Do t-Rec intervals keep their coverage when grid load grows? A pilot on the Rolle hierarchy

A small, reproducible pilot comparing two ways of reconciling probabilistic forecasts on a real distribution-grid hierarchy:

- **MinT-style Gaussian reconciliation** (`reconc_gaussian`)
- **t-Rec** (`reconc_t`; Carrara, Azzimonti, Corani and Zambon, *Modeling uncertainty in the covariance matrix for probabilistic forecast reconciliation*, arXiv 2506.19554)

Both come from the R package [bayesRecon](https://github.com/IDSIA/bayesRecon). The question is what happens to their 95% prediction intervals when load grows in part of the grid.

Author: Misha Parekh · October 2026

## Data

The Rolle distribution-grid benchmark by Nespoli, Medici, Lopatichki and Sossan ([arXiv 1910.03976](https://arxiv.org/abs/1910.03976); data on [Zenodo, 10.5281/zenodo.4549296](https://doi.org/10.5281/zenodo.4549296), CC BY 4.0). The data are **not** included here; download `reduced_dataset.pk` from Zenodo.

- 24 meters, renamed m01–m24 in alphabetical order (`meter_names.csv` gives the mapping), and 7 aggregates.
- Hierarchy, checked on the data: all = S1 + S2; S1 = S11 + S12; S2 = S21 + S22; S11, S12, S21 and S22 each hold six meters (m01–06, m07–12, m13–18, m19–24).
- 13 Jan 2018 – 19 Jan 2019, 10-minute data with no gaps, averaged to 8,928 hourly rows.

## Design (all settings fixed before any result was seen)

| Item | Choice |
|---|---|
| Forecast | One hour ahead, 100 rolling origins, 73 hours apart (so the hour of day rotates) |
| Base models | Seasonal naive; seasonal AR fitted per series, y(t) = a + b·y(t−1) + c·y(t−24) |
| Methods | No reconciliation (Gaussian from each series' residuals), MinT (`reconc_gaussian`), t-Rec (`reconc_t`, `freq = 24`) |
| Load growth | +0, +5, +10, +20% in meters m01–m05 (all in feeder S11), starting 12 hours before the forecast hour. The upper series are recomputed so the data stay coherent. |
| Training window | 48, 96, 168 and 672 hours (2, 4, 7 and 28 days) |
| Scores | Coverage at 80/90/95%, mean width, interval score; for the affected series (all, S1, S11) and per level |

Why a seasonal AR model as well: seasonal-naive forecasts are already coherent, so t-Rec's disagreement term is always zero. Per-series models make the base forecasts disagree.

## Results (affected series, 95% intervals, seasonal-AR base model)

| Training | Growth | Coverage: none / MinT / t-Rec | Interval score, kW: none / MinT / t-Rec |
|---|---|---|---|
| 28 days | 0% | 0.963 / 0.970 / 0.953 | 99.7 / 80.5 / **80.0** |
| 28 days | +20% | 0.927 / 0.863 / 0.803 | 113.7 / **113.6** / 138.6 |
| 2 days | 0% | 0.873 / 0.840 / 0.870 | 119.8 / 125.1 / **111.0** |
| 2 days | +20% | 0.877 / 0.870 / 0.900 | 124.9 / 125.9 / **115.0** |

- **With 28 days of training**, reconciliation cuts the interval score by about 20% when there is no growth. Under +20% growth, t-Rec loses the most coverage: in 18 of 300 forecasts MinT covered the true value and t-Rec did not, and never the reverse.
- **With 2 days of training**, t-Rec has the best interval score at every growth level, at 80, 90 and 95%, and at every level of the hierarchy. With no growth and with +20% growth alike, t-Rec covered and MinT did not in 9 of 300 forecasts, and never the reverse.
- **The likely reason:** t-Rec's degrees of freedom grow with the training length (about 70 at 2 days, 682 at 28 days). With long training it behaves like a Gaussian, and its intervals barely react to growth (71 to 72 kW). For a growing grid, the question becomes the trade-off between short windows that adapt and long windows that estimate the covariance well.

The full tables are in `results_core.csv`, `results_training_length.csv`, `results_upgrade_affected.csv` and `results_upgrade_by_level.csv`. The plot is `coverage_vs_growth.png`.

## What this does not show

- **One dataset, one year, 100 origins.** Differences of 1–3 coverage points are close to noise. No confidence intervals or paired tests were computed beyond the hit counts above.
- **The growth is injected, not observed.** It is a step change in five meters of one feeder, starting 12 hours before the forecast. Real growth (EV chargers, heat pumps) is gradual and spread out.
- **Simple base models only**, one hour ahead. ETS, LightGBM and longer horizons were not tested, and only hourly data were used.
- **Short windows undercover for every method.** At 2 days, all 95% intervals cover only 0.84–0.90. t-Rec is the least bad there, not well calibrated.
- **The aggregates are exact sums of the meters.** It is not yet confirmed whether they were metered or constructed, which may matter for how realistic the errors are.
- **No CRPS or energy score yet**, and t-Rec was used with the package's default prior settings.

## How to run

Requirements: Python 3 with `pandas` and `pytz` (`scipy` for `find_hierarchy.py`), and R 4.6.1 with `bayesRecon` 1.0.2. `session_info.txt` lists the exact versions used.

1. Put `reduced_dataset.pk` (from Zenodo) in this folder.
2. In a terminal, run `python make_hourly.py`. It writes `rolle_hourly.csv` and `meter_names.csv`.
3. Open R or RStudio with this folder as the working directory and run `source("run_all.R")`. It builds and checks A, runs the tests, runs the three pilot scripts and saves every results file. It takes about 10 minutes.

To run only the tests: `source("tests/test_pilot.R")`.

## Files

| File | What it does |
|---|---|
| `inspect_rolle.py` | Prints the dataset's shape, time range, resolution and missing values |
| `find_hierarchy.py` | Recovers which meters sum to each aggregate (non-negative least squares) |
| `make_hourly.py` | Writes the hourly CSV, upper series first, then m01–m24 |
| `build_A.R` | Builds the 7 × 24 matrix A and checks A × meters = aggregates |
| `pilot_core.R` | Run 1: 28-day window, both base models, 95% coverage and width |
| `pilot_training_length.R` | Run 2: training windows of 2, 4, 7 and 28 days |
| `pilot_upgrade.R` | Run 3: 80/90/95%, interval score, per level, and the plot |
| `tests/test_pilot.R` | Five sanity checks: shape of A, complete data, coherence, injected growth, interval score |
| `run_all.R` | Runs everything in order |

## Acknowledgements

Data: Nespoli, Medici et al., Rolle benchmark (CC BY 4.0). Methods: the bayesRecon package (IDSIA) and the t-Rec paper.
