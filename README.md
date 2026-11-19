# Least-squares factor analysis at a Heywood boundary: R code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24278893.svg)](https://doi.org/10.5281/zenodo.24278893)

R code for all numerical results of the article

> Ouhssaine S, El Hadri Z. Least-squares factor analysis at a Heywood boundary: a penalty family and a bounded diagnostic in marginal-variance coordinates. 

Archived release (v1.0.0, the version cited in the article): [doi:10.5281/zenodo.24278893](https://doi.org/10.5281/zenodo.24278893).

Soufiane Ouhssaine and Zouhair El Hadri, Laboratory of Mathematics, Statistics and Applications (LMSA), Department of Mathematics, Faculty of Sciences, Mohammed V University in Rabat, Morocco.

Corresponding author: Soufiane Ouhssaine, soufiane.ouhssaine@um5r.ac.ma

## Contents

| Script | What it computes |
|---|---|
| `Estimators.R` | Populations of the simulation and the four criteria: E1 (unrestricted MINRES), E2 (MINRES with $C_{ii} \leq S_{ii}$), E3 (ULS with nonnegative uniquenesses, two algorithms) and E4 (Gaussian ML with uniquenesses bounded below by 0.005). Also the multiplier of the cap at a capped fit and the Wald statistic of the smallest uniqueness of unrestricted ULS (estimate over standard error, from lavaan). The other scripts read it with `source("Estimators.R")`. |
| `Penalty path.R` | Table 1 (penalty path of the 3 x 3 boundary example, by a refined grid search) and the numerical checks of Sect. 6.7 (322 matrices) and of the multiplier of the cap. |
| `Validation.R` | Checks run before the simulation: the populations, the closed-form values of the 3 x 3 example, agreement of the two E3 algorithms, E1 against `psych::fa`, interior E3 fits equal to E1, and the profile identity. Every row of the printed table should have `pass` equal to `cases`. The last block refits every population matrix from 80 starting values and prints whether they give one common covariance (Sect. 6.2). |
| `Simulation.R` | The simulation of Sects. 6.2 to 6.5: 46 populations, n = 50, 100, 200, 500 and 1000, 1000 samples per cell, E1 to E4 on every sample, and the multiplier of the cap and the Wald statistic on the samples with a Heywood case. Saves one file per design cell in `Simulation_cells/` and, when all 230 are there, `Simulation_results.rds`. |
| `Simulation summary.R` | Reads `Simulation_results.rds` and prints Tables 2, 3 and 4 and the other numbers of Sects. 6.3 to 6.5. Draws `Fig1.pdf`, `Fig2.pdf` and `Fig3.pdf`. The intervals of Table 4 use $1000$ bootstrap resamples (`B_auc` at the top of the script). |
| `Bootstrap.R` | The parametric bootstrap test of the residual (Sect. 6.6), $B = 199$. Saves one file per cell in `Bootstrap_cells/`, prints the rejection rates and draws `Fig4.pdf`. |
| `Real data application.R` | Sect. 7: the two screens of data sets distributed with R, psych, lavaan and MASS, Tables 5 and 6 with the bootstrap check (B = 499), the spread of the fits over 50 starting values, and the smallest uniqueness returned by `psych::fa`. |

## Requirements

R 4.6.1 with the packages psych 2.6.9, lavaan 0.7.3 and MASS 7.3.65:

```r
install.packages(c("psych", "lavaan", "MASS"))
```

## How to run

Clone the repository, make its folder the working directory, and run each script from top to bottom, for example with `source("Simulation.R")`. The scripts are independent, except that `Simulation summary.R` needs the file `Simulation_results.rds` written by `Simulation.R`. A suggested order is

1. `Validation.R`
2. `Penalty path.R`
3. `Simulation.R`, then `Simulation summary.R`
4. `Bootstrap.R`
5. `Real data application.R`

`Validation.R`, `Penalty path.R` and `Real data application.R` set the seed with `set.seed(123)` once. `Simulation.R` and `Bootstrap.R` set the seed 123 + c at the start of design cell c and save each cell in its own file, so their cells can also be run in several R sessions side by side, with the same results. For example, six sessions started from a terminal with

```
Rscript Simulation.R 1 230 6
Rscript Simulation.R 2 230 6
...
Rscript Simulation.R 6 230 6
```

each run every sixth cell (arguments: first cell, last cell, step). Finished cells are skipped when a script is started again, so an interrupted run can be resumed. To start afresh, or after changing `nbsim`, delete the folder `Simulation_cells` (`Bootstrap_cells`). The number of samples is set at the top of each script (`nbsim` in `Simulation.R`, `nbsim` and `B_boot` in `Bootstrap.R`) and can be lowered for a quick trial run.

`Simulation summary.R` and `Bootstrap.R` write their figures with `cairo_pdf()`, which embeds the fonts. If R was built without cairo, replace `cairo_pdf` by `pdf`.

## Relation to the results in the article

All numbers, tables and figures in the article were computed with these scripts and the package versions listed above: `Simulation.R` and `Bootstrap.R` in six R sessions as in the example, the other scripts in one session. A rerun with the same versions should give the same results. Another version of R, of a package or of the linear algebra library can change single fits, because each fit is the best of several local searches, and with them the last digits of the simulation results.

The simulated data are not stored here; the scripts generate them. The published data come with R and with the packages psych, lavaan and MASS.

## License

MIT, see `LICENSE`. To cite the code, use the article or the Zenodo record; `CITATION.cff` has the details.
