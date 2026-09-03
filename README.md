# uniform-over-dimension-location

This repository contains the code for reproducing the results in the paper 
"Uniform-over-dimension location tests for multivariate and high-dimensional data" 
by [anonymized] ([anonymized arxiv link]). The implementation of the proposed 
methods is provided in `code/methods.R` and can be used independently of the simulation 
framework.

## Repository Structure

```text
.
├── unif_over_dim.Rproj         # RStudio project
|
├── run_simulation.R            # Run the complete simulation study
├── make_simulation_tables.R    # Generate Tables "tables/Rtable_p={5,25,50,100}.tex" 
├── make_simulation_plots.R     # Generate Figures "figures/Rplot_p={5,25,50,100}.pdf"
├── real_data_analysis.R        # Perform real data analysis, generate "figures/Rplot_histogram.pdf"
|
├── code/                       # R functions
│   ├── methods.R               # Proposed and competing methods
│   ├── simulation.R            # Data generation and simulation functions
│   └── utils.R                 # Utility functions for result processing
|
├── results/                    # Simulation results
│   └── sim_results.rds
|
├── tables/                     # Tables generated from simulation results
|   ├── Rtable_p=5.tex          # Table 1
|   ├── Rtable_p=25.tex         # Table S1
|   ├── Rtable_p=50.tex         # Table S2
|   └── Rtable_p=100.tex        # Table 2
|
└── figures/                    # Figures generated from simulation results and real data analysis
    ├── Rplot_p=5.pdf           # Figure 1
    ├── Rplot_p=25.pdf          # Figure S1
    ├── Rplot_p=50.pdf          # Figure S2
    ├── Rplot_p=100.pdf         # Figure 2
    └── Rplot_histogram.pdf     # Figure S3
```

## Requirements

The code can be run in RStudio with the `unif_over_dim.Rproj` file or from an R
session with the repository root as the working directory.

The proposed methods can be used by sourcing `code/methods.R` and do not require
any additional R packages.

The simulation study and real data analysis use the following R packages:

`HDNRA`, `mvtnorm`, `highmean`, `DescTools`, `maotai`, `energy`, `future`, 
`future.apply`, `knitr`, `ggplot2`, `patchwork`, `colorspace`, `plsgenomics`.

To install the required packages (except `highmean`), run:

```r
install.packages(c(
  "HDNRA", "mvtnorm", "DescTools", "maotai", "energy",
  "future", "future.apply", "knitr", "ggplot2", "patchwork",
  "colorspace", "plsgenomics"
))
```

The `highmean` package is available in the [CRAN Archive](https://cran.r-project.org/src/contrib/Archive/highmean/).
It can be installed using `remotes`:

```r
install.packages("remotes")
remotes::install_url(
  "https://cran.r-project.org/src/contrib/Archive/highmean/highmean_3.0.tar.gz"
)
```

## Using the Proposed Method

Simply source `code/methods.R`:

```r
source("code/methods.R")
```

The proposed method can then be called directly using the corresponding function
in `code/methods.R`. 

To run the two sample test on the data matrices `dataX` and `dataY` using the 
spatial kernel and 10000 replications for cut-off estimation, run:

```r
result <- kcdg_test(dataX = dataX, dataY = dataY, h = h_spatial, nsim = 10000, estimators = c(1,0))
print(result)
```

To use the test based on the tapering estimator with parameter values 0.1, 0.25 and 0.4, 
run:

```r
result <- kcdg_test(dataX = dataX, dataY = dataY, h = h_spatial, nsim = 10000, estimators = c(0,1), vec_beta = c(0.1, 0.25, 0.4))
print(result)
```

In this case, an array of 3 p-values is returned. To use the difference kernel 
with both the plain and tapering estimators, run:

```r
result <- kcdg_test(dataX = dataX, dataY = dataY, h = h_diff, nsim = 10000, estimators = c(1,1), vec_beta = c(0.1, 0.25, 0.4))
print(result)
```
In this case, an array of 4 p-values is returned with the first element 
corresponding to the test based on plain estimator and the remaining elements 
corresponding to the tests based on the tapering estimators.

## Reproducing the Simulation Results

The complete simulation study can be reproduced from scratch by running:

```r
source("run_simulation.R")
```

The script uses parallel processing via the `future` package and, by default,
uses one fewer worker than the number of available CPU cores. The number of 
workers can be adjusted by changing `num_workers` in the script 
`run_simulation.R` based on the available computational resources.

The simulation results are automatically saved to `results/sim_results.rds`.

### Generating Tables and Figures

The simulation results used in the paper are already provided in
`results/sim_results.rds`. Therefore, the tables and figures can be reproduced
without rerunning the simulation study.

Run:

```r
source("make_simulation_tables.R")
source("make_simulation_plots.R")
```

`make_simulation_tables.R` generates Tables 1 and 2 from the main paper and 
Tables S1 and S2 from the Supplementary Material, and saves the corresponding 
`.tex` files in `tables/`.

`make_simulation_plots.R` generates Figures 1 and 2 from the main paper and
Figures S1 and S2 from the Supplementary Material, and saves the corresponding 
`.pdf` files in `figures/`.

### Obtaining Size and Power Results in Tabular Form

Figures 1, 2, S1, and S2 summarize the size and power results for the hypothesis 
testing methods considered in the simulation study. If the numerical size and power
results are desired in tabular form, the following script can be used.

```r
source("code/utils.R")
results <- readRDS("results/sim_results.rds")

p <- 100          # dimension
model <- 1        # model

if (p == 5) {
  mu_values <- c(0, 0.2, 0.4, 0.6, 0.8, 1, 1.2)
} else {
  mu_values <- 0:6
}

output <- empirical_power(p, model, mu_values, results)
print(output)
```

The `p` parameter specifies the dimension and can be set to `5`, `25`, `50` or
`100`. The `model` parameter specifies the simulation model and can be set to
values from `1` to `9`, corresponding to Models 1.i, 2.i, 3.i, 1.ii, 2.ii, 3.ii,
1.iii, 2.iii, and 3.iii, respectively.

For example, the code above produces the size and power results for dimension
`p = 100` and Model 1.i. To obtain results for a different dimension or model,
simply change the values of `p` and `model` accordingly.


## Reproducing the Real Data Analysis

The real data analysis can be reproduced independently by running:

```r
source("real_data_analysis.R")
```

The script analyzes the human colon tissue dataset from the R package `plsgenomics`. 
It outputs the results corresponding to Table 3 of the main paper in the
R console and generates Figure S3 of the Supplementary Material at `figures/Rplot_histogram.pdf`.