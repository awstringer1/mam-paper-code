# Flexible Marginal Models for Dependent Data

This repository contains scripts to replicate the results in the paper *Flexible Marginal Models for Dependent Data* (**arXiv**). Please follow the below instructions to insatll software and run the analyses.

# Install Packages

## mam

The `mam` package was written to implement the methods developed in this paper. It has one core function, `mam::mam`, that is used within the scripts in this repository.

The `mam` package is in early development and is on [github](https://github.com/awstringer1/mam) but not `CRAN`. To install the `mam` package, please type

```
install.packages('devtools')
devtools::install_github('awstringer1/mam')
```

The `mam` package depends on compiled `C++` code in the form of `TMB` scripts. With apologies, this has not been widely tested, especially on Windows.

## CRAN

The following packages must be installed from `CRAN` in order to run all the scripts:
```
install.packages(c(
  "mgcv",
  "gamm4",
  "tidyverse",
  "parallel"
))
```

## Other Repositories

One package, `geostatsp`, is no longer on `CRAN` at the time of writing. It is required for the `loaloa` example of Section 5.2. It can be installed from `R-forge`:
```
install.packages("geostatsp", repos="http://R-Forge.R-project.org")
```

# Replicate the Analysis

There are three scripts in this repository, each replicating one analysis from the paper. You should not have to make any modifications to any of the scripts. Simply source them each in a clean `R` session, and the results will be saved to temporary storage (`temp.dir()`) on your machine. The location will be printed at the end of each script.

All the data for the examples is available within the `R` packages listed above, so no data needs to be downloaded externally.

The analyses are as follows:

- **Section 4** and **Supplement B**: `01-simulations.R`
- **Section 5.1**: `02-beavers.R`
- **Section 5.2**: `03-loaloa.R`

If any questions, do not hesitate to raise an issue on this repository. Thank you!
