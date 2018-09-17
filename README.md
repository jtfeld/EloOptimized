
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EloOptimized

[![Travis-CI Build
Status](https://travis-ci.org/jtfeld/EloOptimized.svg?branch=master)](https://travis-ci.org/jtfeld/EloOptimized)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/EloOptimized)](https://cran.r-project.org/package=EloOptimized)

[Package website](https://jtfeld.github.io/EloOptimized/)

EloOptimized provides tools to implement the maximum likelihood methods
for deriving Elo scores as published in [Foerster, Franz et al. (2016).
Chimpanzee females queue but males compete for social
status](https://www.nature.com/articles/srep35404). In addition, it
provides functionality to efficiently generate traditional Elo scores
using a simplified procedure that doesn’t require the use of cumbersome
presence matrices. Finally, it quickly generates a number of additional
Elo-based indices (ordinal, normalized, cardinal, and categorical ranks
and rank scores) of potential use to researchers, as outlined in the
linked manuscript.

## Installation

``` r
# Current version on Github:
# install.packages("devtools")
devtools::install_github("jtfeld/EloOptimized")

# CRAN-approved version on CRAN:
install.packages("EloOptimized")
```

## Example

There are two functions of interest. Use eloratingopt() to calculate Elo
scores using optimized Elo parameter values, or eloratingfixed() to
calculate Elo scores using user-defined parameter
values.

``` r
# to generate Elo scores using fixed initial Elo scores (1000) and a ML-fitted value for the K parameter:
nbaelo = eloratingopt(agon_data = nba, fit_init_elo = FALSE)

# to generate Elo scores using fixed default initial Elo scores and default K:
nbaelo = eloratingfixed(agon_data = nbadata, k = 100, init_elo = 1000)
```

To recreate the results from the 2016 manuscript, use the following
code:

``` r
# Males, model type 1:
melo1 = eloratingopt(agon_data = chimpagg_m, pres_data = chimppres_m, fit_init_elo = F)

# Males, model type 3:
melo3 = eloratingopt(agon_data = chimpagg_m[101:nrow(chimpagg_m),], 
                     pres_data = chimppres_m, fit_init_elo = T)

# Females, model type 1: 
felo1 = eloratingopt(agon_data = chimpagg_f, pres_data = chimppres_f, fit_init_elo = F)

# Females, model type 3:
felo3 = eloratingopt(agon_data = chimpagg_f[101:nrow(chimpagg_f),], 
                     pres_data = chimppres_f, fit_init_elo = T)
```
