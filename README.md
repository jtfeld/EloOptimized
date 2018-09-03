
<!-- README.md is generated from README.Rmd. Please edit that file -->
EloOptimized
============

EloOptimized provides tools to implement the maximum likelihood methods for deriving Elo scores as published in [Foerster, Franz et al. (2016). Chimpanzee females queue but males compete for social status](https://www.nature.com/articles/srep35404). In addition, it provides functionality to efficiently generate traditional Elo scores using a simplified procedure that doesn't require the use of cumbersome presence matrices. Finally, it quickly generates a number of additional Elo-based indices (ordinal, normalized, cardinal, and categorical ranks and rank scores) of potential use to researchers, as outlined in the linked manuscript.

Installation
------------

You can install EloOptimized from github with:

``` r
# install.packages("devtools")
devtools::install_github("steffenfoerster/elorating/tree/jtf_devbranch")
```

Example
-------

There are two functions of interest. Use eloratingopt() to calculate Elo scores using optimized Elo parameter values, or eloratingfixed() to calculate Elo scores using user-defined parameter values.

``` r
# to generate Elo scores using fixed initial Elo scores (1000) and a ML-fitted value for the K parameter:
nbaelo = eloratingopt(agon_data = nba, fit_init_elo = FALSE)

# to generate Elo scores using fixed default initial Elo scores and default K:
nbaelo = eloratingfixed(agon_data = nbadata, k = 100, init_elo = 1000)
```
