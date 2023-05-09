# EloOptimized 0.3.1.9000 (unreleased)

* Added a step to check that agon_data are in chronological order, and throw a warning if not.

* Fixed checks for ubuntu and Mac

* Added new and improved informative errors to eloratingopt() that tell you which individuals are filtered out for not having >= 1 win and 1 loss, and which individuals have interactions outside of their presence window.

# EloOptimized 0.3.1

* Mostly a maintenance update

* Removed tcltk, rlist dependencies (holdovers from older code)

* Streamlined some things under the hood to make post-processing more efficient in eloratingopt & eloratingfixed

* EloOrdinal column in output now stored as an integer object for small space saving.

* Added a vignette to demonstrate usage of the package, see [here](https://jtfeld.github.io/EloOptimized/articles/nba-example.html).

* Added cleaning step to avoid error when trying to use tibbles imported with readr::read_csv()

* Added additional informative errors if input data are problematic

* Fixed small bug related to converting character dates to Date class

* Incorporated github actions for CI checks, pkgdown build

# EloOptimized 0.3.0

* Added Travis-ci integration in preparation for submission to CRAN

# EloOptimized 0.2.1

* Updated documentation and parameter names
* Removed AIC element from `eloratingfixed` output

# EloOptimized 0.2.0

* Added a `NEWS.md` file to track changes to the package.
