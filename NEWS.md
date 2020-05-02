# EloOptimized 0.3.0.9000

* Removed tcltk, rlist dependencies (holdovers from older code)

* Streamlined some things under the hood to make post-processing more efficient in eloratingopt & eloratingfixed

* EloOrdinal column in output now stored as an integer object for small space saving.

* Added a vignette to demonstrate usage of the package.

* Added cleaning step to handle tibbles imported with readr::read_csv()

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