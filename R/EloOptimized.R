#' EloOptimized: ML fitting of Elo Scores
#'
#' @description This package implements the maximum likelihood methods for deriving Elo scores as published in Foerster, Franz et al. (2016).
#' Chimpanzee females queue but males compete for social status. Scientific Reports 6, 35404, doi:10.1038/srep35404
#'
#'
#' @section Primary functions:
#' \itemize{
#'   \item{\code{\link{eloratingopt}}: main function}
#'   \item{\code{\link{eloratingtrad}}: traditional Elo scores function}
#'   \item{\code{\link{elo.model1}}: internal function for fitting model type 1}
#'   \item{\code{\link{elo.model3}}: internal function for fitting model type 3}
#'   \item{\code{\link{elo.m3_lik_vect}}: vectorized internal function 
#'     for fitting mod type 3}
#' }
#' 
#' @section Ideas for future development:
#' \itemize{
#'   \item option to include additional rank scores or not (ie leave out ordinal, scaled, etc...)
#'   \item option to go back and recalculate during burn-in period with optimized k for Model 1?
#'   \item find data we can use in vignette
#'   \item create vignette, other package doohickies
#'   \item add fit_k and fit_init_elo options to eloratingopt function
#' }
#'
#' @docType package
#' @name EloOptimized
NULL
