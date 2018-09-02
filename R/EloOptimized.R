#' EloOptimized: ML fitting of Elo Scores
#'
#' @description This package implements the maximum likelihood methods for deriving Elo scores as published in Foerster, Franz et al. (2016).
#' Chimpanzee females queue but males compete for social status. Scientific Reports 6, 35404, doi:10.1038/srep35404
#'
#'
#' @section Primary functions:
#' \itemize{
#'   \item{\code{\link{eloratingopt}}: main function}
#'   \item{\code{\link{eloratingfixed}}: traditional Elo scores function}
#'   \item{\code{\link{elo.model1}}: internal function for fitting model type 1}
#'   \item{\code{\link{elo.model3}}: internal function for fitting model type 3}
#'   \item{\code{\link{elo.m3_lik_vect}}: vectorized internal function 
#'     for fitting mod type 3}
#' }
#' 
#' @section Plans for future development:
#' \itemize{
#'   \item Make package more modular, with a more flexible wrapper function.
#'   \item Option to specify K during burn-in period when fitting only K
#'   \item Add additional example data
#'   \item Create vignette, other package doohickies
#'   \item Add additional user control of the optimization procedure, allowing 
#'     for specification of the burn in period, optimization algorithm, and 
#'     initial values for optimization.
#'   \item Add functionality to plot Elo trajectories from within package.
#' }
#'
#' @docType package
#' @name EloOptimized
NULL
