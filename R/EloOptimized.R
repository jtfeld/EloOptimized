#' EloOptimized: ML fitting of Elo Scores
#'
#' @description This package implements the maximum likelihood methods for deriving Elo scores as published in Foerster, Franz et al. (2016).
#' Chimpanzee females queue but males compete for social status. Scientific Reports 6, 35404, doi:10.1038/srep35404
#'
#'
#' @section Primary functions:
#' \code{\link{eloratingopt}},
#' \code{\link{elo.model1}},
#' \code{\link{elo.model3}},
#' \code{\link{elo.m3_lik_vect}},
#' \code{\link{eloratingoptR}},
#' \code{\link{eloratingopt_simple}}
#' 
#' @section Ideas for future development:
#' \itemize{
#'   \item option to include additional rank scores or not (ie leave out ordinal, scaled, etc...)
#'   \item rather than returning a data frame, we could return a fit object with $daily_elo, 
#'     $optimized_k, $prediction_accuracy, $log_likelihood, $aic, $aicc, etc...
#'   \item maybe name the models by model type, rather than M or F?  That way we're not presupposing
#'     that males all have the same entry value, etc...
#'   \item option to go back and recalculate during burn-in period with optimized k for Model 1?
#'   \item find data we can use in vignette
#'   \item create vignette, other package doohickies
#' }
#'
#' @docType package
#' @name EloOptimized
NULL
