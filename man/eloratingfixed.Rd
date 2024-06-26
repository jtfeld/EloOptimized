% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eloratingfixed.R
\name{eloratingfixed}
\alias{eloratingfixed}
\title{Create daily elo ranks and multiple derivatives with user-defined parameter values}
\usage{
eloratingfixed(agon_data, pres_data, k = 100, init_elo = 1000, outputfile = NULL, 
  returnR = TRUE, p_function = "sigmoid")
}
\arguments{
\item{agon_data}{Input data frame with dominance interactions, should only contain Date, 
Winner, Loser.  Date should be formatted as MONTH/DAY/YEAR, or already as Date class.}

\item{pres_data}{Input data frame with columns "id", "start_date" and "end_date".  Date
columns should be formatted as MONTH/DAY/YEAR, or already as Date class.  If all IDs 
are present the whole time, you ignore this and a pres_data table will be automatically
generated.}

\item{k}{Specified value of the k parameter, default is 100}

\item{init_elo}{The starting Elo value for all individuals, default is 1000}

\item{outputfile}{Name of csv file to save ranks to.  Default is NULL, in which case 
the function will only return a table in R.  If you supply an output file name
the function will save the results as a csv file in your working directory.}

\item{returnR}{whether to return an R object from the function call.  Default is TRUE}

\item{p_function}{function defining probability of winning.  Default "sigmoid" is 
equation (1) from Foerster, Franz et al 2016.  Use "pnorm" to use the 
\code{\link[stats:Normal]{pnorm}}-based method implemented in the EloRating package.}
}
\value{
Returns a list with six elements: 
\describe{
 \item{\strong{elo}}{ Data frame with all IDs and dates they were present, with the following columns:}
   \describe{
     \item{Date}{: Dates of study period}
     \item{Individual}{: the names of each ranked individual, for each date they were present}
     \item{Elo}{: fitted Elo scores for each individual on each day}
     \item{EloOrdinal}{: Daily ordinal rank based on Elo scores}
     \item{EloScaled}{: Daily Elo scores rescaled between 0 and 1 according to 
       \deqn{([individual Elo] - min([daily Elo scores])/(max([daily Elo scores]) - min([daily Elo scores]))}}
     \item{ExpNumBeaten}{: expected number of individuals in the group beaten, which is the sum of 
       winning probabilities based on relative Elo scores of an individual and all others, following 
       equation (4) in Foerster, Franz et al. 2016}
     \item{EloCardinal}{: ExpNumBeaten values rescaled as a percentage of the total number of ranked 
       individuals present in the group on the day of ranking. We encourage the use of this measure.}
     \item{JenksEloCardinal}{: Categorical rank (high, mid, or low) using the Jenks natural breaks 
       classification method implemented in the R package BAMMtools. 
       See \code{\link[BAMMtools]{getJenksBreaks}}}
     }
 \item{\strong{k}}{ User-defined value of the k parameter}
 \item{\strong{init_elo}}{ User-defined initial Elo score when individuals enter the hierarchy}
 \item{\strong{pred_accuracy}}{ Proportion of correctly predicted interactions}
 \item{\strong{logL}}{ The overall log-likelihood of the observed data given the user-supplied parameter 
   values based on winning probabilities (as calculated in equation (1) of Foerster, Franz et al 2016) 
   for all interactions}
 }
}
\description{
Conducts traditional elo rating analyses using specified K value
  and outputs raw, normalized, cardinal, and  categorical ranks as a list object in 
  R or in an output file.  For optimized Elo parameters, use \code{\link{eloratingopt}}.
}
\details{
This function accepts a data frame of date-stamped dominance interactions and 
  (optionally) a data frame of start and end dates for each individual to be ranked, 
  and outputs daily Elo scores with parameters specified by the user.  The default function 
  used to determine probability of winning is equation (1) from Foerster, Franz et al. 2016, 
  but for ease of comparison with the EloRating package, we also added the option to use
  the \code{\link[stats:Normal]{pnorm}}-based method implemented in the EloRating package, and future 
  development will add the option to use the original function from Elo 1978 (as implemented in 
  the elo package).  This function does not require large presence matrices, and efficiently 
  calculates a series of additional indices (described below).  
  
  As opposed to the \code{\link{eloratingopt}} function, this procedure only requires that 
  included individuals have at least one win \emph{or} one loss. 
  
  A detailed description of the function output is given in the \strong{Value} section of 
  this help file:
}
\examples{

nbadata = EloOptimized::nba #nba wins and losses from the 1995-96 season
nbaelo = eloratingfixed(agon_data = nbadata)
# generates traditional Elo scores (with init_elo = 1000 & k = 100) and saves 
#   them as "nbaelo" 

}
