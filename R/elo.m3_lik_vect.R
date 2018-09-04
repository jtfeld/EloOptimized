#' @title optimize k parameter and entry Elo scores, vectorized
#' @description Function to optimize k parameter and entry Elo scores
#' @usage elo.m3_lik_vect(par, IA_data, all_ids)
#' @param par list of parameters, with par[1] being log(k), and par[2:length(par)] 
#'   being the initial elo scores of individuals
#' @param IA_data list of interaction data, with columns "Date", "Winner", and "Loser"
#'   (in that order)
#' @param all_ids list of all ids to rank
#' @examples
#' 
#' # for internal use
#' @export


elo.m3_lik_vect <- function(par, IA_data, all_ids){
  
  k <- par[1]
  
  # init_elo <- par[2:length(par)]
  
  
  # Set intitial elo scores
  currentELO <- c(par[2:length(par)])
  
  names(currentELO) <- all_ids
  
  # Initialize the log likelihood
  L <- 0
  
  
  apply(X = IA_data, MARGIN = 1, function(x){     
    
    # calculate predited winning probablity of the winner
    p_win <- 1/(1+exp(-.01*(currentELO[x[2]] - currentELO[x[3]])))
    
    # Calculation of new ELO scores
    
    currentELO[x[2]] <<- currentELO[x[2]] + exp(k) * (1 - p_win)
    currentELO[x[3]] <<- currentELO[x[3]] - exp(k) * (1 - p_win)
    
    # Update log likelihood
    L <<- L + log(p_win)    
  } )
  
  return(-1*L)
}

