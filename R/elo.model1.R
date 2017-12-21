#' @title Optimize k paramter in Elo rating method
#' @description Function to optimize k parameter in Elo Rating Method
#' @usage elo.model1(par, burn_in=100, init_elo = 1000, IA_data, all_ids, return_likelihood = T)
#' @param par ???
#' @param burn_in ???
#' @param init_elo ???
#' @param IA_data ???
#' @param all_ids ???
#' @param return_likelihood ???
#' @examples
#' #setwd(choose.dir()) # Interactivevly choose a directory where your input files are located
#'
#' #fill in later
#' @export
elo.model1 <- function(par, burn_in=100, init_elo = 1000, IA_data, all_ids, return_likelihood = T)
{ 
  k <- par
  
  # Initialize output columns
  if (!return_likelihood) IA_data$elo_l_before <- IA_data$elo_w_before <- IA_data$elo_l_after <- IA_data$elo_w_after <- NA
  
  # Set intitial elo scores
  currentELO <- rep(init_elo,length(all_ids))
  names(currentELO) <- all_ids
  
  # Initialize the log likelihood
  L <- 0
  
  # Start loop
  for(i in 1:nrow(IA_data))   
  {     
    ind1 <- which(names(currentELO)==IA_data$Winner[i])
    ind2 <- which(names(currentELO)==IA_data$Loser[i])
    
    if (!return_likelihood) 
    {
      IA_data$elo_w_before[i] <- currentELO[ind1]
      IA_data$elo_l_before[i] <- currentELO[ind2]
    }
    
    # calculate predited winning probablity of the winner
    p_win <- 1/(1+exp(-.01*(currentELO[ind1] - currentELO[ind2])))
    
    # Calculation of new ELO scores
    if (i <= burn_in)   # during burn-in period all k values are fixed to 100
    {
      currentELO[ind1] <- currentELO[ind1] + 100 * (1 - p_win)  # new Elo score of the Winner
      currentELO[ind2] <- currentELO[ind2] - 100 * (1 - p_win)  # new Elo score of the Loser
    }
    else  # after the burn-in period fitted k values are used
    {  
      currentELO[ind1] <- currentELO[ind1] + exp(k) * (1 - p_win)  # new Elo score of the Winner
      currentELO[ind2] <- currentELO[ind2] - exp(k) * (1 - p_win)  # new Elo score of the Loser
    }
    
    # write calculated elo scores to output columns
    if (!return_likelihood) 
    {
      IA_data$elo_w_after[i] <- currentELO[ind1]
      IA_data$elo_l_after[i] <- currentELO[ind2]
    }
    
    # Update log likelihood
    if (i > burn_in) L <- L + log(p_win)    
  } 
  
  if (return_likelihood) return(-1*L)
  else return(IA_data)
}