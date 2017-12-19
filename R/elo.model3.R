#' Function to optimize k parameter and entry Elo scores
elo.model3 <- function(par, IA_data, all_ids, return_likelihood = T)
{ 
  k <- par[1]
  
  init_elo <- par[2:length(par)]
  
  # Initialize output columns
  if (!return_likelihood) IA_data$elo_l_before <- IA_data$elo_w_before <- IA_data$elo_l_after <- IA_data$elo_w_after <- NA
  
  # Set intitial elo scores
  currentELO <- c(init_elo)
  
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
    currentELO[ind1] <- currentELO[ind1] + exp(k) * (1 - p_win)  # new Elo score of the Winner
    currentELO[ind2] <- currentELO[ind2] - exp(k) * (1 - p_win)  # new Elo score of the Loser
    
    
    # write calculated elo scores to output columns
    if (!return_likelihood) 
    {
      IA_data$elo_w_after[i] <- currentELO[ind1]
      IA_data$elo_l_after[i] <- currentELO[ind2]
    }
    
    # Update log likelihood
    L <- L + log(p_win)    
  } 
  
  if (return_likelihood) return(-1*L)
  else return(IA_data)
}
