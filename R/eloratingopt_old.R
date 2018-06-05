#' @title Create daily elo ranks and multiple derivatives
#' @description Conducts optimized elo rating analyses as per Foerster, Franz et al
#'   and saves raw, normalized, cardinal, and  categorical ranks in output file.  OLD, DON'T
#'   USE!!!
#' @usage eloratingopt_matrix(agofile, presencefile, sex, dateformat, outputfile)
#' @param agofile Input csv file with dominance interactions, should only contain Date, Winner, Loser
#' @param presencefile Input csv file with date in rows, individuals in columns, and 0/1 for absent/present
#' @param sex Whether data are for males "M" or females "F", no default
#' @param dateformat Need to specify what format dates in data file are, see example
#' @param outputfile Name of file to save ranks to
#' @examples
#' #setwd(choose.dir()) # Interactivevly choose a directory where your input files are located
#'
#' #eloratingopt(agofile="males_ago.csv", presencefile="male_presence.csv", 
#'  #  dateformat="%Y-%m-%d", sex="M", outputfile="rank_data_females.csv")
#' @importFrom stats approx ave optim reshape
#' @importFrom utils head read.csv setWinProgressBar tail winProgressBar write.csv
#' @import reshape2
#' @import BAMMtools
#' @import tcltk
#' @import rlist


eloratingopt_matrix <- function(agofile, presencefile, sex, dateformat, outputfile){
  # Get data
  ago <- read.csv(agofile, header = T, stringsAsFactors = F, sep=",")
  ago$Winner <- as.character(ago$Winner)
  ago$Loser <- as.character(ago$Loser)
  ago$Date <- as.character(ago$Date)
  
  presence <- read.csv(presencefile, header = T, check.names=F, stringsAsFactors = F, sep=",")
  # presence <- read.csv("male_presence.csv", header = T, check.names=F, stringsAsFactors = F, sep=",")
  presence$Date <- as.character(presence$Date)
  presence = presence[complete.cases(presence),]
  
  # Vector with all individuals
  all_inds <- colnames(presence[2:ncol(presence)])
  
  ########################################################################
  #
  #  Make sure ago and presence have same start and end dates
  #
  ########################################################################
  
  presence$Date = lubridate::mdy(presence$Date)
  ago$Date = lubridate::mdy(ago$Date)
  
  tail(ago$Date)
  tail(presence$Date)
  
  min(presence$Date) == min(ago$Date)
  max(presence$Date) == max(ago$Date)
  
  ago = ago[lubridate::year(ago$Date) <= 2013, ] #because there are two seemingly mislabeled 2014 pgs
  
  presence = presence[presence$Date >= min(ago$Date) & presence$Date <= max(ago$Date),]
  
  presence$Date = as.character(presence$Date)
  ago$Date = as.character(ago$Date)
  
  ########################################################################
  #
  #   Filter individuals who do not have at least one win and one loss
  #
  ########################################################################
  # Unique individuals in aggression and pantgrunt data
  indsagg <- unique(c(ago$Winner, ago$Loser))
  
  # Iteration 1
  # Count wins and losses
  win_counts <- rep(0, length(indsagg))
  loss_counts <- rep(0, length(indsagg))
  for (i in 1:length(indsagg)) {
    win_counts[i] <- sum(ago$Winner==indsagg[i])
    loss_counts[i] <- sum(ago$Loser==indsagg[i])
  }
  indsagg <- indsagg[(win_counts>0)*(loss_counts>0)==1]  # keep everybody with at least one win and one loss
  # Get subset of agonistic data that only contains selected individuals
  select <- rep(F, nrow(ago))
  for (i in 1:nrow(ago)) select[i] <- sum(ago$Winner[i]==indsagg)>0 & sum(ago$Loser[i]==indsagg)>0
  # Filtered dataset
  ago <- ago[select,]
  
  # Iteration 2
  # Calculate interaction counts
  ago$CountWinner <- 0
  ago$CountLoser <- 0
  for (i in 1:nrow(ago))
  {
    ago$CountWinner[i] <- sum(ago$Winner[1:i] == ago$Winner[i]) + sum(ago$Loser[1:i] == ago$Winner[i])
    ago$CountLoser[i] <- sum(ago$Winner[1:i] == ago$Loser[i]) + sum(ago$Loser[1:i] == ago$Loser[i])
  }
  # Re-check counts after filtering
  win_counts <- rep(0, length(indsagg))
  loss_counts <- rep(0, length(indsagg))
  for(i in 1:length(indsagg)){
    win_counts[i] <- sum(ago$Winner==indsagg[i])
    loss_counts[i] <- sum(ago$Loser==indsagg[i])
  }
  indsagg2 <- indsagg[(win_counts>0)*(loss_counts>0)==1]
  # Get subset of agonistic data that only contains selected individuals
  select <- rep(F, nrow(ago))
  for (i in 1:nrow(ago)) select[i] <- sum(ago$Winner[i]==indsagg2)>0 & sum(ago$Loser[i]==indsagg2)>0
  # Filtered dataset
  ago <- ago[select,]
  
  # Iteration 3
  # Re-check counts after filtering
  win_counts <- rep(0, length(indsagg2))
  loss_counts <- rep(0, length(indsagg2))
  for(i in 1:length(indsagg2)){
    win_counts[i] <- sum(ago$Winner==indsagg2[i])
    loss_counts[i] <- sum(ago$Loser==indsagg2[i])
  }
  indsagg3 <- indsagg2[(win_counts>0)*(loss_counts>0)==1]
  inds <- indsagg3
  
  # Check number of individuals included in each filtered set
  print(length(indsagg))
  print(length(indsagg2))
  print(length(indsagg3))
  
  
  ##################################
  #
  #   Fit models
  #
  ##################################
  # For testing only
  # Assign existing log table for males to model_log
  ago$Date <- strptime(as.character(ago$Date), dateformat)
  ago$Date <- format(ago$Date, "%Y-%m-%d")
  
  if(sex=="M"){
    # Model 1 (for males)
    model <- optim(par=5, burn_in=100, elo.model1, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='Brent', lower=-10, upper=10)
    model_log <- elo.model1(par=model$par, burn_in=100, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    # model <<- res_m_model1
    # model_log <<- res_m_model1_log
    pred_accuracy <- mean(model_log$elo_w_before[101:nrow(model_log)] > model_log$elo_l_before[101:nrow(model_log)])
  } else if(sex=="F") {
    # Model 3 (for females)
    # model <- optim(par=c(5, rep(0, length(all_inds))), elo.model3, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='BFGS', control = list(maxit = 10000, reltol=1e-10))
    # model_log <- elo.model3(par=model$par, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    model <- res_fem_model3
    model_log <- res_fem_model3_log
    pred_accuracy <- mean(model_log$elo_w_before > model_log$elo_l_before)
  }
  
  # Format dates
  model_log$Date <- as.character(model_log$Date)
  
  
  ######################################
  ######################################
  #
  # Post-processing of elo scores
  #
  # ====================================
  # Step 1:  Normalize elo scores by date
  
  # For females, start here
  if(sex=="F"){
    # Get elo from log object
    elo <- c(model$par[2: length(model$par)])
    names(elo) <- inds
    # Reassign the elo log table to norm to keep the original
    norm <- model_log
    # Create a matrix of elo scores repeated for each date across all individuals in presence data
    elo_matrix <<- t(matrix(rep((elo),nrow(presence)), ncol=nrow(presence)))
    # Filter presence matrix to include only those IDs that are present in filtered indsagg2 table
    presence_matrix <<- presence[,inds]
    # Now fill in NA's for zero in presence matrix
    presence_matrix[presence_matrix==0] <- NA
    # Now multiply elo matrix with all non-NA's in presence matrix
    elo_presence <<- elo_matrix * presence_matrix
    # Convert to dataframe and add date column back
    elo_presence_df <<- as.data.frame(cbind(presence$Date, elo_presence))
    colnames(elo_presence_df)[[1]] <<- "Date"
    # Convert Date into proper date format
    elo_presence_df$Date <<- as.character(elo_presence_df$Date)
    elo_presence_df$Date <<- strptime(as.character(elo_presence_df$Date), "%d.%m.%Y")
    # Go through each row of the dataframe and normalize across all elo scores (ignore NA's)
    elo_norm_presence <<- elo_presence_df
    pb <- winProgressBar(title="Normalizing elo...", label="0%", min=0, max=100, initial=0)
    for(i in 1:nrow(elo_presence_df)){
      elo_norm_presence[i, 2:ncol(elo_norm_presence)] <<- (elo_presence_df[i,2:ncol(elo_presence_df)] - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))/(max(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T) - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))
      info <- sprintf("%d%% completed", round((i/(length(elo_presence_df[[1]]))*100)))
      setWinProgressBar(pb, i/(length(elo_presence_df[[1]]))*100, label=info)
    }
    close(pb)
  } else if(sex=="M") {
    # Reformat elo-after scores of winners and losers into long format
    df <- model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Fomrat columns
    df2 <- seq_long[,c(1,3,4)]
    df2$Date <-strptime(df2$Date, format=dateformat)
    # df2$Individual <- as.numeric(df2$Individual)
    # Order by date and ID
    df2 <- df2[order(df2$Date, df2$Individual),]
    df2$Date <- as.character(df2$Date)
    # Use max achieved score per day
    # df2_daylast <- ddply(df2, .(Date, Individual), summarize, EloScoreAfterLAST=last(EloScoreAfter))
    df2_daymax <- plyr::ddply(df2, plyr::.(Date, Individual), plyr::summarize, EloScoreAfterMax=max(EloScoreAfter))
    df2_daymax$Date <- strptime(df2_daymax$Date, format="%Y-%m-%d")
    # Split dataset by male ID
    elobyid <- split(df2_daymax, df2_daymax$Individual)
    # Filter list objects with only 1 elo score (no interpolation possible)
    # Note: for males this is not necessary as last score is carried forward to end of presence
    # elobyid_f <- list.filter(elobyid, length(which(EloScoreAfterMax!="NA"))>1)
    # Turn presence data into long form to match with interpolated elo scores below
    presence_long <- melt(presence, id.vars = "Date")
    colnames(presence_long) <- c("Date", "Individual", "Presence")
    presence_long$Date <- strptime(presence_long$Date, format=dateformat)
    # Get complete list of dates from presence matrix and add scores for matching dates
    elobyid_full <- list()
    for(i in seq_along(elobyid)){
      elobyid_full[[i]] <- as.data.frame(presence[[1]])
      colnames(elobyid_full[[i]])[1] <- "Date"
      elobyid_full[[i]]$Date <- strptime(elobyid_full[[i]]$Date, format=dateformat)
      # fill in ID info from first row of list object, for length of calendar
      elobyid_full[[i]]$Individual <- rep(elobyid[[i]][1,]$Individual, length(elobyid_full[[i]]$Date))
      # pull in elo scores based on date in elo list objects
      elobyid_full[[i]]$EloScoreAfter <- elobyid[[i]]$EloScoreAfterMax[match(elobyid_full[[i]]$Date, elobyid[[i]]$Date)]
      # merge with presence data table
      elobyid_full[[i]] <- merge(elobyid_full[[i]], presence_long, by.x=c('Individual', 'Date'), by.y=c('Individual', 'Date'))
      # keep only those calendar dates where individual was present as adult
      # elobyid_full[[i]] <- elobyid_full[[i]][elobyid_full[[i]]$Presence==1,]
    }
    
    # Interpolate missing values; carry last score to the end of presence
    interpol_byid <- list()
    elobyid_filled <- list()
    for(i in seq_along(elobyid_full)){
      # Do the interpolation
      interpol_byid[[i]] <- approx(elobyid_full[[i]]$EloScoreAfter, xout=1:length(elobyid_full[[i]]$EloScoreAfter), rule=1:2, method="constant")
      # Fill in interpolated values in orginal dataframe with calendar dates by individual
      interpol_byid[[i]] <- melt(interpol_byid[[i]])
      interpol_byid[[i]] <- interpol_byid[[i]][interpol_byid[[i]]$L1=="y",]
      elobyid_filled[[i]] <- cbind(interpol_byid[[i]], elobyid_full[[i]])[, c(3,4,5,6,1)]
      colnames(elobyid_filled[[i]])[5] <- c('EloScoreAfterInterpolated')
      elobyid_filled[[i]]$Date <- as.character(elobyid_filled[[i]]$Date)
      elobyid_filled[[i]][elobyid_filled[[i]]$Presence==0,]$EloScoreAfterInterpolated <- NA
    }
    # Combine list objects
    elo_data <- melt(elobyid_filled, id.vars=c("Individual", "Date", "Presence"), variable.name="EloType", value.name="EloScore")[,1:5]
    elo_data2 <- elo_data[elo_data$EloType=="EloScoreAfterInterpolated",]
    # elo_data2$Individual <- as.numeric(elo_data2$Individual)
    
    # Fill presence matrix with elo scores from long list object elo_data2
    elo_data2 <- elo_data2[order(elo_data2$Individual,elo_data2$Date),]
    row.names(elo_data2) = NULL
    table(elo_data2$Individual)
    elo_matrix <- matrix(elo_data2$EloScore, ncol=length(unique(elo_data2$Individual))) # ***** wrong number of columns, was ncol(presence)-1
    tail(elo_matrix)
    elo_df <- as.data.frame(elo_matrix)
    colnames(elo_df) <- sort(inds)
    
    # Filter presence matrix to include only those IDs that are present in filtered indsagg2 table
    presence_matrix <- presence[,sort(inds)]
    presence_matrix[presence_matrix==0] <- NA
    elo_presence <- elo_matrix * presence_matrix
    elo_presence_df <- as.data.frame(cbind(presence$Date, elo_presence))
    colnames(elo_presence_df)[[1]] <- "Date"
    elo_presence_df$Date <- as.character(elo_presence_df$Date)
    elo_presence_df$Date <- strptime(as.character(elo_presence_df$Date), dateformat)
    
    # Go through each row of the dataframe and normalize across all elo scores (ignore NA's)
    elo_norm_presence <- elo_presence_df
    pb <- winProgressBar(title="Normalizing elo...", label="0%", min=0, max=100, initial=0)
    for(i in 1:nrow(elo_presence_df)){
      elo_norm_presence[i, 2:ncol(elo_norm_presence)] <- (elo_presence_df[i,2:ncol(elo_presence_df)] - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))/(max(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T) - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))
      info <- sprintf("%d%% completed", round((i/(length(elo_presence_df[[1]]))*100)))
      setWinProgressBar(pb, i/(length(elo_presence_df[[1]]))*100, label=info)
    }
    close(pb)
  }
  
  # save normalized elo scores by date in long format
  elonorm <- elo_norm_presence
  elonorm$Date <- as.character(elonorm$Date)
  elonorm_long <- melt(elonorm, na.rm=T)
  elonorm_long$Date <- strptime(as.character(elonorm_long$Date), "%Y-%m-%d")
  elonorm_long <- elonorm_long[order(elonorm_long$Date, elonorm_long$variable),]
  
  # save non-normalized elo scores by date in long format
  elo_presence_df$Date <- as.character(elo_presence_df$Date)
  elo_long <- melt(elo_presence_df, na.rm=T)
  elo_long$Date <- strptime(as.character(elo_long$Date), "%Y-%m-%d")
  elo_long <- elo_long[order(elo_long$Date, elo_long$variable),]
  colnames(elo_long) <- c("Date", "Individual", "EloScore")
  # Ordinal ranks by day
  elo_long$rank_ord <- ave(elo_long$EloScore, as.character(elo_long$Date), FUN = function(x) rank(-x, ties.method = "first"))
  length(elo_long[[1]])
  
  # ===================================
  # Step 2: Calculate cardinal ranks
  # rel_card_ranks <- elo_presence_df  # data frame for relative cardinal ranks
  pres_scores <- elo_presence_df[,2:ncol(elo_presence_df)]  # data frame for relative cardinal ranks
  pres_ranks <- pres_scores
  
  pb <- winProgressBar(title="Calculate cardinal ranks...", label="0%", min=0, max=100, initial=0)
  for (i in 1:nrow(pres_ranks))
  {
    for (j in 1:ncol(pres_ranks))
    {
      if (!is.na(pres_ranks[i,j])) # if no NA in this cell
      {
        # get relative elo scores
        rel_elo <- pres_scores[i,j] - pres_scores[i,]
        rel_elo[j] <- NA  # exclude difference to own score
        pres_ranks[i,j] <- sum(1 / (1 + exp(-0.01*rel_elo)), na.rm=T)
      }
    }
    info <- sprintf("%d%% completed", round((i/(length(pres_ranks[[1]]))*100)))
    setWinProgressBar(pb, i/(length(pres_ranks[[1]]))*100, label=info)
  }
  close(pb)  # try to speed this upppppppp!  Maybe use nested apply?
  
  # Relativize ranks
  pres_ranks_rel <- pres_ranks
  pb <- winProgressBar(title="Relativize cardinal ranks...", label="0%", min=0, max=100, initial=0)
  for (i in 1:nrow(pres_ranks)) {
    for (j in 1:ncol(pres_ranks)) {
      if (!is.na(pres_ranks[i,j])) {
        pres_ranks_rel[i,j] <- pres_ranks[i,j]/(sum(!is.na(pres_ranks[i,])) - 1)
      }
    }
    info <- sprintf("%d%% completed", round((i/(length(pres_ranks[[1]]))*100)))
    setWinProgressBar(pb, i/(length(pres_ranks[[1]]))*100, label=info)
  }
  close(pb)
  
  # Add date back in
  pres_ranks_rel_date <- cbind(presence[1], pres_ranks_rel)
  head(pres_ranks_rel_date)
  tail(pres_ranks_rel_date)
  
  rel_card_ranks <- pres_ranks_rel_date
  rel_card_ranks_long <- melt(rel_card_ranks, na.rm=T)
  rel_card_ranks_long$Date <- strptime(as.character(rel_card_ranks_long$Date), dateformat)
  rel_card_ranks_long <- rel_card_ranks_long[order(rel_card_ranks_long$Date, rel_card_ranks_long$variable),]
  head(rel_card_ranks_long)
  length(rel_card_ranks_long[[1]])
  
  # ===========================================================
  # Step 3: Find natural breaks in list of elo scores by day
  
  # Split long format cardinal rank scores data table by date
  elobyday <- split(rel_card_ranks_long, as.character(rel_card_ranks_long$Date))
  
  # Iterate through days and calculate breakpoints, then split into 3 categories
  pb <- winProgressBar(title="Calculate cardinal jenks...", label="0%", min=0, max=100, initial=0)
  for(i in seq_along(elobyday)){
    breaks <- getJenksBreaks(elobyday[[i]]$value,4)
    for(j in seq_along(elobyday[[i]][[3]])){
      if (elobyday[[i]][[3]][[j]] <= breaks[[2]]) {
        elobyday[[i]]$cat3[[j]] <- "low"
      } else if (elobyday[[i]][[3]][[j]] > breaks[[3]]) {
        elobyday[[i]]$cat3[[j]] <- "high"
      } else
        elobyday[[i]]$cat3[[j]] <- "mid"
    }
    info <- sprintf("%d%% completed", round((i/(length(elobyday))*100)))
    setWinProgressBar(pb, i/(length(elobyday))*100, label=info)
  }
  close(pb)
  
  # Combine dates into one table again
  elobyday <- do.call("rbind", elobyday)
  
  # Uncomment following section to calculate jenks based on normalized elo
  # # Split long format normalized elo scores data table by date
  # elobyday2 <- split(elonorm_long, as.character(elonorm_long$Date))
  # # Iterate through days and calculate breakpoints, then split into 3 categories
  # for(i in seq_along(elobyday2)){
  #   breaks <- getJenksBreaks(elobyday2[[i]]$value,4)
  #   for(j in seq_along(elobyday2[[i]][[3]])){
  #     if (elobyday2[[i]][[3]][[j]] <= breaks[[2]]) {
  #       elobyday2[[i]]$cat3[[j]] <- "low"
  #     } else if (elobyday2[[i]][[3]][[j]] > breaks[[3]]) {
  #       elobyday2[[i]]$cat3[[j]] <- "high"
  #     } else
  #       elobyday2[[i]]$cat3[[j]] <- "mid"
  #   }
  # }
  # # Combine dates into one table again
  # elobyday2 <- do.call("rbind", elobyday2)
  
  
  ###############################
  #
  #   Compile everything
  #
  ###############################
  # Check if length of objects is the same
  dim(elo_long); dim(elonorm_long); dim(rel_card_ranks_long); dim(elobyday);
  # all(elo_long[,1:2] == elonorm_long[,1:2] & elo_long[,1:2] == rel_card_ranks_long[,1:2]
  #     & elo_long[,1:2] == elobyday[,1:2])
  rank_data <- cbind(elo_long, elonorm_long$value, rel_card_ranks_long$value, elobyday$cat3)
  colnames(rank_data) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "EloCardinal", "JenksEloCardinal")
  row.names(rank_data) = NULL
  head(rank_data)
  
  if(sex=="F"){
    write.csv(rank_data, outputfile, row.names = F)
  } else if(sex=="M"){
    write.csv(rank_data, outputfile, row.names = F)
  }
}



#' @title Create daily elo ranks and multiple derivatives
#' @description Conducts optimized elo rating analyses as per Foerster, Franz et al
#'   and saves raw, normalized, cardinal, and  categorical ranks in output file.  OLD, DON'T USE!!!
#' @usage eloratingoptR(agofile, presencefile, sex, dateformat, outputfile)
#' @param agofile Input csv file with dominance interactions, should only contain Date, Winner, Loser
#' @param presencefile Input csv file with date in rows, individuals in columns, and 0/1 for absent/present
#' @param sex Whether data are for males "M" or females "F", no default
#' @param dateformat Need to specify what format dates in data file are, see example
#' @param outputfile Name of file to save ranks to
#' @examples
#' #setwd(choose.dir()) # Interactivevly choose a directory where your input files are located
#'
#' #eloratingoptR(agofile="males_ago.csv", presencefile="male_presence.csv", 
#' # dateformat="%Y-%m-%d", sex="M", outputfile="rank_data_females.csv")
#' @importFrom stats approx ave optim reshape
#' @importFrom utils head read.csv setWinProgressBar tail winProgressBar write.csv
#' @import reshape2
#' @import BAMMtools
#' @import tcltk
#' @import rlist
#' @importFrom magrittr "%>%"




eloratingoptR <- function(agofile, presencefile, sex, dateformat, outputfile){
  # Get data
  ago <- read.csv(agofile, header = T, stringsAsFactors = F, sep=",")
  ago$Winner <- as.character(ago$Winner)
  ago$Loser <- as.character(ago$Loser)
  ago$Date <- as.character(ago$Date)
  
  presence <- read.csv(presencefile, header = T, check.names=F, stringsAsFactors = F, sep=",")
  # presence <- read.csv("male_presence.csv", header = T, check.names=F, stringsAsFactors = F, sep=",")
  presence$Date <- as.character(presence$Date)
  presence = presence[complete.cases(presence),]
  
  # Vector with all individuals
  # all_inds <- colnames(presence[2:ncol(presence)])
  
# --------------  Make sure ago and presence have same start and end dates -----------------------------

  presence$Date = lubridate::mdy(presence$Date)
  ago$Date = lubridate::mdy(ago$Date)
  
  tail(ago$Date)
  tail(presence$Date)
  
  min(presence$Date) == min(ago$Date)
  max(presence$Date) == max(ago$Date)
  
  ago = ago[lubridate::year(ago$Date) <= 2013, ] #because there are two seemingly mislabeled 2014 pgs
  
  presence = presence[presence$Date >= min(ago$Date) & presence$Date <= max(ago$Date),]
  
  presence$Date = as.character(presence$Date)
  ago$Date = as.character(ago$Date)
  
# -------------   Filter individuals who do not have at least one win and one loss -----------------

  # Unique individuals in aggression and pantgrunt data
  indsagg <- unique(c(ago$Winner, ago$Loser))
  
  # Iteration 1
  # Count wins and losses
  win_counts <- rep(0, length(indsagg))
  loss_counts <- rep(0, length(indsagg))
  for (i in 1:length(indsagg)) {
    win_counts[i] <- sum(ago$Winner==indsagg[i])
    loss_counts[i] <- sum(ago$Loser==indsagg[i])
  }
  indsagg <- indsagg[(win_counts>0)*(loss_counts>0)==1]  # keep everybody with at least one win and one loss
  # Get subset of agonistic data that only contains selected individuals
  select <- rep(F, nrow(ago))
  for (i in 1:nrow(ago)) select[i] <- sum(ago$Winner[i]==indsagg)>0 & sum(ago$Loser[i]==indsagg)>0
  # Filtered dataset
  ago <- ago[select,]
  
  # Iteration 2
  # Calculate interaction counts
  ago$CountWinner <- 0
  ago$CountLoser <- 0
  for (i in 1:nrow(ago))
  {
    ago$CountWinner[i] <- sum(ago$Winner[1:i] == ago$Winner[i]) + sum(ago$Loser[1:i] == ago$Winner[i])
    ago$CountLoser[i] <- sum(ago$Winner[1:i] == ago$Loser[i]) + sum(ago$Loser[1:i] == ago$Loser[i])
  }
  # Re-check counts after filtering
  win_counts <- rep(0, length(indsagg))
  loss_counts <- rep(0, length(indsagg))
  for(i in 1:length(indsagg)){
    win_counts[i] <- sum(ago$Winner==indsagg[i])
    loss_counts[i] <- sum(ago$Loser==indsagg[i])
  }
  indsagg2 <- indsagg[(win_counts>0)*(loss_counts>0)==1]
  # Get subset of agonistic data that only contains selected individuals
  select <- rep(F, nrow(ago))
  for (i in 1:nrow(ago)) select[i] <- sum(ago$Winner[i]==indsagg2)>0 & sum(ago$Loser[i]==indsagg2)>0
  # Filtered dataset
  ago <- ago[select,]
  
  # Iteration 3
  # Re-check counts after filtering
  win_counts <- rep(0, length(indsagg2))
  loss_counts <- rep(0, length(indsagg2))
  for(i in 1:length(indsagg2)){
    win_counts[i] <- sum(ago$Winner==indsagg2[i])
    loss_counts[i] <- sum(ago$Loser==indsagg2[i])
  }
  indsagg3 <- indsagg2[(win_counts>0)*(loss_counts>0)==1]
  inds <- indsagg3
  
  # Check number of individuals included in each filtered set
  print(length(indsagg))
  print(length(indsagg2))
  print(length(indsagg3))
  
  indsagg3
  
  presence = presence[,colnames(presence) %in% c("Date", indsagg3)]
  
  row.names(presence) = NULL
  row.names(ago) = NULL
  
  max(lubridate::ymd(ago$Date)) == max(lubridate::ymd(presence$Date))
  min(lubridate::ymd(ago$Date)) == min(lubridate::ymd(presence$Date))
  
  all_inds <- colnames(presence[2:ncol(presence)])
  
# --------------   Fit models ----------------------------

  # For testing only
  # Assign existing log table for males to model_log
  ago$Date <- strptime(as.character(ago$Date), dateformat)
  ago$Date <- format(ago$Date, "%Y-%m-%d")
  
  if(sex=="M"){
    # Model 1 (for males)
    model <- optim(par=5, burn_in=100, elo.model1, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='Brent', lower=-10, upper=10)
    model_log <- elo.model1(par=model$par, burn_in=100, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    # model <<- res_m_model1
    # model_log <<- res_m_model1_log
    pred_accuracy <- mean(model_log$elo_w_before[101:nrow(model_log)] > model_log$elo_l_before[101:nrow(model_log)])
  } else if(sex=="F") {
    # Model 3 (for females)
    # model <- optim(par=c(5, rep(0, length(all_inds))), elo.model3, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='BFGS', control = list(maxit = 10000, reltol=1e-10))
    # model_log <- elo.model3(par=model$par, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    model <- res_fem_model3
    model_log <- res_fem_model3_log
    pred_accuracy <- mean(model_log$elo_w_before > model_log$elo_l_before)
  }
  
  # Format dates
  model_log$Date <- as.character(model_log$Date)
  
# =================== Post-processing of elo scores ==========================
# ------------- Step 1:  Normalize elo scores by date --------------------------
  
  # For females, start here
  if(sex=="F"){
    # Get elo from log object
    elo <- c(model$par[2: length(model$par)])
    names(elo) <- inds
    # Reassign the elo log table to norm to keep the original
    norm <- model_log
    # Create a matrix of elo scores repeated for each date across all individuals in presence data
    elo_matrix <<- t(matrix(rep((elo),nrow(presence)), ncol=nrow(presence)))
    # Filter presence matrix to include only those IDs that are present in filtered indsagg2 table
    presence_matrix <<- presence[,inds]
    # Now fill in NA's for zero in presence matrix
    presence_matrix[presence_matrix==0] <- NA
    # Now multiply elo matrix with all non-NA's in presence matrix
    elo_presence <<- elo_matrix * presence_matrix
    # Convert to dataframe and add date column back
    elo_presence_df <<- as.data.frame(cbind(presence$Date, elo_presence))
    colnames(elo_presence_df)[[1]] <<- "Date"
    # Convert Date into proper date format
    elo_presence_df$Date <<- as.character(elo_presence_df$Date)
    elo_presence_df$Date <<- strptime(as.character(elo_presence_df$Date), "%d.%m.%Y")
    # Go through each row of the dataframe and normalize across all elo scores (ignore NA's)
    elo_norm_presence <<- elo_presence_df
    pb <- winProgressBar(title="Normalizing elo...", label="0%", min=0, max=100, initial=0)
    for(i in 1:nrow(elo_presence_df)){
      elo_norm_presence[i, 2:ncol(elo_norm_presence)] <<- (elo_presence_df[i,2:ncol(elo_presence_df)] - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))/(max(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T) - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))
      info <- sprintf("%d%% completed", round((i/(length(elo_presence_df[[1]]))*100)))
      setWinProgressBar(pb, i/(length(elo_presence_df[[1]]))*100, label=info)
    }
    close(pb)
  } else if(sex=="M") {
    # Reformat elo-after scores of winners and losers into long format
    df <- model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Fomrat columns
    df2 <- seq_long[,c(1,3,4)]
    df2$Date <-strptime(df2$Date, format=dateformat)
    # df2$Individual <- as.numeric(df2$Individual)
    # Order by date and ID
    df2 <- df2[order(df2$Date, df2$Individual),]
    df2$Date <- as.character(df2$Date)
    # Use max achieved score per day
    # df2_daylast <- ddply(df2, .(Date, Individual), summarize, EloScoreAfterLAST=last(EloScoreAfter))
    df2_daymax <- plyr::ddply(df2, plyr::.(Date, Individual), plyr::summarize, EloScoreAfterMax=max(EloScoreAfter))
    df2_daymax$Date <- strptime(df2_daymax$Date, format="%Y-%m-%d")
    # Split dataset by male ID
    elobyid <- split(df2_daymax, df2_daymax$Individual)
    # Filter list objects with only 1 elo score (no interpolation possible)
    # Note: for males this is not necessary as last score is carried forward to end of presence
    # elobyid_f <- list.filter(elobyid, length(which(EloScoreAfterMax!="NA"))>1)
    # Turn presence data into long form to match with interpolated elo scores below
    presence_long <- melt(presence, id.vars = "Date")
    colnames(presence_long) <- c("Date", "Individual", "Presence")
    presence_long$Date <- strptime(presence_long$Date, format=dateformat)
    # Get complete list of dates from presence matrix and add scores for matching dates
    elobyid_full <- list()
    for(i in seq_along(elobyid)){
      elobyid_full[[i]] <- as.data.frame(presence[[1]])
      colnames(elobyid_full[[i]])[1] <- "Date"
      elobyid_full[[i]]$Date <- strptime(elobyid_full[[i]]$Date, format=dateformat)
      # fill in ID info from first row of list object, for length of calendar
      elobyid_full[[i]]$Individual <- rep(elobyid[[i]][1,]$Individual, length(elobyid_full[[i]]$Date))
      # pull in elo scores based on date in elo list objects
      elobyid_full[[i]]$EloScoreAfter <- elobyid[[i]]$EloScoreAfterMax[match(elobyid_full[[i]]$Date, elobyid[[i]]$Date)]
      # merge with presence data table
      elobyid_full[[i]] <- merge(elobyid_full[[i]], presence_long, by.x=c('Individual', 'Date'), by.y=c('Individual', 'Date'))
      # keep only those calendar dates where individual was present as adult
      # elobyid_full[[i]] <- elobyid_full[[i]][elobyid_full[[i]]$Presence==1,]
    }
    
    # Interpolate missing values; carry last score to the end of presence
    interpol_byid <- list()
    elobyid_filled <- list()
    for(i in seq_along(elobyid_full)){
      # Do the interpolation
      interpol_byid[[i]] <- approx(elobyid_full[[i]]$EloScoreAfter, xout=1:length(elobyid_full[[i]]$EloScoreAfter), rule=1:2, method="constant")
      # Fill in interpolated values in orginal dataframe with calendar dates by individual
      interpol_byid[[i]] <- melt(interpol_byid[[i]])
      interpol_byid[[i]] <- interpol_byid[[i]][interpol_byid[[i]]$L1=="y",]
      elobyid_filled[[i]] <- cbind(interpol_byid[[i]], elobyid_full[[i]])[, c(3,4,5,6,1)]
      colnames(elobyid_filled[[i]])[5] <- c('EloScoreAfterInterpolated')
      elobyid_filled[[i]]$Date <- as.character(elobyid_filled[[i]]$Date)
      elobyid_filled[[i]][elobyid_filled[[i]]$Presence==0,]$EloScoreAfterInterpolated <- NA
    }
    # Combine list objects
    elo_data <- melt(elobyid_filled, id.vars=c("Individual", "Date", "Presence"), variable.name="EloType", value.name="EloScore")[,1:5]
    elo_data2 <- elo_data[elo_data$EloType=="EloScoreAfterInterpolated",]
    # elo_data2$Individual <- as.numeric(elo_data2$Individual)
    
    # Fill presence matrix with elo scores from long list object elo_data2
    elo_data2 <- elo_data2[order(elo_data2$Individual,elo_data2$Date),]
    row.names(elo_data2) = NULL
    elo_data2 = dplyr::filter(elo_data2, !is.na(EloScore))
    # table(elo_data2$Individual)
    # elo_matrix <- matrix(elo_data2$EloScore, ncol=length(unique(elo_data2$Individual))) # ***** wrong number of columns, was ncol(presence)-1
    # tail(elo_matrix)
    # elo_df <- as.data.frame(elo_matrix)
    # colnames(elo_df) <- sort(inds)
    
    # Filter presence matrix to include only those IDs that are present in filtered indsagg2 table
    # presence_matrix <- presence[,sort(inds)]
    # presence_matrix[presence_matrix==0] <- NA
    # elo_presence <- elo_matrix * presence_matrix
    # elo_presence_df <- as.data.frame(cbind(presence$Date, elo_presence))
    # colnames(elo_presence_df)[[1]] <- "Date"
    # elo_presence_df$Date <- as.character(elo_presence_df$Date)
    # elo_presence_df$Date <- strptime(as.character(elo_presence_df$Date), dateformat)
    
    # Go through each row of the dataframe and normalize across all elo scores (ignore NA's)
    # elo_norm_presence <- elo_presence_df
    # pb <- winProgressBar(title="Normalizing elo...", label="0%", min=0, max=100, initial=0)
    # for(i in 1:nrow(elo_presence_df)){
    #   elo_norm_presence[i, 2:ncol(elo_norm_presence)] <- (elo_presence_df[i,2:ncol(elo_presence_df)] - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))/(max(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T) - min(elo_presence_df[i, 2:ncol(elo_presence_df)], na.rm=T))
    #   info <- sprintf("%d%% completed", round((i/(length(elo_presence_df[[1]]))*100)))
    #   setWinProgressBar(pb, i/(length(elo_presence_df[[1]]))*100, label=info)
    # }
    # close(pb)
    
    elo_data2$Date = lubridate::ymd(elo_data2$Date)
    
    elo_long =
      elo_data2 %>%
      dplyr::group_by(Date) %>%
      dplyr::mutate(EloNorm = (EloScore - min(EloScore, na.rm = T))/
                      (max(EloScore, na.rm = T)-min(EloScore, na.rm = T))) %>%
      dplyr::arrange(Date, Individual)
    
  }
  
  # save normalized elo scores by date in long format
  # elonorm <- elo_norm_presence
  # elonorm$Date <- as.character(elonorm$Date)
  # elonorm_long <- melt(elonorm, na.rm=T)
  # elonorm_long$Date <- strptime(as.character(elonorm_long$Date), "%Y-%m-%d")
  # elonorm_long <- elonorm_long[order(elonorm_long$Date, elonorm_long$variable),]
  
  # save non-normalized elo scores by date in long format
  # elo_presence_df$Date <- as.character(elo_presence_df$Date)
  # elo_long <- melt(elo_presence_df, na.rm=T)
  # elo_long$Date <- strptime(as.character(elo_long$Date), "%Y-%m-%d")
  # elo_long <- elo_long[order(elo_long$Date, elo_long$variable),]
  # colnames(elo_long) <- c("Date", "Individual", "EloScore")
  # elo_long =
  #   elo_data2 %>%
  #   filter(!is.na(EloScore)) %>%
  #   arrange(Date, Individual) %>%
  #   as.data.frame()#shit yeah, same thing
  
  # Ordinal ranks by day
  elo_long$rank_ord <- ave(elo_long$EloScore, as.character(elo_long$Date), FUN = function(x) rank(-x, ties.method = "first"))
  # length(elo_long[[1]])
  
# ----------- Step 2: Calculate cardinal ranks -----------------------------
  # rel_card_ranks <- elo_presence_df  # data frame for relative cardinal ranks
  # pres_scores <- elo_presence_df[,2:ncol(elo_presence_df)]  # data frame for relative cardinal ranks
  # pres_ranks <- pres_scores
  
  # pb <- winProgressBar(title="Calculate cardinal ranks...", label="0%", min=0, max=100, initial=0)
  # for (i in 1:nrow(pres_ranks))
  # {
  #   for (j in 1:ncol(pres_ranks))
  #   {
  #     if (!is.na(pres_ranks[i,j])) # if no NA in this cell
  #     {
  #       # get relative elo scores
  #       rel_elo <- pres_scores[i,j] - pres_scores[i,]
  #       rel_elo[j] <- NA  # exclude difference to own score
  #       pres_ranks[i,j] <- sum(1 / (1 + exp(-0.01*rel_elo)), na.rm=T)
  #     }
  #   }
  #   info <- sprintf("%d%% completed", round((i/(length(pres_ranks[[1]]))*100)))
  #   setWinProgressBar(pb, i/(length(pres_ranks[[1]]))*100, label=info)
  # }
  # close(pb)  # try to speed this upppppppp!  Maybe use nested apply?
  
  # elo_long$Date = as.POSIXct(elo_long$Date)
  
  # cardinalize = function(x){
  #   newish = c()
  #   for(i in 1:length(x)){
  #       # if(is.na(x[i])){
  #       #   newish[i] = NA
  #       #   next
  #       # }
  #       rel_elo = x[i] - x
  #       newish[i] = sum(1 / (1 + exp(-0.01*rel_elo)), na.rm=T) - .5
  #   }
  #   return(newish)
  # }
  
  cardinalize = function(x){
    carddat = c()
    carddat = sapply(x, function(y){
      sum(1 / (1 + exp(-0.01*(y - x))), na.rm=T) - .5 #subtracting .5 is equivalent to removing the prob of winning against oneself
      #b/c 1/(1 + exp(-0.01*0)) = 1/(1 + exp(0)) = 1/(1 + 1) = 1/2
    })
    return(carddat)
  }
  
  relativize = function(x){
    reldat = c()
    reldat = sapply(x, function(y){
      y/(length(x) - 1)
    })
    return(reldat)
  }
  
  elo_long =
    elo_long %>%
    dplyr::group_by(Date) %>%
    dplyr::mutate(pct_beaten = cardinalize(EloScore),
                  elo_rel = relativize(pct_beaten))
  
  # Relativize ranks
  # pres_ranks_rel <- pres_ranks
  # pb <- winProgressBar(title="Relativize cardinal ranks...", label="0%", min=0, max=100, initial=0)
  # for (i in 1:nrow(pres_ranks)) {
  #   for (j in 1:ncol(pres_ranks)) {
  #     if (!is.na(pres_ranks[i,j])) {
  #       pres_ranks_rel[i,j] <- pres_ranks[i,j]/(sum(!is.na(pres_ranks[i,])) - 1)
  #     }
  #   }
  #   info <- sprintf("%d%% completed", round((i/(length(pres_ranks[[1]]))*100)))
  #   setWinProgressBar(pb, i/(length(pres_ranks[[1]]))*100, label=info)
  # }
  # close(pb)
  
  # Add date back in
  # pres_ranks_rel_date <- cbind(presence[1], pres_ranks_rel)
  # head(pres_ranks_rel_date)
  # tail(pres_ranks_rel_date)
  #
  # rel_card_ranks <- pres_ranks_rel_date
  # rel_card_ranks_long <- melt(rel_card_ranks, na.rm=T)
  # rel_card_ranks_long$Date <- strptime(as.character(rel_card_ranks_long$Date), dateformat)
  # rel_card_ranks_long <- rel_card_ranks_long[order(rel_card_ranks_long$Date, rel_card_ranks_long$variable),]
  # head(rel_card_ranks_long)
  # length(rel_card_ranks_long[[1]])
  
# ------------- Step 3: Find natural breaks in list of elo scores by day ----------------------------
  
  # Split long format cardinal rank scores data table by date
  # elobyday <- split(rel_card_ranks_long, as.character(rel_card_ranks_long$Date))
  #
  # # Iterate through days and calculate breakpoints, then split into 3 categories
  # pb <- winProgressBar(title="Calculate cardinal jenks...", label="0%", min=0, max=100, initial=0)
  # for(i in seq_along(elobyday)){
  #   breaks <- getJenksBreaks(elobyday[[i]]$value,4)
  #   for(j in seq_along(elobyday[[i]][[3]])){
  #     if (elobyday[[i]][[3]][[j]] <= breaks[[2]]) {
  #       elobyday[[i]]$cat3[[j]] <- "low"
  #     } else if (elobyday[[i]][[3]][[j]] > breaks[[3]]) {
  #       elobyday[[i]]$cat3[[j]] <- "high"
  #     } else
  #       elobyday[[i]]$cat3[[j]] <- "mid"
  #   }
  #   info <- sprintf("%d%% completed", round((i/(length(elobyday))*100)))
  #   setWinProgressBar(pb, i/(length(elobyday))*100, label=info)
  # }
  # close(pb)
  
  jenksify = function(x){
    breaks = getJenksBreaks(x, 4)
    # cats = c()
    cats = ifelse(x <= breaks[[2]], "low",
                  ifelse(x > breaks[[3]], "high", "mid"))
    return(cats)
  }
  
  elo_long =
    elo_long %>%
    dplyr::group_by(Date) %>%
    dplyr::mutate(JenksEloCardinal = jenksify(elo_rel))
  
  
  # Combine dates into one table again
  # elobyday <- do.call("rbind", elobyday)
  
  # Uncomment following section to calculate jenks based on normalized elo
  # # Split long format normalized elo scores data table by date
  # elobyday2 <- split(elonorm_long, as.character(elonorm_long$Date))
  # # Iterate through days and calculate breakpoints, then split into 3 categories
  # for(i in seq_along(elobyday2)){
  #   breaks <- getJenksBreaks(elobyday2[[i]]$value,4)
  #   for(j in seq_along(elobyday2[[i]][[3]])){
  #     if (elobyday2[[i]][[3]][[j]] <= breaks[[2]]) {
  #       elobyday2[[i]]$cat3[[j]] <- "low"
  #     } else if (elobyday2[[i]][[3]][[j]] > breaks[[3]]) {
  #       elobyday2[[i]]$cat3[[j]] <- "high"
  #     } else
  #       elobyday2[[i]]$cat3[[j]] <- "mid"
  #   }
  # }
  # # Combine dates into one table again
  # elobyday2 <- do.call("rbind", elobyday2)
  
  
# ---------------   Compile everything ----------------------

  # Check if length of objects is the same
  # dim(elo_long); dim(elonorm_long); dim(rel_card_ranks_long); dim(elobyday);
  # all(elo_long[,1:2] == elonorm_long[,1:2] & elo_long[,1:2] == rel_card_ranks_long[,1:2]
  #     & elo_long[,1:2] == elobyday[,1:2])
  # rank_data <- cbind(elo_long, elonorm_long$value, rel_card_ranks_long$value, elobyday$cat3)
  elo_long =
    elo_long %>%
    dplyr::select(Date, Individual, EloScore, rank_ord, EloNorm, pct_beaten, elo_rel, JenksEloCardinal) %>%
    as.data.frame()
  
  colnames(elo_long) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "ExpNumBeaten", "EloCardinal", "JenksEloCardinal")
  
  head(elo_long)
  
  if(sex=="F"){
    write.csv(elo_long, outputfile, row.names = F)
  } else if(sex=="M"){
    write.csv(elo_long, outputfile, row.names = F)
  }
}



