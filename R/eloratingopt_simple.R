#' @title Create daily elo ranks and multiple derivatives (simplified)
#' @description Conducts optimized elo rating analyses as per Foerster, Franz et al
#'   and saves raw, normalized, cardinal, and  categorical ranks in output file.
#' @usage eloratingopt_simple(agon_data, pres_data, sex, outputfile = NULL, returnR = TRUE)
#' @param agon_data Input data frame with dominance interactions, should only contain Date, Winner, Loser
#' @param pres_data Input data frame with columns "id", "start_date" and "end_date".  Date
#'   columns should be formatted as MONTH/DAY/YEAR, or already as Date class
#' @param sex Whether data are for males "M" or females "F", no default
#' @param outputfile Name of csv file to save ranks to.  Default is NULL, in which case 
#'   the function will only return a table in R.  If you supply an output file name
#'   the function will save the results as a csv file
#' @param returnR whether to return an R object from the function call.  Default is TRUE
#' @examples
#'
#' #eloratingoptR(agon_data="males_ago.csv", pres_data="male_presence.csv", 
#' #  dateformat="%Y-%m-%d", sex="M", outputfile="rank_data_females.csv")
#' @export
#' @importFrom stats approx ave optim reshape
#' @importFrom utils head read.csv setWinProgressBar tail winProgressBar write.csv
#' @import reshape2
#' @import BAMMtools
#' @import tcltk
#' @import rlist
#' @importFrom magrittr "%>%"




eloratingopt_simple <- function(agon_data, pres_data, sex, outputfile = NULL, returnR = TRUE){
  # Get data
  
  if(length(outputfile) == 0 & returnR == FALSE){
    stop("supply an outputfile name or set returnR to TRUE (or both)")
  }
  
  ago = agon_data
  names(ago) <- tools::toTitleCase(tolower(names(ago)))
  if(!all(names(ago) %in% c("Date", "Winner", "Loser"))){
    stop("colnames in agonistic data should be 'Date', 'Winner', 'Loser' (not case sensitive)")
  }
  if(class(ago$Date) != "Date"){
    ago$Date = lubridate::mdy(ago$Date)
  }
  
  
  presence <- pres_data
  names(presence) = tolower(names(presence))
  if(!all(names(presence) %in% c("id", "start_date", "end_date"))){
    stop("colnames in presence data should be 'id', 'start_date', 'end_date' (not case sensitive)")
  }
  if(class(pres_data$start_date) != "Date"){
    pres_data$start_date = lubridate::mdy(pres_data$start_date)}
  if(class(pres_data$end_date) != "Date"){
    pres_data$end_date = lubridate::mdy(pres_data$end_date)}
  
  presence$id = as.character(presence$id)
  
  
  # ---------------  Make sure ago and presence have same start and end dates ----------------------
  # **** maybe we don't need to do this anymore??? 
  # or maybe move this after the filtering in case the ago file is altered much
  
  # min(presence$start_date) == min(ago$Date)
  # min(presence$start_date)
  # min(ago$Date)
  
  presence$start_date[presence$start_date < min(ago$Date)] = min(ago$Date)
  
  if(max(presence$end_date) > (max(ago$Date) + lubridate::days(31))){
    warning("careful, ranks extend more than a month after final agonistic interaction")
  }
  
  if(nrow(presence[presence$start_date >= presence$end_date,]) > 0){
    
    badid = presence[presence$start_date >= presence$end_date, "id"]
    
    stop(paste("some start dates are later than end dates. ID's with this problem:\n", paste(badid, collapse = ", ")))
    
  }
  
  
  # ---------- Filter individuals who do not have at least one win and one loss ----------------
  
  presence$wl = 0 #add dummy column to count wins and losses
  
  # vectorized loop to remove individuals from presence and ago data with 0 wins or 0 losses
  repeat{
    
    oldnum = nrow(presence)
    
    presence$wl = sapply(X = presence$id, function(x) sum(ago$Winner == x) * sum(ago$Loser == x))
    
    presence = presence %>% dplyr::filter(wl != 0)
    
    ago = ago %>% 
      dplyr::filter(Winner %in% presence$id & 
                      Loser %in% presence$id)
    
    if(nrow(presence) == oldnum) break
    
  }
  
  presence = presence[,-4] # remove dummy variable
  
  all_inds = sort(presence$id)
  
  
  # ---------------   Fit models  --------------------------------
  
  
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
    model <- optim(par=c(5, rep(0, length(all_inds))), elo.m3_lik_vect, all_ids = all_inds, IA_data = ago, method='BFGS', control = list(maxit = 10000, reltol=1e-10))
    ### USE SAVED "../data prep code/fem_mod_kk_2013.RData" TO SAVE TIME!
    model_log <- elo.model3(par=model$par, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    # model <- res_fem_model3
    # model_log <- res_fem_model3_log
    pred_accuracy <- mean(model_log$elo_w_before > model_log$elo_l_before)
  }
  
  
  # ================ Post-processing of elo scores =================================
  
  # ---------------- Step 1:  Normalize elo scores by date -------------------------
  
  # For females, start here
  if(sex=="F"){
    
    # Get elo from log object
    initelo <- data.frame(Date = presence[order(presence$id), "start_date"],
                          Individual = all_inds,
                          EloScoreAfter = model$par[2: length(model$par)], 
                          stringsAsFactors = F)
    # names(elo) <- all_inds #pretty sure this should be "all_inds", but DOUBLE CHECK!!! (
    # changed from "inds" to "all_inds")
    
    df <- model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Format columns
    df2 <- seq_long[,c(1,3,4)]
    row.names(df2) = NULL
    
    df2 = rbind.data.frame(initelo, df2) #combine starting elo scores with elo scores after interactions
    
  } else if(sex=="M") {
    # Reformat elo-after scores of winners and losers into long format
    df <- model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Format columns
    df2 <- seq_long[,c(1,3,4)]
    
    #skipping step with the male models to add starting elo scores. Thus in male models 
    # individuals start being ranked after their first interaction, whereas
    # in female models individuals start being ranked immediately upon entry.
    
  }
  
  # Order by date and ID
  df2 <- df2[order(df2$Date, df2$Individual),]  
  
  # Use max achieved score per day
  df2_daymax <-
    df2 %>%
    dplyr::group_by(Date, Individual) %>%
    dplyr::summarise(EloScoreAfterMax=max(EloScoreAfter)) %>%
    as.data.frame()
  
  #create list of all days each individual is present
  temp = list()
  for(i in 1:nrow(presence)){
    temp[[i]] = cbind.data.frame(Individual = presence[i, "id"], 
                                 Date = seq(from = presence[i, "start_date"], 
                                            to = presence[i, "end_date"], 
                                            by = 1))
  }
  presence_long <- do.call(rbind.data.frame, temp)
  presence_long$Presence = 1
  presence_long$Individual = as.character(presence_long$Individual)
  
  
  presence_long$Elo = df2_daymax$EloScoreAfterMax[match(paste0(presence_long$Individual, presence_long$Date),
                                                        paste0(df2_daymax$Individual, df2_daymax$Date))]
  
  presence_long = 
    presence_long %>%
    dplyr::group_by(Individual) %>%
    dplyr::mutate(Elo_interpol = approx(Elo, xout = 1:length(Elo), 
                                        rule = 1:2, method = "constant")$y) %>%
    as.data.frame()
  
  elo_data2 = dplyr::filter(presence_long, !is.na(Elo_interpol)) %>% 
    dplyr::select(-Elo, -Presence) %>% 
    dplyr::rename(EloScore = Elo_interpol) %>%
    as.data.frame()
  
  elo_long =
    elo_data2 %>%
    dplyr::group_by(Date) %>%
    dplyr::mutate(EloNorm = (EloScore - min(EloScore, na.rm = T))/
                    (max(EloScore, na.rm = T)-min(EloScore, na.rm = T))) %>%
    dplyr::arrange(Date, Individual)
  
  
  
  # ----------------- Step 2: Ordinal ranks by day -----------------------------------
  
  elo_long$rank_ord <- ave(elo_long$EloScore, as.character(elo_long$Date), FUN = function(x) rank(-x, ties.method = "first"))
  
  # ----------------- Step 3: Calculate cardinal ranks -------------------------------
  
  
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
  
  
  # --------------------- Step 4: Find natural breaks in list of elo scores by day --------------------
  
  jenksify = function(x){
    breaks = BAMMtools::getJenksBreaks(x, 4)
    # cats = c()
    cats = ifelse(x <= breaks[[2]], "low",
                  ifelse(x > breaks[[3]], "high", "mid"))
    return(cats)
  }
  
  elo_long =
    elo_long %>%
    dplyr::group_by(Date) %>%
    dplyr::mutate(JenksEloCardinal = jenksify(elo_rel))
  
  
  # ------------------------ Step 5: prettify -----------------------------------
  elo_long =
    elo_long %>%
    dplyr::select(Date, Individual, EloScore, rank_ord, EloNorm, pct_beaten, elo_rel, JenksEloCardinal) %>%
    as.data.frame()
  
  colnames(elo_long) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "ExpNumBeaten", "EloCardinal", "JenksEloCardinal")
  
  # head(elo_long)
  
  cat(paste0("k = ", round(exp(model$par[1]), 3)))
  cat(paste0("prediction accuracy = ", round(pred_accuracy, 3)))
  
  if(length(outputfile) > 0){
    write.csv(elo_long, outputfile, row.names = F)
  } 
  
  if(returnR == TRUE){
    return(elo_long)
  }
  
}