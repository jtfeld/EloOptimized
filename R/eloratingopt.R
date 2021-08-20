#' @title Create daily ML fitted Elo ranks and multiple derivatives
#' @description Conducts \strong{optimized} elo rating analyses as per Foerster, Franz et al
#'   and outputs raw, normalized, cardinal, and  categorical ranks as a list object in 
#'   R or in an output file. For non-optimized Elo score calculation, use 
#'   \code{\link{eloratingfixed}}.
#' @usage eloratingopt(agon_data, pres_data, fit_init_elo = FALSE, outputfile = NULL, 
#'   returnR = TRUE)
#' @param agon_data Input data frame with dominance interactions, should only contain Date, 
#'   Winner, Loser.  Date should be formatted as MONTH/DAY/YEAR, or already as Date class.
#' @param pres_data Input data frame with columns "id", "start_date" and "end_date".  Date
#'   columns should be formatted as MONTH/DAY/YEAR, or already as Date class.  If all IDs 
#'   are present the whole time, you can ignore this and a pres_data table will be 
#'   automatically generated.
#' @param fit_init_elo If FALSE (the default), fits only the K parameter, with a default 
#'   starting Elo score of 1000 for each individual.  If TRUE, fits K and starting Elo for 
#'   each individual.  The latter option is \emph{much} slower.
#' @param outputfile Name of csv file to save ranks to.  Default is NULL, in which case 
#'   the function will only return a table in R.  If you supply an output file name
#'   the function will save the results as a csv file in your working directory.
#' @param returnR whether to return an R object from the function call.  Default is TRUE
#' 
#' @details This function accepts a data frame of date-stamped dominance interactions and 
#'   (optionally) a data frame of start and end dates for each individual to be ranked, 
#'   and outputs daily Elo scores with K parameter, and optionally initial elo scores, fitted using 
#'   a maximum likelihood approach.  The optimization procedure uses the \code{optim()} function, 
#'   with a burn in period of 100 interactions.  We use the "Brent" method when fitting only the K 
#'   parameter, and the "BFGS" method for fitting both K and initial Elo scores.  See 
#'   \code{\link[stats]{optim}} for more details.  Future package development will add additional 
#'   user control of the optimization procedure, allowing for specification of the burn in period, 
#'   optimization algorithm, and initial values for optimization.  
#'   
#'   Note also that the fitting procedure requires each individual to have at least one win and 
#'   one loss, so any individual that doesn't meet those criteria is automatically removed.  
#'   Additionally, any instance of an individual winning against itself is cleaned from the data,
#'   and several other checks of the data are performed before the optimization procedure is run.
#'   
#'   A detailed description of the function output is given in the \strong{Value} section of 
#'   this help file:
#'   
#' 
#' @return Returns a list with five or six elements (depending on input): 
#' \itemize{
#'  \item{\strong{elo}}{ Data frame with all IDs and dates they were present, with the following columns:}
#'    \itemize{
#'      \item{Date}{: Dates of study period}
#'      \item{Individual}{: the names of each ranked individual, for each date they were present}
#'      \item{Elo}{: fitted Elo scores for each individual on each day}
#'      \item{EloOrdinal}{: Daily ordinal rank based on Elo scores}
#'      \item{EloScaled}{: Daily Elo scores rescaled between 0 and 1 according to 
#'        \deqn{([individual Elo] - min([daily Elo scores])/(max([daily Elo scores]) - min([daily Elo scores]))}}
#'      \item{ExpNumBeaten}{: expected number of individuals in the group beaten, which is the sum of 
#'        winning probabilities based on relative Elo scores of an individual and all others, following 
#'        equation (4) in Foerster, Franz et al. 2016}
#'      \item{EloCardinal}{: ExpNumBeaten values rescaled as a percentage of the total number of ranked 
#'        individuals present in the group on the day of ranking. We encourage the use of this measure.}
#'      \item{JenksEloCardinal}{: Categorical rank (high, mid, or low) using the Jenks natural breaks 
#'        classification method implemented in the R package BAMMtools. 
#'        See \code{\link[BAMMtools]{getJenksBreaks}}}
#'      }
#'  \item{\strong{k}}{ The maximum-likelihood fitted k parameter value}
#'  \item{\strong{pred_accuracy}}{ Proportion of correctly predicted interactions}
#'  \item{\strong{maxLogL}}{ The overall log-likelihood of the observed data given the fitted parameter values 
#'    based on winning probabilities (as calculated in equation (1) of Foerster, Franz et al 2016) for all 
#'    interactions}
#'  \item{\strong{AIC}}{ Akaike's Information Criterion value as a measure of model fit}
#'  \item{\strong{init_elo}}{ (\emph{Only returned if you fit initial Elo scores}) initial Elo for each individual}
#'  } 
#' 
#' @examples
#'
#' nbadata = EloOptimized::nba #nba wins and losses from the 1995-96 season
#' nbaelo = eloratingopt(agon_data = nbadata, fit_init_elo = FALSE)
#' # generates optimized elo scores (optimizing only K) and saves them as "nbaelo" 
#' 
#' @export
#' @importFrom stats approx ave optim reshape
#' @importFrom utils read.csv write.csv
#' @importFrom rlang .data
#' @import reshape2
#' @import BAMMtools
#' @importFrom magrittr "%>%"




eloratingopt <- function(agon_data, pres_data, fit_init_elo = FALSE, outputfile = NULL, returnR = TRUE){
  # Get data
  
  if(length(outputfile) == 0 & returnR == FALSE){
    stop("supply an outputfile name or set returnR to TRUE (or both)")
  }
  
  ago = agon_data
  names(ago) <- tools::toTitleCase(tolower(names(ago)))
  if(!all(names(ago) %in% c("Date", "Winner", "Loser"))){
    stop("colnames in agonistic data should be 'Date', 'Winner', 'Loser' (not case sensitive)")
  }
  
  # to clean up data from readr::read_csv()
  attr(ago, "spec") = NULL
  ago = as.data.frame(ago)
  
  ago$Winner = as.character(ago$Winner)
  ago$Loser = as.character(ago$Loser)
  
  if(any(ago$Winner == ago$Loser)){
    stop("can't have same ID win and lose in one interaction")
  }
  
  if(class(ago$Date) != "Date"){
    ago$Date = lubridate::mdy(ago$Date)
  }
  
  if(any(ago$Date < dplyr::lag(x = ago$Date, 
                               n = 1, 
                               # to avoid NA in first value
                               default = min(ago$Date) - 
                                lubridate::years(1)))){
    stop("agon_data dates should be in chronological order")
  }
  
  
  if(missing(pres_data)){
    
    startids = sort(unique(c(ago$Winner, ago$Loser)))
    presence = data.frame(id = startids,
                          start_date = lubridate::ymd(sapply(startids, function(x) 
                            as.character(min(ago[ago$Winner == x | 
                                                   ago$Loser == x, "Date"])))),                                                                               
                          end_date = max(ago$Date), 
                          stringsAsFactors = F) 
    rm(startids)
    
  } else {
    
    presence <- pres_data
    # clean up data from readr::read_csv()
    attr(presence, "spec") = NULL
    presence = as.data.frame(presence)
    names(presence) = tolower(names(presence))
    if(!all(names(presence) %in% c("id", "start_date", "end_date"))){
      stop("colnames in presence data should be 'id', 'start_date', 'end_date' (not case sensitive)")
    }
    if(class(presence$start_date) != "Date"){
      presence$start_date = lubridate::mdy(presence$start_date)}
    if(class(presence$end_date) != "Date"){
      presence$end_date = lubridate::mdy(presence$end_date)}
    
    presence$id = as.character(presence$id)
    
  }
  
  if(fit_init_elo == FALSE){mod_type = 1} else {mod_type = 3}
  
  
  # ---------------  Make sure ago and presence have same start and end dates ----------------------
  # **** maybe we don't need to do this anymore??? 
  # or maybe move this after the filtering in case the ago file is altered much
  
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
    
    presence = presence %>% dplyr::filter(.data$wl != 0)
    
    ago = ago %>% 
      dplyr::filter(.data$Winner %in% presence$id & 
                      .data$Loser %in% presence$id)
    
    if(nrow(presence) == oldnum) break
    
  }
  
  presence = presence[,-4] # remove dummy variable
  
  all_inds = sort(presence$id)
  
  
  if(mod_type == 1 & nrow(ago) <= 100){
    stop("Currently you can't fit only K with less than 100 interactions (after filtering) due to burn in")
  }
  
  # error in case all interactions fall outside presence window:
  if(any(apply(presence, MARGIN = 1, function(x){
      
      sum(ago$Date >= x[2] & ago$Date <= x[3] & (ago$Winner == x[1] | ago$Loser == x[1]))
      
    }) == 0)){
    
    bad = sum((apply(presence, MARGIN = 1, function(x){
      
      sum(ago$Date >= x[2] & ago$Date <= x[3] & (ago$Winner == x[1] | ago$Loser == x[1]))
      
    })) == 0)
    
    stop(paste(bad, "individual(s) have no interactions within their presence window after filtering"))
    
  }
  
  
  # ---------------   Fit models  --------------------------------
  
  
  if(mod_type == 1){
    # Model 1 (for males)
    model <- optim(par=5, burn_in=100, elo.model1, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='Brent', lower=-10, upper=10)
    model_log <- elo.model1(par=model$par, burn_in=100, all_ids = all_inds, IA_data = ago, return_likelihood=F)
    # model <<- res_m_model1
    # model_log <<- res_m_model1_log
    pred_accuracy <- mean(model_log$elo_w_before[101:nrow(model_log)] > model_log$elo_l_before[101:nrow(model_log)])
  } else if(mod_type == 3) {
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
  if(mod_type == 3){
    
    # Get elo from log object
    initelo <- data.frame(Date = presence[order(presence$id), "start_date"],
                          Individual = all_inds,
                          EloScoreAfter = model$par[2: length(model$par)], 
                          stringsAsFactors = F)
    # names(elo) <- all_inds #pretty sure this should be "all_inds", but DOUBLE CHECK!!! (
    # changed from "inds" to "all_inds")
    
    df <- model_log[, names(model_log) %in% c("Date", "Winner", "Loser", "elo_w_after", "elo_l_after")]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Format columns
    df2 <- seq_long[,c(1,3,4)]
    row.names(df2) = NULL
    
    df2 = rbind.data.frame(initelo, df2) #combine starting elo scores with elo scores after interactions
    
    row.names(df2) = NULL
    
  } else if(mod_type == 1) {
    # Reformat elo-after scores of winners and losers into long format
    df <- model_log[, names(model_log) %in% c("Date", "Winner", "Loser", "elo_w_after", "elo_l_after")] #model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Format columns
    df2 <- seq_long[,c(1,3,4)] # make this a dplyr::select, then dplyr::arrange 
    
    #skipping step with the male models to add starting elo scores. Thus in male models 
    # individuals start being ranked after their first interaction, whereas
    # in female models individuals start being ranked immediately upon entry.
    
  }
  
  # Order by date and ID
  df2 <- df2[order(df2$Date, df2$Individual),]  # maybe can remove this, appears default behavior is to reorder...
  
  # Use max achieved score per day
  df2_daymax <-
    df2 %>%
    dplyr::group_by(.data$Date, .data$Individual) %>%
    dplyr::summarise(EloScoreAfterMax = max(.data$EloScoreAfter)) %>%
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
  presence_long$Individual = as.character(presence_long$Individual)
  
  # 
  # presence_long$Elo = df2_daymax$EloScoreAfterMax[match(paste0(presence_long$Individual, presence_long$Date),
  #                                                       paste0(df2_daymax$Individual, df2_daymax$Date))]
  
  # add elo scores to presence data and interpolate:
  
  presence_long = dplyr::left_join(x = presence_long, y = df2_daymax, 
                                    by = c("Individual" = "Individual", "Date" = "Date")) %>% 
    dplyr::rename(Elo = .data$EloScoreAfterMax) %>%
    dplyr::group_by(.data$Individual) %>%
    dplyr::mutate(Elo_interpol = approx(.data$Elo, xout = 1:length(.data$Elo), 
                                        rule = 1:2, method = "constant")$y) %>%
    as.data.frame()
  
  # post-processing:
  
  elo_long = 
    dplyr::filter(presence_long, !is.na(.data$Elo_interpol)) %>% 
    dplyr::select(-.data$Elo) %>% 
    dplyr::rename(EloScore = .data$Elo_interpol) %>%
    dplyr::group_by(.data$Date) %>%
    dplyr::mutate(EloNorm = (.data$EloScore - min(.data$EloScore, na.rm = T))/
                    (max(.data$EloScore, na.rm = T) - min(.data$EloScore, na.rm = T)),
                  rank_ord = dplyr::row_number(dplyr::desc(.data$EloScore)),
                  pct_beaten = cardinalize(.data$EloScore),
                  elo_rel = relativize(.data$pct_beaten),
                  JenksEloCardinal = jenksify(.data$elo_rel)) %>%
    dplyr::arrange(.data$Date, .data$Individual) %>%
    dplyr::select(.data$Date, .data$Individual, Elo = .data$EloScore, 
                  EloOrdinal = .data$rank_ord, EloScaled = .data$EloNorm, 
                  ExpNumBeaten = .data$pct_beaten, EloCardinal = .data$elo_rel, 
                  JenksEloCardinal = .data$JenksEloCardinal) %>%
    as.data.frame()
  
  
    # as.data.frame()
    # 
    # colnames(elo_long) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "ExpNumBeaten", "EloCardinal", "JenksEloCardinal")
    # 
    # 
  # ----------------- Step 2: Ordinal ranks by day -----------------------------------
  
  # elo_long$rank_ord <- ave(elo_long$EloScore, as.character(elo_long$Date), FUN = function(x) rank(-x, ties.method = "first"))
  
  # elo_long <- # don't forget to add the .data$ stuff in, but now we can combine this with the below steps!
  #   elo_long %>% 
  #   dplyr::group_by(.data$Date) %>% 
  #   dplyr::mutate(rank_ord = dplyr::row_number(dplyr::desc(.data$EloScore))) %>%
  #   as.data.frame()
  
  # ----------------- Step 3: Calculate cardinal ranks -------------------------------
  
  
  # cardinalize = function(x){
  #   carddat = sapply(x, function(y){
  #     sum(1 / (1 + exp(-0.01*(y - x))), na.rm=T) - .5 #subtracting .5 is equivalent to removing the prob of winning against oneself
  #     #b/c 1/(1 + exp(-0.01*0)) = 1/(1 + exp(0)) = 1/(1 + 1) = 1/2
  #   })
  #   return(carddat)
  # }
  # 
  # relativize = function(x){
  #   reldat = sapply(x, function(y){
  #     y/(length(x) - 1)
  #   })
  #   return(reldat)
  # }
  
  # elo_long =
  #   elo_long %>%
  #   dplyr::group_by(.data$Date) %>%
  #   dplyr::mutate(rank_ord = dplyr::row_number(dplyr::desc(.data$EloScore)),
  #                 pct_beaten = cardinalize(.data$EloScore),
  #                 elo_rel = relativize(.data$pct_beaten),
  #                 JenksEloCardinal = jenksify(.data$elo_rel)) %>% # later combine this with the ordinal rank step!
  #   as.data.frame()
  
  
  # --------------------- Step 4: Find natural breaks in list of elo scores by day --------------------
  
  # jenksify = function(x){
  #   breaks = BAMMtools::getJenksBreaks(x, 4)
  #   cats = ifelse(x <= breaks[[2]], "low",
  #                 ifelse(x > breaks[[3]], "high", "mid"))
  #   return(cats)
  # }
  
  # elo_long =
  #   elo_long %>%
  #   dplyr::group_by(.data$Date) %>%
  #   dplyr::mutate(JenksEloCardinal = jenksify(.data$elo_rel)) %>%
  #   as.data.frame()
  
  
  # ------------------------ Step 5: prettify -----------------------------------
  # elo_long =
  #   elo_long %>%
  #   dplyr::select(.data$Date, .data$Individual, .data$EloScore, .data$rank_ord, .data$EloNorm, 
  #                 .data$pct_beaten, .data$elo_rel, .data$JenksEloCardinal) %>%
  #   as.data.frame()
  # 
  # colnames(elo_long) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "ExpNumBeaten", "EloCardinal", "JenksEloCardinal")
  # 
  # k = exp(model$par[1])
  # pred_accuracy
  # AIC = 2*as.numeric(model$value) + 2*length(model$par)
  
  cat(paste0("k = ", round(exp(model$par[1]), 3), "\n"))
  cat(paste0("prediction accuracy = ", round(pred_accuracy, 3), "\n"))
  
  if(length(outputfile) > 0){
    write.csv(elo_long, outputfile, row.names = F)
  } 
  
  if(returnR == TRUE){
    
    res = list()
    res$elo = elo_long
    res$k = exp(model$par[1])
    res$pred_accuracy = pred_accuracy
    res$maxLogL = unname(-model$value)
    res$AIC = 2*as.numeric(model$value) + 2*length(model$par)
    
    if(mod_type == 3){
      
      temp = elo_long[!duplicated(elo_long$Individual),]
      row.names(temp) = NULL
      
      res$init_elo = temp
      
      rm(temp)
      
    }
    
    return(res)
  }
  
}