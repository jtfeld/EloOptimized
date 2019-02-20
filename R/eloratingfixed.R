#' @title Create daily elo ranks and multiple derivatives with user-defined parameter values
#' @description Conducts traditional elo rating analyses using specified K value
#'   and outputs raw, normalized, cardinal, and  categorical ranks as a list object in 
#'   R or in an output file.  For optimized Elo parameters, use \code{\link{eloratingopt}}.
#' @usage eloratingfixed(agon_data, pres_data, k = 100, init_elo = 1000, outputfile = NULL, 
#'   returnR = TRUE, p_function = "sigmoid")
#' @param agon_data Input data frame with dominance interactions, should only contain Date, 
#'   Winner, Loser.  Date should be formatted as MONTH/DAY/YEAR, or already as Date class.
#' @param pres_data Input data frame with columns "id", "start_date" and "end_date".  Date
#'   columns should be formatted as MONTH/DAY/YEAR, or already as Date class.  If all IDs 
#'   are present the whole time, you ignore this and a pres_data table will be automatically
#'   generated.
#' @param k Specified value of the k parameter, default is 100
#' @param init_elo The starting Elo value for all individuals, default is 1000
#' @param outputfile Name of csv file to save ranks to.  Default is NULL, in which case 
#'   the function will only return a table in R.  If you supply an output file name
#'   the function will save the results as a csv file in your working directory.
#' @param returnR whether to return an R object from the function call.  Default is TRUE
#' @param p_function function defining probability of winning.  Default "sigmoid" is 
#'   equation (1) from Foerster, Franz et al 2016.  Use "pnorm" to use the 
#'   \code{\link[stats:Normal]{pnorm}}-based method implemented in the EloRating package. 
#'   
#' @details This function accepts a data frame of date-stamped dominance interactions and 
#'   (optionally) a data frame of start and end dates for each individual to be ranked, 
#'   and outputs daily Elo scores with parameters specified by the user.  The default function 
#'   used to determine probability of winning is equation (1) from Foerster, Franz et al. 2016, 
#'   but for ease of comparison with the EloRating package, we also added the option to use
#'   the \code{\link[stats:Normal]{pnorm}}-based method implemented in the EloRating package, and future 
#'   development will add the option to use the original function from Elo 1978 (as implemented in 
#'   the elo package).  This function does not require large presence matrices, and efficiently 
#'   calculates a series of additional indices (described below).  
#'   
#'   As opposed to the \code{\link{eloratingopt}} function, this procedure only requires that 
#'   included individuals have at least one win \emph{or} one loss. 
#'   
#'   A detailed description of the function output is given in the \strong{Value} section of 
#'   this help file:
#'   
#' @return Returns a list with six elements: 
#' \itemize{
#'  \item{\strong{elo}}{ Data frame with all IDs and dates they were present, with the following columns:}
#'    \itemize{
#'      \item{Date}{: Dates of study period}
#'      \item{Individual}{: the names of each ranked individual, for each date they were present}
#'      \item{Elo}{: fitted Elo scores for each individual on each day}
#'      \item{EloOrdinal}{: Daily ordinal rank based on Elo scores}
#'      \item{EloScaled}{: Daily Elo scores rescaled between 0 and 1 according to 
#'        \code{([individual Elo] - min([daily Elo scores])/(max([daily Elo scores]) - min([daily Elo scores]))}}
#'      \item{ExpNumBeaten}{: expected number of individuals in the group beaten, which is the sum of 
#'        winning probabilities based on relative Elo scores of an individual and all others, following 
#'        equation (4) in Foerster, Franz et al. 2016}
#'      \item{EloCardinal}{: ExpNumBeaten values rescaled as a percentage of the total number of ranked 
#'        individuals present in the group on the day of ranking. We encourage the use of this measure.}
#'      \item{JenksEloCardinal}{: Categorical rank (high, mid, or low) using the Jenks natural breaks 
#'        classification method implemented in the R package BAMMtools. 
#'        See \code{\link[BAMMtools]{getJenksBreaks}}}
#'      }
#'  \item{\strong{k}}{ User-defined value of the k parameter}
#'  \item{\strong{init_elo}}{ User-defined initial Elo score when individuals enter the hierarchy}
#'  \item{\strong{pred_accuracy}}{ Proportion of correctly predicted interactions}
#'  \item{\strong{logL}}{ The overall log-likelihood of the observed data given the user-supplied parameter 
#'    values based on winning probabilities (as calculated in equation (1) of Foerster, Franz et al 2016) 
#'    for all interactions}
#'  } 
#'   
#' @examples
#'
#' nbadata = EloOptimized::nba #nba wins and losses from the 1995-96 season
#' nbaelo = eloratingfixed(agon_data = nbadata)
#' # generates traditional Elo scores (with init_elo = 1000 & k = 100) and saves 
#' #   them as "nbaelo" 
#' 
#' @export
#' @importFrom stats approx ave optim reshape
#' @importFrom utils read.csv write.csv
#' @importFrom rlang .data
#' @import reshape2
#' @import BAMMtools
#' @importFrom magrittr "%>%"




eloratingfixed <- function(agon_data, pres_data, k = 100, init_elo = 1000, 
                          outputfile = NULL, returnR = TRUE, p_function = "sigmoid"){
  # Get data
  
  if(length(outputfile) == 0 & returnR == FALSE){
    stop("supply an outputfile name or set returnR to TRUE (or both)")
  }
  
  ago = agon_data
  names(ago) <- tools::toTitleCase(tolower(names(ago)))
  if(!all(names(ago) %in% c("Date", "Winner", "Loser"))){
    stop("colnames in agonistic data should be 'Date', 'Winner', 'Loser' (not case sensitive)")
  }
  
  ago$Winner = as.character(ago$Winner)
  ago$Loser = as.character(ago$Loser)
  
  if(any(ago$Winner == ago$Loser)){
    stop("can't have same ID win and lose in one interaction")
  }
  
  if(class(ago$Date) != "Date"){
    ago$Date = lubridate::mdy(ago$Date)
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
    names(presence) = tolower(names(presence))
    if(!all(names(presence) %in% c("id", "start_date", "end_date"))){
      stop("colnames in presence data should be 'id', 'start_date', 'end_date' (not case sensitive)")
    }
    if(class(pres_data$start_date) != "Date"){
      pres_data$start_date = lubridate::mdy(pres_data$start_date)}
    if(class(pres_data$end_date) != "Date"){
      pres_data$end_date = lubridate::mdy(pres_data$end_date)}
    
    presence$id = as.character(presence$id)
    
  }
  
  
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
  
  
  # ---------- Filter individuals who do not have at least one win or one loss ----------------
  
  presence$wl = 0 #add dummy column to count wins and losses
  
  # vectorized loop to remove individuals from presence and ago data with 0 wins AND 0 losses
  # this should work fine but 
  # repeat{
  #   
  #   oldnum = nrow(presence)
    
    presence$wl = sapply(X = presence$id, function(x) sum(ago$Winner == x) + sum(ago$Loser == x))
    
    presence = presence %>% dplyr::filter(.data$wl != 0) %>% dplyr::select(-.data$wl)
    
  #   ago = ago %>% 
  #     dplyr::filter(.data$Winner %in% presence$id & 
  #                     .data$Loser %in% presence$id)
  #   
  #   if(nrow(presence) == oldnum) break
  #   
  # }
  
  # presence = presence[,-4] # remove dummy variable
  
  all_inds = sort(presence$id)
  
  
  # ---------------   Fit models  --------------------------------
  
  
  # if(mod_type == 1){
  #   # Model 1 (for males)
  #   model <- optim(par=5, burn_in=100, elo.model1, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='Brent', lower=-10, upper=10)
  #   model_log <- elo.model1(par=model$par, burn_in=100, all_ids = all_inds, IA_data = ago, return_likelihood=F)
  #   # model <<- res_m_model1
  #   # model_log <<- res_m_model1_log
  #   pred_accuracy <- mean(model_log$elo_w_before[101:nrow(model_log)] > model_log$elo_l_before[101:nrow(model_log)])
  # } else if(mod_type == 3) {
  #   # Model 3 (for females)
  #   # model <- optim(par=c(5, rep(0, length(all_inds))), elo.model3, all_ids = all_inds, IA_data = ago, return_likelihood=T, method='BFGS', control = list(maxit = 10000, reltol=1e-10))
  #   model <- optim(par=c(5, rep(0, length(all_inds))), elo.m3_lik_vect, all_ids = all_inds, IA_data = ago, method='BFGS', control = list(maxit = 10000, reltol=1e-10))
  #   ### USE SAVED "../data prep code/fem_mod_kk_2013.RData" TO SAVE TIME!
  #   model_log <- elo.model3(par=model$par, all_ids = all_inds, IA_data = ago, return_likelihood=F)
  #   # model <- res_fem_model3
  #   # model_log <- res_fem_model3_log
  #   pred_accuracy <- mean(model_log$elo_w_before > model_log$elo_l_before)
  # }
  
  model = elo.model1(par = log(k), burn_in = 0, init_elo = init_elo, IA_data = ago, all_ids = all_inds, 
                     return_likelihood = T, p_function = p_function)
  model_log = elo.model1(par = log(k), burn_in = 0, init_elo = init_elo, IA_data = ago, all_ids = all_inds, 
                         return_likelihood = F, p_function = p_function)
  pred_accuracy <- mean(model_log$elo_w_before > model_log$elo_l_before)
  
  # ================ Post-processing of elo scores =================================
  
  # ---------------- Step 1:  Normalize elo scores by date -------------------------
  
  # For females, start here
  # if(mod_type == 3){
  #   
  #   # Get elo from log object
  #   initelo <- data.frame(Date = presence[order(presence$id), "start_date"],
  #                         Individual = all_inds,
  #                         EloScoreAfter = model$par[2: length(model$par)], 
  #                         stringsAsFactors = F)
  #   # names(elo) <- all_inds #pretty sure this should be "all_inds", but DOUBLE CHECK!!! (
  #   # changed from "inds" to "all_inds")
  #   
  #   df <- model_log[, names(model_log) %in% c("Date", "Winner", "Loser", "elo_w_after", "elo_l_after")]
  #   seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
  #   # Format columns
  #   df2 <- seq_long[,c(1,3,4)]
  #   row.names(df2) = NULL
  #   
  #   df2 = rbind.data.frame(initelo, df2) #combine starting elo scores with elo scores after interactions
  #   
  #   row.names(df2) = NULL
  #   
  # } else if(mod_type == 1) {
    # Reformat elo-after scores of winners and losers into long format
    df <- model_log[, names(model_log) %in% c("Date", "Winner", "Loser", "elo_w_after", "elo_l_after")] #model_log[, c(1:3, 6:7)]
    seq_long <- reshape(df, varying=list(c(2:3), c(4:5)), v.names=c("Individual", "EloScoreAfter"), direction="long")
    # Format columns
    df2 <- seq_long[,c(1,3,4)]
    
    #skipping step with the male models to add starting elo scores. Thus in male models 
    # individuals start being ranked after their first interaction, whereas
    # in female models individuals start being ranked immediately upon entry.
    
  # }
  
  # Order by date and ID
  df2 <- df2[order(df2$Date, df2$Individual),]  
  
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
  
  # elo_data2 = dplyr::filter(presence_long, !is.na(.data$Elo_interpol)) %>% 
  #   dplyr::select(-.data$Elo) %>% 
  #   dplyr::rename(EloScore = .data$Elo_interpol) %>%
  #   as.data.frame()
  # 
  # elo_long =
  #   elo_data2 %>%
  #   dplyr::group_by(.data$Date) %>%
  #   dplyr::mutate(EloNorm = (.data$EloScore - min(.data$EloScore, na.rm = T))/
  #                   (max(.data$EloScore, na.rm = T) - min(.data$EloScore, na.rm = T))) %>%
  #   dplyr::arrange(.data$Date, .data$Individual) %>%
  #   as.data.frame()
  # 
  # 
  # 
  # # ----------------- Step 2: Ordinal ranks by day -----------------------------------
  # 
  # elo_long$rank_ord <- ave(elo_long$EloScore, as.character(elo_long$Date), FUN = function(x) rank(-x, ties.method = "first"))
  # 
  # # ----------------- Step 3: Calculate cardinal ranks -------------------------------
  # 
  # 
  # cardinalize = function(x){
  #   carddat = c()
  #   carddat = sapply(x, function(y){
  #     sum(1 / (1 + exp(-0.01*(y - x))), na.rm=T) - .5 #subtracting .5 is equivalent to removing the prob of winning against oneself
  #     #b/c 1/(1 + exp(-0.01*0)) = 1/(1 + exp(0)) = 1/(1 + 1) = 1/2
  #   })
  #   return(carddat)
  # }
  # 
  # relativize = function(x){
  #   reldat = c()
  #   reldat = sapply(x, function(y){
  #     y/(length(x) - 1)
  #   })
  #   return(reldat)
  # }
  # 
  # elo_long =
  #   elo_long %>%
  #   dplyr::group_by(.data$Date) %>%
  #   dplyr::mutate(pct_beaten = cardinalize(.data$EloScore),
  #                 elo_rel = relativize(.data$pct_beaten)) %>%
  #   as.data.frame()
  # 
  # 
  # # --------------------- Step 4: Find natural breaks in list of elo scores by day --------------------
  # 
  # jenksify = function(x){
  #   breaks = BAMMtools::getJenksBreaks(x, 4)
  #   # cats = c()
  #   cats = ifelse(x <= breaks[[2]], "low",
  #                 ifelse(x > breaks[[3]], "high", "mid"))
  #   return(cats)
  # }
  # 
  # elo_long =
  #   elo_long %>%
  #   dplyr::group_by(.data$Date) %>%
  #   dplyr::mutate(JenksEloCardinal = jenksify(.data$elo_rel)) %>%
  #   as.data.frame()
  # 
  # 
  # # ------------------------ Step 5: prettify -----------------------------------
  # elo_long =
  #   elo_long %>%
  #   dplyr::select(.data$Date, .data$Individual, .data$EloScore, .data$rank_ord, .data$EloNorm, 
  #                 .data$pct_beaten, .data$elo_rel, .data$JenksEloCardinal) %>%
  #   as.data.frame()
  # 
  # colnames(elo_long) <- c("Date", "Individual", "Elo", "EloOrdinal", "EloScaled", "ExpNumBeaten", "EloCardinal", "JenksEloCardinal")
  # 
  cat(paste0("k = ", k, "\n"))
  cat(paste0("prediction accuracy = ", round(pred_accuracy, 3), "\n"))
  
  if(length(outputfile) > 0){
    write.csv(elo_long, outputfile, row.names = F)
  } 
  
  if(returnR == TRUE){
    
    res = list()
    res$elo = elo_long
    res$k = k
    res$init_elo = init_elo
    res$pred_accuracy = pred_accuracy
    res$logL = -as.numeric(model)
    # res$AIC = 2*as.numeric(model) + 0 # becausee no fitted parameters, AIC not appropriate...
    
    return(res)
  }
  
}