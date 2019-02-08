
#' @title internal fn to create cardinal rank scores
#' @description internal function for generating cardinal ranks
#' @usage cardinalize(x)
#' @param x input vector
#' @details converts raw Elo scores into predicted number of individuals beaten
#'   (using Equation 1 from paper)
#' @return returns new vector of cardinal rank scores
#' @details subtracting .5 is equivalent to removing the prob of winning against oneself
#'   because 1/(1 + exp(-0.01*0)) = 1/(1 + exp(0)) = 1/(1 + 1) = 1/2

cardinalize = function(x){
  carddat = sapply(x, function(y){
    sum(1 / (1 + exp(-0.01*(y - x))), na.rm=T) - .5 
  })
  return(carddat)
}

#' @title internal fn to relativize rank scores
#' @description internal function for generating scaled cardinal ranks
#' @usage relativize(x)
#' @param x input vector
#' @details scales cardinal Elo scores between 0 and 1
#' @return returns new vector of scaled rank scores

relativize = function(x){
  reldat = sapply(x, function(y){
    y/(length(x) - 1)
  })
  return(reldat)
}

#' @title internal fn to generate categorical ranks
#' @description internal function for generating categorical ranks using 
#'   jenks natural breaks algorithm
#' @usage jenksify(x)
#' @param x input vector
#' @details creates categorical ranks using jenks natural breaks algorithm
#' @return returns new vector of categorical ranks (high/medium/low)

jenksify = function(x){
  breaks = BAMMtools::getJenksBreaks(x, 4)
  cats = ifelse(x <= breaks[[2]], "low",
                ifelse(x > breaks[[3]], "high", "mid"))
  return(cats)
}