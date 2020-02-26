#' get W in supervised way
#'
#' @param Zbar summary of various blocks of ceofficient  matrices
#' @param Ybar summary of various block of  observations
#' @param group group for each individual
#'
#' @return  Wc matrix W
#' @export
#' @importFrom stats lm
#' @importFrom dplyr %>%
#'
get.W.supervised <- function(Zbar, Ybar, group) {
  Wc <-  lapply(1:nrow(Ybar), FUN=solvew.group, Ybar=Ybar, Zbar=Zbar, group=group) %>% simplify2array %>% t
  return(Wc)
}

solvew.group <- function(ind, Ybar, Zbar, group){
  g <- group[ind]
  ybar <- Ybar[ind,]
  p <- ncol(Zbar)
  w <- numeric(p)
  w[g] <-  abs(cov(ybar, Zbar[,g])/var(ybar)) ## Fit linear regression
  return(round(w, 1))
}
