#' SolveInt
#' @title SolveInt solve equation for H and W
#'
#' @param Y observations, list of blocks of data
#' @param p The number of latent profiles
#' @param max.it by default 20. Maximum iteration of the algorithm, else until convergence
#'   of matrix W
#' @param group default is NULL (unsupervised way) otherwise the W matrix is
#'   fixed and algorithm is run in a supervised way.
#' @param type type of data (by default all are "none"),
#'   can take values into \code{ c("none", "methylation", "mutation")},
#'   this argument allows to adapt modelisation by adding constraints.
#'   For exemple methylation data takes values between 0 and 1, then H block corresponding to
#'   methylation data takes values between 0 and 1.
#' @param verbose A logical value indicating whether to print extra information.
#'   Defaults to FALSE
#' @param flavor_mod glmnet or (oem not implemented yet)
#' @param init_flavor flavor to initialize W matrix (hclust, SNF, random, pca or svd)
#' @param ... other arguments
#'
#' @return A list that contains :
#'   \code{H} is a list of matrix (latent profiles for each block of data)
#'   \code{W} is a matrix (weight matrix).
#'   \code{loss} similarity of the clustering between two iterations
#'   \code{pve} Percentage of variance explained
#' @export
#' @importFrom glmnet glmnet
#' @importFrom dplyr bind_cols
#' @importFrom future.apply future_lapply
#' @importFrom future plan
#' @importFrom nnet class.ind
#' @importFrom pcaMethods pca
#' @import stats
SolveInt <- function(Y, p, max.it=20, flavor_mod="glmnet", group=NULL, type=rep("none", length(Y)), init_flavor="snf", verbose=FALSE,...
) {
  if (!is.list(Y)) {
    stop("Y is not a list")
  }
  family = c("mutation", "none", "methylation")
  if (!all(type %in% family)) {
    stop(sprintf("Type is not in the list, choose type in %s", paste(family, collapse = ", ")))
  }

  if(verbose) message(sprintf("Run algorithm for p=%s\n", p))

  if(verbose) message(sprintf("type of initialization : %s\n", init_flavor))
  gaussian_type <- which(type == "none")
  Yt <- lapply(1:length(Y),
               function(yy) {
                 if (yy %in% gaussian_type) {
                   if(verbose) message (sprintf("Center cols of Y\n"))
                   y <- sweep(Y[[yy]],MARGIN = 2,colMeans(Y[[yy]]) )
                   # y <- scale(Y[[yy]])

                 }else{
                   y <- Y[[yy]]
                 }
                 return(y)
               })
  Yt <- Y
  if(verbose) message(sprintf("method used is %s:",flavor_mod))

  ############################################################
  ############ Group
  ############################################################
  if(!is.null(group)){
    if(!is.numeric(group)){
      if(is.character(group)){
        warning(sprintf("group must be a numeric vector \n variable group has been transformed to a numeric variable"))
        group = as.numeric(as.factor(group))
      }else{
        stop(sprintf("group must be a numeric or character vector of length %s", nrow(Y[[1]])))
        }
    }
    if (nrow(Y[[1]])!=length(group)) {
      stop(sprintf("group must be a numeric vector of length %s", nrow(Y[[1]])))
    }
  }

  ############################################################
  ############ initialization
  ############################################################
  initH_clust <- function(yyt, hhc, pp) {
    cluster <- cutree(hhc, k = pp)
    ## Averaged profiles by cluster
    t(sapply(split(as.data.frame(yyt), f = cluster), FUN = colMeans))
  }
  initH_random <- function(yyt, pp) {
    idx <- sample(1:nrow(yyt),pp)
    yyt[idx,]
  }
  initH_pca <- function(yyt, p) {
    lapply(yyt, function (h) pca(h, nPcs = p)@loadings %>% t)
  }
  initH_svd <- function(yyt, p) {
    lapply(yyt, function (h) svd(h, nv=p)$v %>% t)
  }

  if(verbose) message (sprintf("Init Hc \n"))

  dd <- lapply(Yt, dist)
  hc <- lapply(dd, hclust, method = "ward.D")

  if(init_flavor=="hclust"){
    Hc <- mapply(initH_clust,Yt,hc, rep(p, length(Y)), SIMPLIFY = FALSE)
  }
  if(init_flavor=="svd"){
    Hc <- initH_svd(Yt, p )
  }
  if(init_flavor=="pca"){
    Hc <- initH_pca(Yt, p )
  }
  if(init_flavor=="random"){
    Hc <- mapply(initH_random,Yt, rep(p, length(Y)), SIMPLIFY = FALSE)
  }

  it <- 1
  loss <- numeric(0)
  pve <- numeric(0)
  Wc <- matrix(0, nrow(Y[[1]]), p)
  while (it <= max.it) {
    if(verbose) message("Solve W\n")
    if(is.null(group)){

      if(it==1 & init_flavor=="snf"){

        Wc <- class.ind(init_SNF(Yt, K=p)$clust)
        W.old <- Wc

      }else{
        Hc_norm <- lapply(Hc, function(hh) hh/sqrt(ncol(hh)))
        Zbar <- t(bind_cols(lapply(Hc_norm, data.frame))) %>% as.matrix()
        Ybar <- bind_cols(lapply(Yt, data.frame)) %>% as.matrix
        Wc <- get.W(Zbar=Zbar,
                    Ybar=Ybar)
        if(it==1){
          W.old <- Wc
        }
      }
    }else{
      if(it==1){
        Wc <- class.ind(group)
        W.old <- Wc
      }else{
        Hc_norm <- lapply(Hc, function(hh) hh/sqrt(ncol(hh)))
        Zbar <- t(bind_cols(lapply(Hc_norm, data.frame))) %>% as.matrix()
        Ybar <- bind_cols(lapply(Yt, data.frame)) %>% as.matrix
        Wc <- get.W.supervised(Zbar=Zbar,
                    Ybar=Ybar, group=group)
      }
    }

    loss[it ] <- Wc %>% dist %>% hclust(method="ward.D2") %>% cutree(p) %>%
      mclust::adjustedRandIndex(W.old %>% dist %>% hclust(method="ward.D2") %>% cutree(p))
    if(verbose) message(sprintf("Loss equal to %s", loss[it]))
    if(loss[it] >0.1){
      if(verbose) message ("Solve Hc\n")
      Hc <- future_lapply(1:length(Yt), function(yy){
        if(verbose){
          message(sprintf("Solve data number %s:",yy))
        }
        get.H(y = Yt[[yy]], W = Wc, flavor_mod=flavor_mod, verbose)}
      )

      pve[it] <- PVE(Yt, Hc, Wc)
      if(verbose) (sprintf("loss : %s\n", round(loss[it + 1],2)))
      it <- it + 1
      if(verbose) message(sprintf("iteration number %s", cat(it, "\n")))
    }else{
      message(sprintf("The clusterings are too different between two iterations: stop here"))
      Wc <- W.old
      it <- max.it+1
    }
  }
  ## Define names of variables
  var_names <- lapply(Y, colnames)
  Hc <- lapply(1:length(Y), function(ii) {
    h <- Hc[[ii]]
    colnames(h) = var_names[[ii]]
    return(h)
  })
  return(list(H = Hc, W = Wc,  loss = loss, pve = pve))
}
