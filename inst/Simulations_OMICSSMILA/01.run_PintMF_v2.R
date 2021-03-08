library("parallel")
library(PintMF)
library(future.apply)
library(dplyr)
source("inst/Simulations_OMICSSMILA/00.setup_v2.R")


ComputepAUC <- function(roc) {
  TPR <- roc$TPR
  FPR <- roc$FPR
  aucs <- sapply(1:length(FPR), function (ii){
    x <- FPR[[ii]]
    y <- TPR[[ii]]
    y[is.infinite(y)] <- NaN
    auc <- sum(tis::lintegrate(c(x,1), c(y,1), xint=c(x,1)))
    denom <- sum(tis::lintegrate(c(0, 0, max(x)), c(0, 1, 1), xint=c(0, 0, max(x))))
    res <- auc
    res
  })
  names(aucs) <- c("dataset 1", "dataset 2")
  aucs
}


TPR_compute <- function(truth, selected_var,nvar=NULL){
  ndat <- length(truth)
  denom_tp <- sapply(truth, length)
  tp <- future_lapply(1:ndat, function(ii){
    if(is.null(nvar)){nvar= length(selected_var[[ii]])}
    ho <- selected_var[[ii]][1:nvar]
    sapply(1:length(ho), function (tt){
      t <- 1:tt
      tpr <- (ho[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })
  return(tp)
}

FPR_compute <- function(truth, selected_var, J,nvar=NULL){
  ndat <- length(truth)
  denom_tp <- sapply(truth, length)
  fp <- future_lapply(1:ndat, function(ii){
    if(is.null(nvar)){nvar= length(selected_var[[ii]])}
    ho <- selected_var[[ii]][1:nvar]
    sapply(1:length(ho), function (tt){
      t <- 1:tt
      fpr <- (ho[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(fp)
}

for(ii in 1:n_batch){
  pathMeth_sub <- R.utils::Arguments$getWritablePath("../data/extdata_PintMF/res_v2")
  data <- file[[ii]][1:3]
  remove_zero <- function (dat){
    lapply(dat, function(dd){
      idx <- which(colSums(dd)==0)
      if(length(idx)!=0){
        return(dd[, -idx])
      }else{
        return(dd)
      }
    })
  }
  data_filter <- data %>% remove_zero
  data_filter_t <- data_filter
  data_filter_t[["meth"]] <- log2((data_filter[["meth"]]+0.0001)/(1-(data_filter[["meth"]]+0.0001)))%>% t %>% na.omit %>%t
  print("my_meth")

  my_meth_results_2 <- data_filter_t %>% SolveInt(p=2, max.it=5, flavor_mod = "glmnet", init_flavor = "snf")
#  my_meth_results_2 <- data_filter_t[c(1,3)] %>% SolveInt(p=2, max.it=5, flavor_mod = "glmnet", init_flavor = "snf")

  H_coef <- lapply(my_meth_results_2$H, function(hh) apply(hh, 2, diff) %>% abs)
  H_sorted <- lapply(H_coef, function(h) sort(h,decreasing = TRUE) %>% names)
  t_meth <- file[[ii]]$methpos %>% intersect( colnames(data_filter_t[["meth"]]))
  truth <- list(posMeth= t_meth, posRNA= file[[ii]][[2]] %>% data.frame %>% dplyr::select(contains("DE")) %>% colnames)
  J = sapply(data_filter_t, ncol)
  FPR = FPR_compute(truth, H_sorted[1:2], J[1:2])
  TPR = TPR_compute(truth, H_sorted[1:2])
  roc_pint <- list(FPR=FPR, TPR=TPR)
  auc <- ComputepAUC(roc_pint)
  saveRDS(auc,file=file.path(pathMeth_sub, sprintf("my_meth_ROC_batch%s_v2.rds",ii)))

  print(auc)
  saveRDS(my_meth_results_2, file=file.path(pathMeth_sub, sprintf("my_meth_res_batch%s_v2.rds",ii)))
}
