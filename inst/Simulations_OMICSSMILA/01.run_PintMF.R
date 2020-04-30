library("parallel")
library(PintMF)
library(future.apply)
source("inst/Simulations_OMICSSMILA/00.setup.R")
for(ii in 1:n_batch){
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
  data_filter_t[["meth"]] <- log2((data_filter[["meth"]]+0.001)/(1-(data_filter[["meth"]]+0.001)))
  print("my_meth")
  my_meth_results <- data_filter %>% SolveInt(p=2, max.it=5, flavor_mod = "glmnet", init_flavor = "snf")
  clust <- my_meth_results$W %>% dist %>% hclust(method="ward.D2") %>% cutree(2)

    H_coef <- lapply(my_meth_results$H, function(hh) apply(hh, 2, max))
    H_sorted <- lapply(H_coef, function(h) sort(h,decreasing = TRUE) %>% names)

    t_meth <- intersect(file[[ii]][4][[1]], my_meth_results$H[[1]] %>% colnames())

  truth <- list(posMeth= t_meth, posRNA= file[[ii]][[2]] %>% data.frame %>% dplyr::select(contains("DE")) %>% colnames)
  J = sapply(data_filter, ncol)
  FPR = FPR_compute(truth, H_sorted[1:2], J[1:2])
  TPR = TPR_compute(truth, H_sorted[1:2])
  plot(FPR[[1]], TPR[[1]])
  plot(FPR[[2]], TPR[[2]])

  saveRDS(my_meth_results, file=file.path(pathMeth_sub, sprintf("my_meth_res.rds")))
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

