### TEst 20/11
## Remplacement du colSUms par le colSds



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

doROC_Moa <- function(truth, fit){
  K <- fit@data %>% length
  a <- fit@loading
  J <- sapply(fit@data, nrow)
  selectVars_1 <- which(a %>% rowSums !=0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  regexp <- "[[:digit:]]+"
  test <-  apply(a, 1, function(aa) abs(aa) %>% sum)
  idx <- which(test!=0)
  test_o <- sort(test[idx], decreasing = TRUE) %>% names
  test_o <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), test_o)
    gsub(sprintf("_dat%s", kk), "", test_o[idx])
  })
  denom_tp <- sapply(truth, length)

  TPR_list <- lapply(1:K, function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:K, function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}

doROC_SGCCA <- function(truth, fit){
  a <- fit$a
  J <- sapply(a, nrow)
  test <- lapply(a, function(aa)  apply(aa, 1, function(aa) sum(abs(aa))))
  denom_tp <- sapply(truth, length)

  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names
  })
  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}


doROC_intNMF <- function(truth, fit){
  a <- fit$H
  J <- sapply(a, ncol)
  test <- lapply(a, function(aa)   apply(aa, 2, sd))

  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names
  })
  denom_tp <- sapply(truth, length)

  TPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })
  return(list(TPR = TPR_list, FPR = FPR_list))
}


doROC_iCluster <- function(truth, fit, data_filter_t){
  a <- lapply(1:length(fit$beta), function(ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    names(rowsum) <- 1:length(rowsum)
    upper=quantile(rowsum,prob=0.85)
    sigfeatures=names(which(rowsum>upper))
  })
  # process string
  test_o <- lapply(1:length(fit$beta), function (ii){
    rowsum=rowSums(abs(fit$beta[[ii]]))
    names(rowsum) <- colnames(data_filter_t[[ii]])
    names(sort(rowsum, decreasing = TRUE)) })
  J <- sapply(fit$beta, nrow)

  denom_tp <- sapply(truth, length)
  TPR_list <-lapply(1:length(test_o), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      tpr <- (tt_o[t] %>% intersect(truth[[ii]]) %>% length)/denom_tp[ii]
    })
  })

  FPR_list <- lapply(1:length(a), function(ii){
    tt_o <- test_o[[ii]]
    sapply(1:length(tt_o), function (tt){
      t <- 1:tt
      fpr <- (tt_o[t]%>% setdiff(truth[[ii]]) %>% length)/(J[ii]-denom_tp[ii])
    })
  })

  return(list(TPR = TPR_list, FPR = FPR_list))
}

