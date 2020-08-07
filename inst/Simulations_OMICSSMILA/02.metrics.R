library("parallel")
library(PintMF)
library(future.apply)
library(dplyr)
source("inst/Simulations_OMICSSMILA/00.setup.R")


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

ROCValues <- list()

ARIValues <- numeric(n_batch)
for(ii in 1:n_batch){
  pathMeth_sub <- "inst/extdata/"
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
  data_filter_t[["meth"]] <- log2((data_filter[["meth"]])/(1-(data_filter[["meth"]])))%>% t %>% na.omit %>%t
  my_meth_results_2 <- readRDS(file=file.path(pathMeth_sub, sprintf("my_meth_res_batch%s_v2.rds",ii)))
  print(file.path(pathMeth_sub, sprintf("my_meth_res_batch%s_v2.rds",ii)))
  clust <- my_meth_results_2$W %>% dist %>% hclust(method="ward.D2") %>% cutree(2)
  true.clust = rownames(data[[1]]) %>% stringr::str_detect("CASE")+1
  ARI = mclust::adjustedRandIndex(clust, true.clust)
  ARIValues[ii] <- ARI
}

library(CrIMMix)
ARIValues_mo  <- ARIValues_sgcca <- numeric(n_batch)

for(ii in 1:n_batch){
  pathMeth_sub <- "inst/extdata/"
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
  mo_results <- data_filter_t %>% IntMultiOmics(K=2,method="Mocluster", k=0.1)

  saveRDS(mo_results, file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))

}


for(ii in 1:n_batch){
  pathMeth_sub <- "inst/extdata/"
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
  sgcca_res <- data_filter_t %>% IntMultiOmics(K=2,method="SGCCA")
  saveRDS(sgcca_res, file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))
}


for(ii in 1:n_batch){
  pathMeth_sub <- "inst/extdata/"
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

  mo_results <- readRDS(file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))
  sgcca_results <- readRDS(file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))


  true.clust = rownames(data[[1]]) %>% stringr::str_detect("CASE")+1
  ARIValues_mo[[ii]] = mclust::adjustedRandIndex(mo_results$clust, true.clust)
  ARIValues_sgcca[[ii]] = mclust::adjustedRandIndex(sgcca_results$clust, true.clust)
}


## ROC
ROCValues_sgcca <- list()
ROCValues_mo <- list()

for(ii in 1:n_batch){
  pathMeth_sub <- "inst/extdata/"
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

  mo_results <- readRDS(file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))
  sgcca_results <- readRDS(file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))

  t_meth <- file[[ii]]$methpos %>% intersect( colnames(data_filter_t[["meth"]]))
  truth <- list(posMeth= t_meth, posRNA= file[[ii]][[2]] %>% data.frame %>% dplyr::select(contains("DE")) %>% colnames)

  fit <- mo_results$fit
  K <- 2
  a <- fit@loading
  J <- sapply(fit@data, nrow)
  selectVars_1 <- which(a %>% rowSums !=0) %>% names
  selectVars <- lapply(1:K, function (kk){
    idx <- grep(sprintf("dat%s", kk), selectVars_1)
    gsub(sprintf("_dat%s", kk), "", selectVars_1[idx])
  })
  test <- rowSums(abs(a))
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
  roc <- list(FPR=FPR_list, TPR=TPR_list)
  ROCValues_mo[[ii]] <- ComputepAUC(roc)


  fit <- sgcca_results$fit

  a <- fit$a
  J <- sapply(a, nrow)
  selectVars <- lapply(a, function(aa) which(rowSums(aa) != 0) %>% names)
  test <- lapply(a, function(aa) rowSums(abs(aa)))
  regexp <- "[[:digit:]]+"
  denom_tp <- sapply(truth, length)

  test_o <- lapply(test, function (tt){
    idx <- which(tt!=0)
    sort(tt[idx], decreasing = TRUE) %>% names
  })
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

  roc_sgcca <- list(FPR=FPR_list, TPR=TPR_list)
  ROCValues_sgcca[[ii]] <- ComputepAUC(roc_sgcca)
}

saveRDS(ROCValues_sgcca, file.path(pathMeth_sub, "sGCCA_ROC.rds"))
saveRDS(ROCValues_mo, file.path(pathMeth_sub, "mocluster_ROC.rds"))

ROCValues<- lapply(list.files(pathMeth_sub, pattern ="my_meth_ROC", full.names=TRUE), readRDS)
ROCValues_sgcca <- readRDS(file.path(pathMeth_sub, "sGCCA_ROC.rds"))
ROCValues_mo <- readRDS(file.path(pathMeth_sub, "mocluster_ROC.rds"))

roc_Pint <- do.call(rbind, ROCValues) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="PIntMF")
roc_Moclust <- do.call(rbind, ROCValues_mo) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="MoCluster")
label_names <- c(
  `dataset 1` = "Methylation",
  `dataset 2` = "Gene expression"
)

roc_sgcca <- do.call(rbind, ROCValues_sgcca) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="SGCCA")

ROC_df <- rbind(roc_Pint, roc_Moclust, roc_sgcca)

library(ggplot2)
g <- ggplot(ROC_df, aes(x=method, y=AUC))+geom_violin()+facet_wrap(.~dataset,labeller = as_labeller(label_names))+theme_bw()+ylim(c(0,1))
g
ggsave(filename="../../papers/NAR-LaTeX-20170208/OMICSsimla_new_sim_roc.eps", height=5)


d <- ROC_df %>% filter(dataset=="dataset 1")

lm(AUC~method, data=d) %>% summary
d <- ROC_df %>% filter(dataset=="dataset 2")

lm(AUC~method, data=d) %>% summary

t1 <- ROC_df %>% filter(dataset=="dataset 1", method=="SGCCA") %>% dplyr::select(AUC) %>% pull
t2 <- ROC_df %>% filter(dataset=="dataset 1", method=="PIntMF") %>% dplyr::select(AUC) %>% pull
t3 <- ROC_df %>% filter(dataset=="dataset 1", method=="MoCluster") %>% dplyr::select(AUC) %>% pull

p_val_d1 <- c(t.test(t1, t2, paired=TRUE)$p.value, t.test(t1, t3, paired=TRUE)$p.value, t.test(t2, t3, paired=TRUE)$p.value)



t1 <- ROC_df %>% filter(dataset=="dataset 2", method=="SGCCA") %>% dplyr::select(AUC) %>% pull
t2 <- ROC_df %>% filter(dataset=="dataset 2", method=="PIntMF") %>% dplyr::select(AUC) %>% pull
t3 <- ROC_df %>% filter(dataset=="dataset 2", method=="MoCluster") %>% dplyr::select(AUC) %>% pull

p_val_d2 <- c(t.test(t1, t2, paired=TRUE)$p.value, t.test(t1, t3, paired=TRUE)$p.value, t.test(t2, t3, paired=TRUE)$p.value)

d_pval<- rbind(p_val_d1, p_val_d2)
colnames(d_pval) <- c("SGCCA vs PIntMF", "SGCCA vs MoCluster", "PIntMF vs MoCluster")

d_pval*6
