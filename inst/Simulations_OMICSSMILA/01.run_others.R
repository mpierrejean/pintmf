# Run Mocluster 6/10/2020
library("parallel")
library(future.apply)
library(dplyr)
library(CrIMMix)
source("inst/Simulations_OMICSSMILA/00.setup.R")
print(n_batch)
for(ii in 1:n_batch){
  pathMeth_sub <- "../data/extdata_PintMF/res_v2/"
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

  #mo_results <- data_filter_t %>% IntMultiOmics(K=2,method="Mocluster", k=0.1)
  #saveRDS(mo_results, file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))
  sgcca_res <- data_filter_t %>% IntMultiOmics(K=2,method="SGCCA", ncomp = rep(2, 3))
  saveRDS(sgcca_res, file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))
 #intnmf_res <- data_filter %>% IntMultiOmics(K=2,method="intNMF")
 #saveRDS(intnmf_res, file=file.path(pathMeth_sub, sprintf("intNMF_res_batch%s.rds",ii)))

  #icluster_res <- data_filter_t %>% IntMultiOmics(K=2,method="iCluster")
  #saveRDS(icluster_res, file=file.path(pathMeth_sub, sprintf("icluster_res_batch%s.rds",ii)))

}


