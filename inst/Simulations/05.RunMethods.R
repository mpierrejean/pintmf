library("parallel")
library(PintMF)
library(future.apply)
source("inst/Simulations/00.setup.R")
listBenchmark <- list.files(pathDat)
nbCPU <- 2
for(ii in 2:8){
  b <- listBenchmark[ii]
  K <- nclust[ii]
  pathDat_sim <- (sprintf("%s/%s", pathDat, b))
  pathMeth_sub <- (sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  data <- lapply(list.sim, function (ll) ll$data)
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
  data_filter <- data %>% lapply(remove_zero)
  data_t_filter <- lapply(data_filter, function(dd){
    dd[[3]] <- log2(dd[[3]]/(1-dd[[3]]))
    dd
  })
  print("my_meth")
  true.clusters <- list.sim[[1]]$true.clust

  my_meth_results <- lapply(data_t_filter, SolveInt, group=true.clusters, p=K, max.it=5, flavor_mod = "glmnet", init_flavor = "snf")
  saveRDS(my_meth_results, file=file.path(pathMeth_sub, sprintf("my_meth_res_supervised.rds")))
}
