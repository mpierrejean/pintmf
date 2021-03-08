library("parallel")
library("PintMF")
library(future.apply)
source("inst/Simulations/00.setup.R")
listBenchmark <- list.files(pathDat)
nbCPU <- 2
ARI <- lapply(seq(along=listBenchmark)[1:8], function(ii){
  b <- listBenchmark[ii]
  K <- nclust[ii]
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  print("my_meth")
  true.clusters <- list.sim[[1]]$true.clust

  my_meth_results <- readRDS(file=file.path(pathMeth_sub, sprintf("my_meth_res.rds")))
  ARIs_my_meth <- sapply(my_meth_results, function (R){
    clust <- R$W %>% dist %>% hclust(method="ward.D2") %>% cutree(K)
    mclust::adjustedRandIndex(clust, true.clusters)

  })
  data.frame(ARI=ARIs_my_meth, Benchmark= b)
})


ARI_df <- do.call(rbind, ARI)


ARI_supervised <- lapply(seq(along=listBenchmark)[1:8], function(ii){
  b <- listBenchmark[ii]
  K <- nclust[ii]
  pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))
  print(pathMeth_sub)
  list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
  print("my_meth")
  true.clusters <- list.sim[[1]]$true.clust

  my_meth_results <- readRDS(file=file.path(pathMeth_sub, sprintf("my_meth_res_supervised.rds")))
  ARIs_my_meth <- sapply(my_meth_results, function (R){
    clust <- R$W %>% dist %>% hclust(method="ward.D2") %>% cutree(K)
    mclust::adjustedRandIndex(clust, true.clusters)

  })
  data.frame(ARI=ARIs_my_meth, Benchmark= b)
})

ARI_supervised_df <- do.call(rbind, ARI_supervised)


library("tidyverse")
methods <- c("CIMLR", "NMF","SGCCA", "MoCluster","iCluster")

true.clusters <- mapply(function(n, nnC) {
  rep(1:n, nnC)
}, nclust,n_by_Clust)

listBenchmark <- list.files(pathDat)


ARI_dat <- do.call(rbind, lapply(1:8, function(ii){
  b <- listBenchmark[ii]
  pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))

  adjustedRIs <- do.call(rbind, lapply(methods, function(mm){
    print(mm)
    ff <- file.path(pathMeth_sub, sprintf("%s_res.rds", mm))
    r <- readRDS(ff)
    adjRI <- sapply(r, function(ss) {
      if (BBmisc::is.error(ss)) {
        return(NA)
      }
      if (is.na(ss)||is.null(ss)) {
        return(NA)
      } else{
        if(mm=="coca"){
          ss  %>% mclust::adjustedRandIndex(true.clusters[[ii]])

        } else{
          ss$clust  %>% mclust::adjustedRandIndex(true.clusters[[ii]])
        }
      }
    })
    data.frame(method = mm, ARI = adjRI)
  }))
  adjustedRIs$Benchmark = b
  return(adjustedRIs)
}))

ARI_dat$method <- gsub("NMF", "intNMF", ARI_dat$method)
ARI_dat$method <- gsub("iCluster", "iClusterPlus", ARI_dat$method)

ARI_dat$method <- ARI_dat$method %>% factor(levels =
                                              c(
                                                "intNMF",
                                                "SGCCA",
                                                "MoCluster",
                                                "iClusterPlus",
                                                "CIMLR"

                                              ))
ARI_df <- ARI_df %>% mutate(method="PIntMF")
#ARI_supervised_df <- ARI_supervised_df %>% mutate(method="PIntMF supervised")

ARI_dat <- rbind(ARI_dat,ARI_df)

g_adj <- ARI_dat %>%
  mutate(method=factor(method,  levels=c("PIntMF", "intNMF",
                                         "SGCCA", "MoCluster",
                                         "iClusterPlus", "CIMLR")))%>%
  ggplot(aes(x = Benchmark, y = ARI)) +
  geom_boxplot(fill="#696969") + ylab("Adjusted Rand Index") +
  theme_bw() +
facet_wrap(.~method, ncol=2,
           labeller = label_wrap_gen(multi_line=TRUE))+
  theme(axis.text.x = element_text(vjust = 0.5, size = 10,
                                   angle = 90,
                                   ), strip.text.x = element_text(size = 10), legend.position = "none", axis.title= element_blank(), axis.text.y = element_text(size=12))

g_adj

# define your own tagss
ggsave(filename=file.path("../Figs/", "Clust_eval_benchmark1-4.pdf"),
       g_adj, width=7, height=7)


ARI_dat <- rbind(ARI_df,ARI_supervised_df)
g_adj <- ARI_dat %>% ggplot(aes(x = method, y = ARI, fill = method)) +
  geom_boxplot() + ylab("Adjusted Rand Index") + theme_bw() +
  scale_fill_brewer(palette = "Set2")+facet_wrap(.~Benchmark, ncol=4, labeller = label_wrap_gen(multi_line=TRUE))+
  theme(axis.text.x = element_text(size = 12, angle = 90, color=c("#993333","#993333" )), strip.text.x = element_text(size = 15), legend.position = "none", axis.title= element_blank())+scale_x_discrete(labels=c("PIntMF"=expression(bold(PIntMF)), parse=TRUE))

g_adj
