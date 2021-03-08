library("parallel")
library(PintMF)
library(future.apply)
library(dplyr)
source("inst/Simulations_OMICSSMILA/00.setup.R")
source("inst/Simulations_OMICSSMILA/useful_function.R")

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

  mo_results <- readRDS(file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))
  sgcca_results <- readRDS(file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))


  true.clust = rownames(data[[1]]) %>% stringr::str_detect("CASE")+1
  ARIValues_mo[[ii]] = mclust::adjustedRandIndex(mo_results$clust, true.clust)
  ARIValues_sgcca[[ii]] = mclust::adjustedRandIndex(sgcca_results$clust, true.clust)
}


## ROC
ROCValues_sgcca <- list()
ROCValues_mo <- list()
ROCValues_intNMF <- list()
ROCValues_icluster <- list()
pathMeth_sub <- "../data/extdata_PintMF/res_v2"

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
  data_filter_t[["meth"]] <- log2((data_filter[["meth"]]+0.0001)/(1-(data_filter[["meth"]]+0.0001)))%>% t %>% na.omit %>%t

  mo_results <- readRDS(file=file.path(pathMeth_sub, sprintf("Mocluster_res_batch%s.rds",ii)))
  sgcca_results <- readRDS(file=file.path(pathMeth_sub, sprintf("SGCCA_res_batch%s.rds",ii)))
  intNMF_results <- readRDS(file=file.path(pathMeth_sub, sprintf("intNMF_res_batch%s.rds",ii)))
  icluster_results <- readRDS(file=file.path(pathMeth_sub, sprintf("icluster_res_batch%s.rds",ii)))

  t_meth <- file[[ii]]$methpos %>% intersect(colnames(data_filter_t[["meth"]]))
  truth <- list(posMeth= t_meth, posRNA= file[[ii]][[2]] %>% data.frame %>% dplyr::select(contains("DE")) %>% colnames, posProt=NA)

  ##### MOCLUSTER
  fit <- mo_results$fit
  roc_mo  <- doROC_Moa(truth, fit)
  ROCValues_mo[[ii]] <- ComputepAUC(roc_mo)
  ##### SGCCA

  fit <- sgcca_results$fit
  #roc_sgcca <- doROC_SGCCA(truth, fit)
  #ROCValues_sgcca[[ii]] <- ComputepAUC(roc_sgcca)

  fit <- intNMF_results$fit
  #roc_nmf <- doROC_intNMF(truth, fit)
  #ROCValues_intNMF[[ii]] <- ComputepAUC(roc_nmf)

  fit <- icluster_results$fit
  #roc_iclust <- doROC_iCluster(truth, fit, data_filter_t)
  #ROCValues_icluster[[ii]] <- ComputepAUC(roc_iclust)
}

saveRDS(ROCValues_sgcca, file.path(pathMeth_sub, "2011sGCCA_ROC.rds"))
saveRDS(ROCValues_mo, file.path(pathMeth_sub, "2011mocluster_ROC.rds"))
saveRDS(ROCValues_intNMF, file.path(pathMeth_sub, "2011intNMF_ROC.rds"))
saveRDS(ROCValues_icluster, file.path(pathMeth_sub, "icluster_ROC.rds"))

pathMeth_sub <- "../data/extdata_PintMF/res_v2"

ROCValues<- lapply(list.files(pathMeth_sub, pattern ="my_meth_ROC", full.names=TRUE), readRDS)
ROCValues_sgcca <- readRDS(file.path(pathMeth_sub, "2011sGCCA_ROC.rds"))
ROCValues_mo <- readRDS(file.path(pathMeth_sub, "2011mocluster_ROC.rds"))
ROCValues_intNMF <- readRDS(file.path(pathMeth_sub, "2011intNMF_ROC.rds"))
ROCValues_icluster <- readRDS(file.path(pathMeth_sub, "icluster_ROC.rds"))

roc_Pint <- do.call(rbind, ROCValues) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="PIntMF")
roc_Moclust <- do.call(rbind, ROCValues_mo) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="MoCluster")
label_names <- c(
  `dataset 1` = "Methylation",
  `dataset 2` = "Gene expression"
)

roc_sgcca <- do.call(rbind, ROCValues_sgcca) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="SGCCA")
roc_nmf<- do.call(rbind, ROCValues_intNMF) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="intNMF")
roc_iclust <- do.call(rbind, ROCValues_icluster) %>% as.data.frame %>% tidyr::gather(key="dataset",value="AUC") %>% mutate(method="icluster")

ROC_df <- rbind(roc_Pint, roc_Moclust, roc_sgcca, roc_nmf, roc_iclust)

library(ggplot2)
my.cols <- c(RColorBrewer::brewer.pal(7, name = "Set2")[-1], "red")

ROC_df <- ROC_df %>% mutate(method=factor(method, levels=c("MoCluster", "SGCCA", "intNMF", "icluster", "PIntMF")))
g <- ggplot(na.omit(ROC_df), aes(x=method, y=AUC))+geom_boxplot()+facet_wrap(.~dataset,labeller = as_labeller(label_names),  scales = "free_y")+theme_bw()+xlab("")+theme( legend.position = "none",axis.text.x = element_text(angle=90, color=c("#000000" , "#000000" , "#000000","#000000", "#000000", size=15)),axis.text.y = element_text(size=10))+scale_x_discrete(labels=c("PIntMF"=expression(bold(PIntMF))))
g+ylim(c(0.50,1))
g %>% ggsave(filename="../Figs/OMICSsimla_new_sim_roc.eps", height=3, width=6)

library(tsutils)
ROC_df %>% filter(dataset=='dataset 1') %>% mutate(id= rep(1:50, 5)) %>%
  dplyr::select(method, AUC, id) %>% tidyr::spread(key=method, value=AUC) %>%
  select(-id) %>% as.matrix %>%
  nemenyi(plottype = "mcb")

library(rstatix)
library(ggpubr)
my_comparisons <- list( c("PIntMF", "MoCluster"), c("PIntMF", "SGCCA"), c("PIntMF", "intNMF"),c("PIntMF", "icluster") )

ROC_df_1 <- ROC_df %>% filter(dataset=='dataset 1')
g1_pwc <- ggboxplot(ROC_df_1, x = "method", y = "AUC", x.text.angle = 45,
                    main = "", xlab = FALSE)+
  stat_compare_means(comparisons = my_comparisons) # Add pairwise comparisons p-value

ROC_df_2 <- ROC_df %>% filter(dataset=='dataset 2')
g2_pwc <- ggboxplot(ROC_df_2, x = "method", y = "AUC", x.text.angle = 45,
                    main = "", xlab = FALSE, ylab=FALSE) +
  stat_compare_means(comparisons = my_comparisons)

g1_pwc%>% ggsave(filename="../Figs/OMICSsimla_new_sim_boxplot_pval_meth.pdf", height=5, width=4)
g2_pwc%>% ggsave(filename="../Figs/OMICSsimla_new_sim_boxplot_pval_expr.pdf", height=5, width=4)

### Best model 20/11
ii=1
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
pint_p <- lapply(2:7, function(p) {
  SolveInt( data_filter_t , p=p, max.it=5, flavor_mod = "glmnet", init_flavor = "snf")
}
)
pves <- sapply(pint_p, function (rr){
  rr$pve[length(rr$pve)]
})

bics <- mapply( function(res, Y) (computeLL(Y=Y, res=res) %>% compute_BIC(pen=0)),
                pint_p, list(data_filter_t))

sub_bic <-(t(bics) %>% as.data.frame)
g_bic <- sub_bic %>% ggplot(aes(x=p.p, y=BIC.n)) +geom_point()+ theme(legend.position="bottom")+scale_x_continuous()+theme_bw()+theme_bw()+ theme(legend.position="bottom")+scale_x_continuous(
  breaks = 1:10)+xlab("P")+ylab("MSE")+geom_line()
g_bic
coph <- sapply(pint_p, function (rr) compute_coph(rr))


df_pve <- data.frame(pve=pves, iteration=2:(length(pves)+1))
g_pve <- ggplot(df_pve, aes(x=iteration, y=pve))+geom_point()+theme_bw()+
  geom_line()+scale_x_continuous(breaks=1:10)+xlab("P")

df_coph <- data.frame(cophenetic=coph, iteration=2:(length(coph)+1))
g_c <- ggplot(df_coph, aes(x=iteration, y=cophenetic))+geom_point()+
  theme_bw()+geom_line()+scale_x_continuous(breaks=1:10)+xlab("P")

g <- gridExtra::grid.arrange(g_bic, g_pve,g_c, ncol=3)
ggsave("../Figs/OMICSSMILA_bics_pve.pdf", g, width =7, height=3.5)

