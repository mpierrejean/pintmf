##########################################################################
## ROC evaluation
##########################################################################
library(dplyr)
library(CrIMMix)
library(future)
library(purrr)
library(stringr)
plan(multiprocess)
library(tis)
source("inst/Simulations/00.setup.R")

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
  names(aucs) <- c("dataset 1", "dataset 2", "dataset 3")
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
pathDat <- "../data/extdata/Data_sim_20181012/"
listBenchmark <- list.files(pathDat)
pathMeth <- "../data/extdata/Data_Results_20181012/"
library(future.apply)
forceEVAl= FALSE
if(forceEVAl){
  auc_eval_dat_my_meth <- do.call(rbind,lapply(1:8, function(ii){
    b <- listBenchmark[ii]
    pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))

    pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))

    print(pathMeth_sub)
    print(pathDat_sim)

    list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
    regexp <- "[[:digit:]]+"

    trueDat1 <- sapply(list.sim, function(ss) ss$biomark$dat1 %>% unlist %>% unique %>% str_extract(pattern=regexp))
    trueDat2 <-  sapply(list.sim, function(ss) ss$biomark$dat2 %>% unlist %>% unique %>% str_extract(pattern=regexp))
    trueDat3 <-  sapply(list.sim, function(ss)  ss$biomark$dat3 %>% unlist %>% unique %>% str_extract(pattern=regexp))

    truth <- lapply(1:S, function(ss) list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]]) )


    mm <- "my_meth"
    print(mm)
    pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
    ff <- readRDS(pp)
    H_list <- lapply(ff, function(f){
      H_coef <- lapply(f$H, function(hh) apply(abs(hh), 2, sum))
      H_sorted <- lapply(H_coef, function(h) sort(h,decreasing = TRUE) %>% names %>% str_extract(pattern=regexp))
    } )
    l <- list(t = truth, h = H_list)
    J <- sapply(list.sim[[1]]$data, ncol)
    roc_funct <- function(t, h, J){
      roc_eval <- list(FPR = FPR_compute(t, h, J), TPR = TPR_compute(t, h))
    }
    auc_eval_my_meth <- purrr::pmap(l, roc_funct, J = J) %>% sapply(ComputepAUC) %>% t
  }))

  auc_eval_dat_my_meth %>% head
  auc_eval_dat_my_meth  <- auc_eval_dat_my_meth %>% as.data.frame %>% mutate(method="PIntMF", noise=rep(sprintf("Benchmark%s", 1:8), each=50))
  saveRDS(auc_eval_dat_my_meth, "../data/extdata/Data_Results_20181012/AUC_eval_PintMF.rds")

}

forceEVAlsupervised <- FALSE
if(forceEVAlsupervised){
  auc_eval_dat_my_meth <- do.call(rbind,lapply(1:8, function(ii){
    b <- listBenchmark[ii]
    pathMeth_sub <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathMeth, b))

    pathDat_sim <- R.utils::Arguments$getWritablePath(sprintf("%s/%s", pathDat, b))

    print(pathMeth_sub)
    print(pathDat_sim)

    list.sim <- list.files(pathDat_sim, full.names = TRUE) %>% lapply(readRDS)
    regexp <- "[[:digit:]]+"

    trueDat1 <- sapply(list.sim, function(ss) ss$biomark$dat1 %>% unlist %>% unique %>% str_extract(pattern=regexp))
    trueDat2 <-  sapply(list.sim, function(ss) ss$biomark$dat2 %>% unlist %>% unique %>% str_extract(pattern=regexp))
    trueDat3 <-  sapply(list.sim, function(ss)  ss$biomark$dat3 %>% unlist %>% unique %>% str_extract(pattern=regexp))

    truth <- lapply(1:S, function(ss) list(trueDat1[[ss]], trueDat2[[ss]], trueDat3[[ss]]) )


    mm <- "my_meth_res_supervised"
    print(mm)
    pp <- list.files(pathMeth_sub, pattern = mm, full.names = TRUE)
    ff <- readRDS(pp)
    H_list <- lapply(ff, function(f){
      H_coef <- lapply(f$H, function(hh) apply(abs(hh), 2, sum))
      H_sorted <- lapply(H_coef, function(h) sort(h,decreasing = TRUE) %>% names %>% str_extract(pattern=regexp))
    } )
    l <- list(t = truth, h = H_list)
    J <- sapply(list.sim[[1]]$data, ncol)
    roc_funct <- function(t, h, J){
      roc_eval <- list(FPR = FPR_compute(t, h, J), TPR = TPR_compute(t, h))
    }
    auc_eval_my_meth <- purrr::pmap(l, roc_funct, J = J) %>% sapply(ComputepAUC) %>% t
  }))

  auc_eval_dat_my_meth %>% head
  auc_eval_dat_my_meth  <- auc_eval_dat_my_meth %>% as.data.frame %>% mutate(method="PIntMF supervised", noise=rep(sprintf("Benchmark%s", 1:8), each=50))
  saveRDS(auc_eval_dat_my_meth, "../data/extdata/Data_Results_20181012/AUC_eval_PintMF_supervised.rds")

}

auc_eval_dat_my_meth_supervised <- readRDS("../data/extdata_Crimmix/Data_Results_20181012/AUC_eval_PintMF_supervised.rds")
auc_eval_dat_my_meth <- readRDS("../data/extdata_Crimmix/Data_Results_20181012/AUC_eval_PintMF.rds")
auc_eval_dat <- readRDS("../data/extdata_Crimmix/Data_Results_20181012/AUC_eval.rds")
auc_eval_dat <- rbind(auc_eval_dat, auc_eval_dat_my_meth)
auc_eval_dat_2 <- auc_eval_dat %>% group_by(noise, method) %>% summarize(MeanD1 = mean(`dataset 1`), MeanD2 = mean(`dataset 2`), MeanD3= mean(`dataset 3`))

x1 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD1) %>% spread(method, MeanD1)
x2 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD2) %>% spread(method, MeanD2)
x3 <- auc_eval_dat_2 %>% dplyr::select(method, noise, MeanD3) %>% spread(method, MeanD3)
rbind(x1,x2,x3) %>% xtable::xtable()
auc_eval_dat %>% group_by(noise, method) %>% filter(noise%in% paste("Benchmark", 1:8, sep=""),method%in%c("CIMLR", "icluster", "iNMF", "MoCluster", "PIntMF", "SGCCA", "PIntMF supervised"))%>% summarize(Gaussian= mean(`dataset 1`) %>% round(2), Binary= mean(`dataset 2`)%>% round(2), "Beta-like"= mean(`dataset 3`)%>% round(2)) %>%
  pivot_longer(
    cols = Gaussian:`Beta-like`,
    names_to = "Dataset",
    values_to = "AUC",
    values_drop_na = TRUE
  )%>%spread(key = method, value=AUC) %>% xtable::xtable()
auc_eval_dat %>% group_by(noise, method) %>% filter(noise%in% paste("Benchmark", 1:8, sep=""),method%in%c("CIMLR", "icluster", "iNMF", "MoCluster", "PIntMF", "SGCCA", "PIntMF supervised"))%>% summarize(Gaussian= sd(`dataset 1`) %>% round(2), Binary= sd(`dataset 2`)%>% round(2), "Beta-like"= sd(`dataset 3`)%>% round(2)) %>%
  pivot_longer(
    cols = Gaussian:`Beta-like`,
    names_to = "Dataset",
    values_to = "AUC",
    values_drop_na = TRUE
  )%>%spread(key = method, value=AUC) %>% xtable::xtable()

auc_eval_dat_2$method <- gsub("iNMF", "intNMF", auc_eval_dat_2$method)
auc_eval_dat_2$method <- gsub("icluster", "iClusterPlus", auc_eval_dat_2$method)

library(ggplot2)
my.cols <- c(RColorBrewer::brewer.pal(8, name = "Set2")[-c(1, 7)])
g_ROC <- auc_eval_dat_2 %>% filter(noise%in% paste("Benchmark", 1:8, sep=""),method%in%c("CIMLR", "iClusterPlus", "intNMF", "MoCluster", "PIntMF", "SGCCA", "SNF"))%>%
  mutate(method=factor(method, levels =
                         c("intNMF",
                           "SGCCA",
                           "MoCluster",
                           "iClusterPlus",
                           "CIMLR","PIntMF"
                         ))) %>% rename(Gaussian = MeanD1, Binary=MeanD2, "Beta-like"=MeanD3)   %>%
  tidyr::gather(key="data", value="AUC", "Gaussian":"Beta-like") %>%mutate(data=factor(data, levels=c("Gaussian", "Binary", "Beta-like"))) %>%
  ggplot(aes(x=method, y=AUC, fill=method))+geom_bar(stat='identity', show.legend = FALSE)+facet_grid(data~noise)+theme_bw()+scale_y_continuous(breaks=c(0, 0.5, 1))+
  scale_fill_brewer(palette = "Set2")+
  theme(axis.text.x = element_text(size = 12, angle = 90, color=c("#000000" , "#000000" , "#000000" , "#000000" , "#000000" , "#993333" )), strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 18), legend.position = "none", axis.title= element_blank(),  axis.text.y=element_text(size=13))+scale_x_discrete(labels=c("PIntMF"=expression(bold(PIntMF)), "PIntMF supervised"=expression(bold(PIntMFsup)), parse=TRUE))
g_ROC


ggsave(filename=file.path("../../papers/NAR-LaTeX-20170208/", "auc_b1to4.eps"),g_ROC, width=12, height=7)



auc_eval_dat %>% group_by(method) %>% filter(noise%in% paste("Benchmark", 1:4, sep=""),method%in%c("CIMLR", "icluster", "iNMF", "MoCluster", "PIntMF", "SGCCA"))%>% summarize(MeanD1= mean(`dataset 1`) %>% round(2), MeanD2= mean(`dataset 2`)%>% round(2), MeanD3= mean(`dataset 3`)%>% round(2))



g_ROC <- auc_eval_dat_2 %>% filter(noise%in% paste("Benchmark", 1:8, sep=""),method%in%c("PIntMF", "PIntMF supervised"))%>%
  mutate(method=factor(method, levels =
                         c("PIntMF", "PIntMF supervised"
                         ))) %>% rename(Gaussian = MeanD1, Binary=MeanD2, "Beta-like"=MeanD3)   %>%
  tidyr::gather(key="data", value="AUC", "Gaussian":"Beta-like") %>%mutate(data=factor(data, levels=c("Gaussian", "Binary", "Beta-like"))) %>%
  ggplot(aes(x=method, y=AUC, fill=method))+geom_bar(stat='identity', show.legend = FALSE)+facet_grid(data~noise)+theme_bw()+scale_fill_manual(values=my.cols)+scale_y_continuous(breaks=c(0.8,0.85,0.90,0.95, 1))+
  theme(axis.text.x = element_text(size = 12, angle = 90, color=c("#993333", "#993333" )), strip.text.x = element_text(size = 15), strip.text.y = element_text(size = 15), legend.position = "none", axis.title= element_blank())+scale_x_discrete(labels=c("PIntMF"=expression(bold(PIntMF)), "PIntMF supervised"=expression(bold(PIntMFsup)), parse=TRUE))+ coord_cartesian(ylim=c(0.85,1.0))
g_ROC

ggsave(filename=file.path("../../papers/NAR-LaTeX-20170208/", "auc_supervised.eps"),g_ROC, width=15, height=10)

