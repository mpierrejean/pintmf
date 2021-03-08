library(dplyr)
library(CrIMMix)
library(future)
library(purrr)
library(stringr)
library(PintMF)
library(ggplot2)
source("inst/Simulations/00.setup.R")

pathDat <- "../data/extdata_Crimmix/Data_sim_20181012/"
listBenchmark <- list.files(pathDat)

df_all_bench <- lapply(listBenchmark, function (bb){
  pp <- file.path(pathDat, bb)
  list_simu <- sprintf("simu%s.rds", 1:25)
  l_sim <- file.path(pp, list_simu)
  do.call(rbind, lapply(l_sim, function(ss){
    dat <- readRDS(ss)
    num_sim <- str_extract(str_remove(ss,".*simu"), "[1-2]")
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
    data_filter <- dat$data %>% remove_zero
    data_t_filter <- data_filter
    data_t_filter[[3]] <- log2(data_t_filter[[3]]/(1-data_t_filter[[3]]))
    print("my_meth")
    my_meth_results <- lapply(2:10, function(K) data_t_filter %>% SolveInt( p=K, max.it=5, flavor_mod = "glmnet", init_flavor = "snf"))
    df_crit_ms <- rbind(sapply(my_meth_results, function(res){
      pve = res$pve[length(res$pve)]
      RSS = computeLL(Y=data_t_filter, res=res)
      coph = cor(hclust(dist(res$W), method="ward.D2") %>% cophenetic, dist(res$W))
      return(list(pve=pve, rss=RSS[1], coph=coph))
    }), p=2:10) %>% t %>% as.data.frame() %>% tidyr::gather(key = crit, value = val, pve:coph) %>%
      mutate(p=as.numeric(p), crit=as.factor(crit), val=as.numeric(val), simu=num_sim)
  })) %>% mutate(Benchmark= bb)
})

rr <- do.call(rbind,df_all_bench)
rr$simu <- rep(1:25, each=9*3)
saveRDS(rr, "../data/extdata_PintMF/model_selection.rds")
rr <- readRDS("../data/extdata_PintMF/model_selection.rds")

crit_names <- c(
  `coph` = "Cophenetic coefficient",
  `pve` = "PVE",
  `rss` = "MSE"
)
library(dplyr)
sapply(listBenchmark, function (bb) {
  xinter=4
  if(str_extract(bb, "[0-9]")==6){
    xinter <- 2
  }
  if(str_extract(bb, "[0-9]")==7){
    xinter <- 3
  }
  gg <- rr %>% filter(Benchmark==bb)%>% ggplot(aes(x=p, y=val, group=simu))+
    geom_line(alpha=0.5)+
  facet_wrap(.~crit, scale="free", labeller = as_labeller(crit_names))+
    theme_bw()+geom_vline(xintercept=xinter, col = "red")+
  xlim(c(2,10))+xlab('P')+ylab('Value')+theme(
    strip.text.x =element_text(size=15),
    axis.title.x = element_text(size=15),
    axis.text.x = element_text(size=13),
  )
  ggsave(filename = sprintf("../Figs/%s_crit_test.pdf", bb), width=10, height=5)
})
