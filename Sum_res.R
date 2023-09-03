########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# Summarize simulation results across a total of 10 scenarios
# !!! please make sure the working directory is the folder where 
# R script named "FUNs.R" is saved. 
########################################################################

# Load in packages
library(writexl)
library(splines)
library(dplyr)

# Load in self-defined functions
source("FUNs.R")

##################################################
################ SUMMARIZE RESULTS ###############
##################################################
### 1. waic 
### 2. betas (constant and tvarying)
### 3. hyper parameters
### 4. counts (yobs, y0, AD)

beta0.true <- -0.7
xi.true <- 200
beta1.true.cons <- 1
beta2.true.cons <- 0.8
parm.true <- c(0.05^2, 0.999, 200, 0.5, 0.9)
parm.ns <- readRDS("Coefs_ns.RData")
X.design = cbind(1, ns(11:110, df = 5))
beta1.true.tvary <- as.numeric(X.design %*% 
                                 matrix(parm.ns$coef.beta1, ncol = 1))
beta2.true.tvary <- as.numeric(X.design %*% 
                                 matrix(parm.ns$coef.beta2, ncol = 1))

 
count = 1
waic.freq <- data.frame(dependency = rep(NA, 10),
                        beta = rep(NA, 10),
                        freq = rep(NA, 10))
count.sum.all <- list()
beta.bar.list <- list()
parm.list <- list()
beta.rand.pointwise.list <- list()
beta.rand.sum.list <- list()
names.vec <- c()
for(k in c("inde", "pos_s", "pos_m", "nega_s", "nega_m")){
  for(j in c("cons", "tvary")){
    res.all <- list()
    for(i in 1:100){
      if(k == "inde"){
        load(paste0("./inter_res/SIMs_res_", "inde", "_",
                    j, "_ns_st_", i, ".rda"))
      }else{
        load(paste0("./inter_res/SIMs_res_corre_", k, "_",
                    j, "_ns_st_", i, ".rda"))
      }
      res.all[[i]] <- re.curr
    }
    
    waic.mat <- do.call(rbind, lapply(res.all, 
                                      function(x) x[["waic"]][["waic"]]))
    pwaic.mat <- do.call(rbind, lapply(res.all, 
                                       function(x) x[["waic"]][["pwaic"]]))
    
    # get frequency of selecting models #
    waic.freq[count, 1:2] <- c(k, j)
    waic.freq[count, 3] <- mean(waic.mat[,1] > waic.mat[,2])
    
    beta.bar.mat <- get_sumparms(res.all = res.all, 
                                 true.vec = c(beta0.true, 
                                              beta1.true.cons, 
                                              beta2.true.cons, 
                                              beta1.true.cons + beta2.true.cons), 
                                 s1 = "beta")
    beta.bar.list[[count]] <- beta.bar.mat
    
    parm.mat <- get_sumparms(res.all = res.all, 
                             true.vec = parm.true, 
                             s1 = "parm")
    parm.list[[count]] <- parm.mat
    
    beta.rand.true.mat <- cbind(beta1.true.tvary, beta2.true.tvary,
                                beta1.true.tvary + beta2.true.tvary)
    beta.rand.mat0 <- get_sumbetarands(res.all = res.all,
                                       beta.rand.true.mat = beta.rand.true.mat)
    beta.rand.mat.pointwise <- beta.rand.mat0$pointwise
    beta.rand.mat.sum <- beta.rand.mat0$sum
    beta.rand.pointwise.list[[count]] <- beta.rand.mat.pointwise
    beta.rand.sum.list[[count]] <- beta.rand.mat.sum
    
    s1.1.vec <- c("total", "timeseries", "map", "st")
    s1.2.vec <- c("cons", "tvary")
    s1.3.vec <- c("pred.counts", "pred.null.counts", "AD",
                  paste0(c("pred.counts", "pred.null.counts", "AD"), ".randdraw"))
    counts.sum <- NULL
    for(s1.1 in s1.1.vec){
      for(s1.3 in s1.3.vec){
        for(s1.2 in s1.2.vec){
          counts.sum[[s1.1]][[s1.3]][[s1.2]] <- get_sumcounts(res.all = res.all,
                                                              s1.1 = s1.1, 
                                                              s1.2 = s1.2, 
                                                              s1.3 = s1.3)
        }
      }
    }
    count.sum.all[[count]] <- counts.sum
    names.vec <- c(names.vec, paste(k, j, sep = "_"))
    count = count + 1
    cat(count)
  } # end loop over beta coefs
} # end loop over different dependency 
names(count.sum.all) <- names(beta.bar.list) <- 
  names(parm.list) <- names(beta.rand.pointwise.list) <- 
  names(beta.rand.sum.list) <- names.vec


### Generate Table 1 
tab.print0 <- do.call(rbind, lapply(1:10, function(x) do.call(rbind, 
                      lapply(count.sum.all[[x]]$timeseries$AD.randdraw,
                             function(y) y$sum))))
tab.print0.pred <- do.call(rbind, lapply(1:10, function(x) 
  do.call(rbind, lapply(count.sum.all[[x]]$timeseries$pred.counts.randdraw,
                        function(y) y$sum))))

tab.print1 <- data.frame(coef.true = rep(c("cons", "cons", "tvary", "tvary"), 5),
                         dependency = rep(c("inde", "pos_s", "pos_m", 
                                            "nega_s", "nega_m"), each = 4),
                         coef.fit = rep(c("cons", "tvary"), 10),
                         rel.bias = formatC(tab.print0$avg.bias/tab.print0$avg.true,
                                            digits = 3, format = "f"),
                         cover = formatC(tab.print0$avg.cover*100,
                                         digits = 2, format = "f"),
                         rse = formatC(tab.print0$rse,
                                       digits = 3, format = "f"))
tab.print1$dependency <- factor(tab.print1$dependency, 
                                level = c("inde", 
                                          "pos_s", "pos_m",
                                          "nega_s", "nega_m"))
tab.print1 <- tab.print1[order(tab.print1$coef.true, 
                               tab.print1$dependency), ]


tab.beta.print0 <- do.call(rbind, lapply(beta.rand.sum.list, function(x) x[3, ]))
tab.beta.print1 <- data.frame(coef.true = rep(c("cons", "tvary"), 5),
                              dependency = rep(c("inde", "pos_s", "pos_m", 
                                                 "nega_s", "nega_m"), each = 2),
                              coef.fit = rep(c("tvary"), 10),
                              rel.bias = formatC(tab.beta.print0[,2]/tab.beta.print0[,1],
                                                 digits = 3, format = "f"),
                              cover = formatC(tab.beta.print0[,3]*100,
                                              digits = 2, format = "f"),
                              rse = formatC(tab.beta.print0[,5],
                                            digits = 3, format = "f"))
tab.beta.print1$dependency <- factor(tab.beta.print1$dependency, 
                                     level = c("inde", "pos_s", "pos_m",
                                               "nega_s", "nega_m"))
tab.beta.print1 <- tab.beta.print1[order(tab.beta.print1$coef.true, 
                                         tab.beta.print1$dependency), ]

tab.print.all <- left_join(tab.print1, tab.beta.print1, 
                           by = c("dependency" = "dependency",
                                  "coef.true" = "coef.true",
                                  "coef.fit" = "coef.fit")) 
colnames(tab.print.all)[-c(1:3)] <- 
  c(paste0("AD:", c("rel.bias", "cover", "RSE")),
    paste0("beta(t):", c("rel.bias", "cover", "RSE")))


# Create a directory to save simulated datasets under different scenarios
dir.create("./table")
writexl::write_xlsx(tab.print.all,
                    path = paste0("./table/Table1.xlsx"))

