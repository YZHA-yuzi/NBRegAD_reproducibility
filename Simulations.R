########################################################################
# --------------- R codes used in the simulation study --------------- #
# Fit the proposed model to the simulated data under different scenarios
# NOTE:
# !!! please make sure the working directory is the folder where 
# R script named "FUNs.R" is located at. 
# !!! please make sure the R script "Simulation_data.R" was ran first.
# !!! please make sure there is a sub-folder named "inter_res" in the
# working directory
########################################################################

args <-  commandArgs(trailingOnly = TRUE)
# the index for simulation 1 - 100
isim = eval( parse(text=args[1]) ) 
# the index for the strength of temporal correlations between the percent
# positive and observed counts 
# (1 = independent, 
# 2 = strong positive, 3 = moderate positive,
# 4 = strong negative, 5 = moderate negative)
sceindex = eval( parse(text=args[2]) ) 
# the index for the type of coefficients 
# (1 = constant, 2 = time-varying)
j.index = eval( parse(text=args[3]) ) 

# Load in packages
library(MASS) # mvrnorm 
library(igraph); # get adjacency matrix for temporal trend
library(Matrix) # Matrix, convert to sparse matrix
library(sparseMVN) # compute MVN density with sparse cov-var matrix
library(BayesLogit) # rpg, draw samples from PG 
library(truncnorm) # for truncated normal dist
library(dplyr)
library(mvtnorm)
library(Rcpp)
library(splines)

# Read in self-defined R functions
source("FUNs.R")
sourceCpp("FUNs_cpp.cpp")

# Read in the adjacency for 48 states 
W.s <- readRDS(file = "AdjMatrix_states.RData")

# Load in simulated data
sce.vec <- c("inde", "pos_s", "pos_m", "nega_s", "nega_m")
sce.curr <- sce.vec[sceindex]
load(paste0("./data/dat_", sce.curr, "_ns_st.rda"))

j = c("cons", "tvary")[j.index]
dat.sim.fit <- dat.sim.list[[isim]][[j]]
ntotal = nrow(dat.sim.fit)

##################################################
##### 1. fit model assuming constant effects #####
##################################################
## !!! ONE simulation takes about 100s !!!
## on a computer with Apple M2 Pro chip and 32 Gb memory 
set.seed(12345)
# start = proc.time()[3]
re.cons.1 <- fit.NB.st(formula = y.obs ~ lag1 + lag2,
                       data = dat.sim.fit,
                       offset = NULL,
                       niter = 10000, burn_in = 5000,
                       rand.int = TRUE,
                       str.int = list("week.index", "AR1"),
                       slope.tvarying = FALSE,
                       covariate.AD = c("lag1", "lag2"),
                       spatial.str = list("states", "CAR"),
                       adj.mat.s = W.s)
# proc.time()[3] - start
## complete the model fitting 
burn_in = 5000
parm.mat <- re.cons.1$parm
pred.counts.mat <- re.cons.1$pred.counts
pred.counts.null.mat <- re.cons.1$pred.null.counts
xi.post = as.numeric(parm.mat[-c(1:burn_in), "xi"])
pred.counts.nb.mat <- do.call(cbind, lapply(1:length(xi.post), 
                                            function(x) 
                       rnbinom(n = ntotal, size = xi.post[x], 
                               mu = pred.counts.mat[,x])))
pred.counts.null.pois.mat <- do.call(cbind, 
                                     lapply(1:length(xi.post), 
                                     function(x) rpois(n = ntotal, 
                                     lambda = pred.counts.null.mat[,x])))
re.cons.1$pred.counts.randdraw <- pred.counts.nb.mat
re.cons.1$pred.null.counts.randdraw <- pred.counts.null.pois.mat
re.cons.1$AD.randdraw <- pred.counts.nb.mat - pred.counts.null.pois.mat

###########################################################
##### 2. fit model assuming time-varying coefficients #####
###########################################################
## !!! ONE simulation takes about 220s !!!
## on a computer with Apple M2 Pro chip and 32 Gb memory
# start = proc.time()[3]
re.tvary.1 <- fit.NB.st(formula = y.obs ~ lag1 + lag2,
                        data = dat.sim.fit,
                        offset = NULL,
                        niter = 10000, burn_in = 5000,
                        rand.int = TRUE,
                        str.int = list("week.index", "AR1"),
                        slope.tvarying = TRUE, 
                        str.tvarying = list(c("lag1", "week.index", "AR1"),
                                            c("lag2", "week.index", "AR1")),
                        covariate.AD = c("lag1", "lag2"),
                        spatial.str = list("states", "CAR"), 
                        adj.mat.s = W.s)
# proc.time()[3] - start

### complete the model fitting ###
burn_in = 5000
parm.mat <- re.tvary.1$parm
pred.counts.mat <- re.tvary.1$pred.counts
pred.counts.null.mat <- re.tvary.1$pred.null.counts
xi.post = as.numeric(parm.mat[-c(1:burn_in), "xi"])
pred.counts.nb.mat <- do.call(cbind, lapply(1:length(xi.post), 
                                            function(x) 
                              rnbinom(n = ntotal, size = xi.post[x], 
                                      mu = pred.counts.mat[,x])))
pred.counts.null.pois.mat <- do.call(cbind, lapply(1:length(xi.post), 
                                                   function(x) 
                                                     rpois(n = ntotal, 
                             lambda = pred.counts.null.mat[,x])))
re.tvary.1$pred.counts.randdraw <- pred.counts.nb.mat
re.tvary.1$pred.null.counts.randdraw <- pred.counts.null.pois.mat
re.tvary.1$AD.randdraw <- pred.counts.nb.mat - pred.counts.null.pois.mat

## combine results from assuming both constant effects
## and time-varying effects
re.curr <- get_sumall(re1 = re.cons.1, re2 = re.tvary.1, 
                      dat.curr = dat.sim.fit,
                      burn_in = 5000)
if(sce.curr == "inde"){
  save(re.curr, file = paste0("./inter_res/SIMs_res_", sce.curr, "_",
                               j, "_ns_st_", isim, ".rda"))
}else{
  save(re.curr,
       file = paste0("./inter_res/SIMs_res_corre_", sce.curr, "_",
                     j, "_ns_st_", isim, ".rda"))
}

# if(sce.curr == "inde"){
#   ## combine results from assuming both constant effects 
#   ## and time-varying effects
#   res.inde <- get_sumall(re1 = re.cons.1, re2 = re.tvary.1, 
#                          dat.curr = dat.sim.fit,
#                          burn_in = 5000)
#   save(res.inde, file = paste0("./inter_res/SIMs_res_", sce.curr, "_", 
#                                j, "_ns_st_", isim, ".rda")) 
# }else{
#   ## combine results from assuming both constant effects 
#   ## and time-varying effects
#   res.corre <- get_sumall(re1 = re.cons.1, re2 = re.tvary.1, 
#                           dat.curr = dat.sim.fit,
#                           burn_in = 5000)
#   save(res.corre, 
#        file = paste0("./inter_res/SIMs_res_corre_", sce.curr, "_", 
#                      j, "_ns_st_", isim, ".rda"))
# }






