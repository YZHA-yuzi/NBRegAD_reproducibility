########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# generating data used in simulation studies
# please make sure the working directory is the folder where 
# the R script "FUNs.R" and 
# adjacency matrix for 48 states "AdjMatrix_states.RData" are saved.
########################################################################

# Load packages
library(MASS)
library(igraph)
library(Matrix)
library(dplyr)
library(splines)

# source self-defined R functions
source("FUNs.R")
# Read in basis coefs for natural cubic splines to generate data assuming
# time-varying associations
parm.ns <- readRDS("Coefs_ns.RData")
# Read in the adjacency for 48 states 
W.s.true <- readRDS(file = "AdjMatrix_states.RData")
Dw.s.true <- Diagonal(x = apply(W.s.true, 1, sum))
states.sel <- rownames(W.s.true)

# Create a directory to save simulated datasets under different scenarios
dir.create("./data")

##########################################################################
# I. GENERATE DATA WHILE ASSUMING THE PERCENT POSITIVE AND OBSERVED COUNTS 
# ARE INDEPENDENT
##########################################################################
# Set random seed
set.seed(12345)

# The number of replicates under each scenario
nsims = 100
# The number of locations = 48
nlocs = nrow(W.s.true)

# Initialize a list to store simulated datasets
dat.sim.list <- list()
for(isim in 1:nsims){
  
  dat.sim.list.i <- list()
  counti = 1
  ### -------- Step 1: generate percent positive ------- ###
  ### additive residuals ###
  ## time-specific random effects 
  ntimes.true = 100
  tau2.p.true <- 0.15
  rho.p.true <- 0.999
  W.p.true <- compute.adjmat(ntimes = ntimes.true + 2, str = "AR1")
  Dw.p.true <- Matrix::Diagonal(x = apply(W.p.true, 1, sum))
  e.p.t.true <- as.numeric( MASS::mvrnorm(1, mu = rep(0, ntimes.true + 2),
                                    Sigma = tau2.p.true*
                            solve(Dw.p.true - rho.p.true*W.p.true)) )
  
  ## state-specific random effects 
  tau2.p.s.true <- 0.1
  rho.p.s.true <- 0.9
  e.p.s.true <- as.numeric( MASS::mvrnorm(1, mu = rep(0, nlocs),
                                          Sigma = tau2.p.s.true*
                           solve(Dw.s.true - rho.p.s.true*W.s.true)) )
  
  ## combine together to get spatial-temporal data ##
  e.p.t.mat <- data.frame(week.index = 9:110,
                          e.t = e.p.t.true)
  e.p.s.mat <- data.frame(states = states.sel,
                          e.s = e.p.s.true)
  parm.prob.sel0 <- data.frame(states = rep(states.sel, length(9:110)),
                    week.index = as.numeric(rep(9:110, each = nlocs)))
  parm.prob.sel0 <- dplyr::left_join(parm.prob.sel0, e.p.t.mat, 
                              by = c("week.index" = "week.index")) %>% 
    dplyr::left_join(e.p.s.mat, by = c("states" = "states"))
  
  
  expit <- function(x){exp(x)/(1 + exp(x))}
  parm.prob.sel <- data.frame(alpha = -1, parm.prob.sel0)
  parm.prob.sel$perc <- expit(parm.prob.sel$alpha + 
                                parm.prob.sel0$e.t + 
                                parm.prob.sel0$e.s)
  
  ### get lag1 and lag2 ###
  parm.prob.sel <- parm.prob.sel %>% 
    group_by(states) %>% 
    mutate(lag1 = dplyr::lag(perc, n=1),
           lag2 = dplyr::lag(perc, n=2))
  parm.prob.sel <- parm.prob.sel[!is.na(parm.prob.sel$lag2), ]

  ### -------- Step 2: generate baseline death counts ------- ##
  ntimes.true = 100
  W.t.true <- compute.adjmat(ntimes = ntimes.true, str = "AR1")
  Dw.t.true <- Diagonal(x = apply(W.t.true, 1, sum))
  beta0.true <- -0.7 
  tau2.beta0.t.true <- 0.05^2
  rho.beta0.t.true <- 0.999
  e.t.mat <- data.frame(week.index = 11:110, 
                        e.count.t = 
                          as.numeric( MASS::mvrnorm(1, mu = rep(0, ntimes.true),
                                    Sigma = tau2.beta0.t.true*
                          solve(Dw.t.true - rho.beta0.t.true*W.t.true)) ))
  tau2.beta0.s.true = 0.5
  rho.beta0.s.true = 0.9
  e.s.mat <- data.frame(states = states.sel, 
                        e.count.s = 
                          as.numeric( MASS::mvrnorm(1, mu = rep(0, nlocs),
                                     Sigma = tau2.beta0.s.true*
                          solve(Dw.s.true - rho.beta0.s.true*W.s.true)) ))
  
  dat.sim0 <- parm.prob.sel %>% 
    dplyr::left_join(e.t.mat, by = c("week.index" = "week.index")) %>%
    dplyr::left_join(e.s.mat, by = c("states" = "states"))
  
  eta0.true <- as.numeric(beta0.true + dat.sim0$e.count.t + dat.sim0$e.count.s)
  q0.true = 1/(1+exp(eta0.true)) # 1 - Pr(success)
  xi.true = 200
  ## simulate y0 from Poisson ##
  y0.true <- rpois(n = length(eta0.true),
                   lambda = xi.true*exp(eta0.true))

  for(j in c("cons", "tvary")){
    
    ### Step 3: generate baseline deaths + under-reported COVID-19 deaths ###
    ### beta1 and beta2 are parametric functions ###
    if(j == "cons"){
      beta1.true <- 1
      beta2.true <- 0.8
    }else if(j == "tvary"){
      X.design = cbind(1, ns(11:110, df = 5))
      beta1.true <- as.numeric(X.design %*% 
                                 matrix(parm.ns$coef.beta1, ncol = 1))
      beta2.true <- as.numeric(X.design %*% 
                                 matrix(parm.ns$coef.beta2, ncol = 1))
    }
    
    dat.sim1 <- dat.sim0 %>% left_join(data.frame(week.index = 11:110,
                                                  beta1 = beta1.true,
                                                  beta2 = beta2.true),
                                       by = c("week.index" = "week.index"))
    dat.sim1 <- dat.sim1 %>% mutate(pred.vals = beta1*lag1 + beta2*lag2)
    
    y.obs <- rnbinom(n = nrow(dat.sim1), 
                     size = xi.true, 
                     prob = 1/(1+exp(eta0.true + dat.sim1$pred.vals)))
    
    ## Step 4: get unobserved death counts ##
    y.CV19.true = y.obs - y0.true
    dat.sim <- data.frame(dat.sim1, 
                          y.obs = round(y.obs),
                          y0.true = y0.true,
                          y.CV19.true = y.CV19.true)
    
    dat.sim.list.i[[counti]] <- dat.sim
    counti  = counti + 1
  }
  names(dat.sim.list.i) <- c("cons", "tvary")
  dat.sim.list[[isim]] <- dat.sim.list.i
} # end loop over simulations
save(dat.sim.list, file = paste0("./data/dat_inde_ns_st.rda"))



##########################################################################
# II. GENERATE DATA WHILE ASSUMING THE PERCENT POSITIVE AND OBSERVED COUNTS 
# ARE CORRELATED
# Levels of the strength of temporal correlation between the percent positive
# and observed counts are:
# (1). strong positive ("pos_s")
# (2). moderate positive ("pos_m")
# (3). strong negative ("nega_s")
# (4). moderate negative ("nega_m")
##########################################################################

# Set random seed
set.seed(12345)
nsims = 100
nlocs = nrow(W.s.true)

for(k in c("pos_s", "pos_m", "nega_s", "nega_m")){
  
  dat.sim.list <- list()
  for(isim in 1:nsims){
    
    ntimes.true = 100
    ### Step 0: generate two correlated temporal process ###
    ### Kim, Sun, and Tsutakawa (2001) 
    ### multivariate CAR model in the bivariate case 
    rho1 = 0.999; rho2 = 0.999
    if(k == "pos_s"){
      rho0 = 0.9; rho3 = 0.9 # strong positively correlated
    }else if(k == "pos_m"){
      rho0 = 0.5; rho3 = 0.5 # moderate positively correlated 
    }else if(k == "nega_s"){
      rho0 = -0.9; rho3 = -0.9 # strong negatively correlated
    }else if(k == "nega_m"){
      rho0 = -0.5; rho3 = -0.5 # moderate negatively correlated
    } 
    
    W.p.true <- compute.adjmat(ntimes = ntimes.true + 2, str = "AR1")
    Dw.p.true <- Matrix::Diagonal(x = apply(W.p.true, 1, sum))
    tau1 = 1; tau2 = 25
    I <- Matrix::Diagonal(x = rep(1, ntimes.true + 2))
    comp11 <- (2*Dw.p.true + I - rho1*W.p.true)*tau1
    comp12 <- -(rho0*I + rho3*W.p.true)*sqrt(tau1*tau2)
    comp22 <- (2*Dw.p.true + I - rho2*W.p.true)*tau2
    preMCAR <- solve(rbind(cbind(comp11, comp12),
                           cbind(comp12, comp22)))
    
    # residuals for percent cent and residuals for deaths counts 
    e.MCAR <- as.numeric( MASS::mvrnorm(1, mu = rep(0, 2*(ntimes.true+2)),
                                  Sigma = preMCAR) )
    
    
    ### -------- Step 1: generate percent positive ------- ###
    ### additive residuals ###
    ## time-specific random effects 
    e.p.t.true <- as.numeric( e.MCAR[1:(ntimes.true + 2)] )
    
    ## state-specific random effects 
    tau2.p.s.true <- 0.1
    rho.p.s.true <- 0.9
    e.p.s.true <- as.numeric( MASS::mvrnorm(1, mu = rep(0, nlocs),
                              Sigma = tau2.p.s.true*
                              solve(Dw.s.true - rho.p.s.true*W.s.true)) )
    
    ## combine together to get spatial-temporal data ##
    e.p.t.mat <- data.frame(week.index = 9:110,
                            e.t = e.p.t.true)
    e.p.s.mat <- data.frame(states = states.sel,
                            e.s = e.p.s.true)
    parm.prob.sel0 <- data.frame(states = rep(states.sel, length(9:110)),
                      week.index = as.numeric(rep(9:110, each = nlocs)))
    parm.prob.sel0 <- dplyr::left_join(parm.prob.sel0, e.p.t.mat, 
                                by = c("week.index" = "week.index")) %>% 
      dplyr::left_join(e.p.s.mat, by = c("states" = "states"))
    
    
    expit <- function(x){exp(x)/(1 + exp(x))}
    parm.prob.sel <- data.frame(alpha = -1, parm.prob.sel0)
    parm.prob.sel$perc <- expit(parm.prob.sel$alpha + 
                                  parm.prob.sel0$e.t + 
                                  parm.prob.sel0$e.s)
    
    ### get lag1 and lag2 ###
    parm.prob.sel <- parm.prob.sel %>% 
      group_by(states) %>% 
      mutate(lag1 = dplyr::lag(perc, n=1),
             lag2 = dplyr::lag(perc, n=2))
    parm.prob.sel <- parm.prob.sel[!is.na(parm.prob.sel$lag2), ]
    
    ### -------- Step 2: generate baseline death counts ------- ##
    ntimes.true = 100
    W.t.true <- compute.adjmat(ntimes = ntimes.true, str = "AR1")
    Dw.t.true <- Diagonal(x = apply(W.t.true, 1, sum))
    beta0.true <- -0.7 
    # log(1000) - log(xi), where 1000 is the desired mean of the counts 
    tau2.beta0.t.true <- 0.05^2
    rho.beta0.t.true <- 0.999
    e.t.mat <- data.frame(week.index = 11:110, 
                          e.count.t = as.numeric( 
                e.MCAR[(ntimes.true + 2 + 1 + 2):(2*(ntimes.true + 2))] ))
    tau2.beta0.s.true = 0.5
    rho.beta0.s.true = 0.9
    e.s.mat <- data.frame(states = states.sel, 
                          e.count.s = 
                            as.numeric( MASS::mvrnorm(1, mu = rep(0, nlocs),
                                      Sigma = tau2.beta0.s.true*
                          solve(Dw.s.true - rho.beta0.s.true*W.s.true)) ))
    
    dat.sim0 <- parm.prob.sel %>% 
      dplyr::left_join(e.t.mat, by = c("week.index" = "week.index")) %>%
      dplyr::left_join(e.s.mat, by = c("states" = "states"))
    
    eta0.true <- as.numeric(beta0.true + dat.sim0$e.count.t + dat.sim0$e.count.s)
    q0.true = 1/(1+exp(eta0.true)) # 1 - Pr(success)
    xi.true = 200
    ## simulate y0 from Poisson ##
    y0.true <- rpois(n = length(eta0.true),
                     lambda = xi.true*exp(eta0.true))
    
    dat.sim.list.i <- list()
    counti = 1
    for(j in c("cons", "tvary")){
      
      ### Step 3: generate baseline deaths + under-reported COVID-19 deaths ###
      ### beta1 and beta2 are parametric functions ###
      if(j == "cons"){
        beta1.true <- 1
        beta2.true <- 0.8
      }else if(j == "tvary"){
        X.design = cbind(1, ns(11:110, df = 5))
        beta1.true <- as.numeric(X.design %*% 
                                   matrix(parm.ns$coef.beta1, ncol = 1))
        beta2.true <- as.numeric(X.design %*% 
                                   matrix(parm.ns$coef.beta2, ncol = 1))
      }
      
      dat.sim1 <- dat.sim0 %>% dplyr::left_join(data.frame(week.index = 11:110,
                                                    beta1 = beta1.true,
                                                    beta2 = beta2.true),
                                         by = c("week.index" = "week.index"))
      dat.sim1 <- dat.sim1 %>% mutate(pred.vals = beta1*lag1 + beta2*lag2)
      y.obs <- rnbinom(n = nrow(dat.sim1), 
                       size = xi.true, 
                       prob = 1/(1+exp(eta0.true + dat.sim1$pred.vals)))
      
      ## Step 4: get unobserved death counts ##
      y.CV19.true = y.obs - y0.true
      dat.sim <- data.frame(dat.sim1, 
                            y.obs = round(y.obs),
                            y0.true = y0.true,
                            y.CV19.true = y.CV19.true)
      dat.sim.list.i[[counti]] <- dat.sim
      counti = counti + 1
    }
    
    names(dat.sim.list.i) <- c("cons", "tvary")
    dat.sim.list[[isim]] <- dat.sim.list.i
  } # end loop over simulations
  
  save(dat.sim.list, 
       file = paste0("./data/dat_", k, "_ns_st.rda"))
  cat(paste(k, j, sep = ","))
  
} # end loop over temporal dependency





