########################################################################
# ------ Self-defined R functions used in the simulation study --------#
########################################################################

# A function to get adjacency matrix 
compute.adjmat <- function(ntimes, str){
  if(str == "AR1"){
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = igraph::get.adjacency(graph.edgelist(adj.list, directed=FALSE))
  }else if(str == "AR2"){
    adj.list <- as.matrix(rbind(cbind(c(1:(ntimes-1)), c(2:ntimes)),
                                cbind(c(1:(ntimes-2)), c(3:ntimes))))
    adj.list <- adj.list[order(adj.list[,1]), ]
    W.t = igraph::get.adjacency(graph.edgelist(adj.list, directed=FALSE))
  }
  return(W.t)
}

get.randcovariate.list <- function(X.design.rand.slop,
                                   beta.rand.list,
                                   curr.index){
  num.ele = length(beta.rand.list)
  comp.rand.slop.list = lapply(1:num.ele,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][curr.index, ], ncol = 1)))
  # re <- Reduce("+", comp.rand.slop.list)
  return(comp.rand.slop.list)
}

### A function to fit space-time additive negative-binomial regression model
### to estimate attributable deaths (AD) based on spatial-temporal data
fit.NB.st <- function(formula, data, offset = NULL, 
                      niter = 5000, 
                      burn_in = 2500,
                      rand.int = TRUE, str.int = NULL, 
                      slope.tvarying = FALSE, str.tvarying = NULL, 
                      covariate.AD,
                      countsonly = FALSE,
                      spatial.str, adj.mat.s = NULL, 
                      ID.spacetime = NULL, ...){
  
  ### INPUTS:
  ## formula: y ~ x1 + x2
  ## data: a dataframe containing data
  ## offset: specify a variable that will be used as an offset term in the model
  ## niter: the number of MCMC iterations
  ## burn_in: the number of burn-in samples
  ## str.int: a list containing: (1) a variable name indicating 
  ###         the clusters of time-specific intercepts; (2) structure: exchangeable/AR-1/AR-2 
  ## slope.tvarying: a logic value indicating whether time-varying coefs are specifieds
  ## str.tvarying: a list list containing: (1) a variable name indicating 
  ###         the clusters of time-varying coefs; 
  ### (2) structure: exchangeable/AR-1/AR-2 can only be "exch", "AR1", "AR2" 
  ## covariate.AD: a vector of strings which specifies variable names corresponding 
  ### covariates that need to set to zero when computing attributable deaths  
  ## countsonly: a logic value indicating whether posterior samples of 
  ###   random effects (i.e., time-specific intercepts or time-varying coefs) are returned
  ## spatial.str: a list containing: (1) a variable name indicating 
  ### the clusters of location-specific intercepts; 
  ### (2) structure: exchangeable/CAR can only be "exch" or "CAR"
  ## adj.mat.s: an adjacency matrix for location-specific intercepts, must be provided when 
  ## spatial.str = c("", "CAR")
  
  
  
  pars.all <- as.list(match.call()[-1]) # remove the function name
  
  outcome.name <- all.vars(formula)[1]
  pred.names <- all.vars(formula)[-1]
  
  offset.name <- as.character(pars.all$offset)
  if(length(offset.name) != 0){
    offset.vec <- as.matrix(data[,offset.name])
    offset.vec <- as.numeric(offset.vec)
  }else{
    offset.vec <- NULL
  }
  outcome.vec <- as.matrix(data[,outcome.name])
  outcome.vec <- as.numeric(outcome.vec)
  num.obs <- length(outcome.vec)
  
  if(!all(covariate.AD %in% pred.names)){
    stop("variables must be included in the formula")
  }
  
  if(spatial.str[[2]] == "CAR" & is.null(adj.mat.s)){
    stop("an adjacency matrix for location-specific random effects must be specified")
  }
  
  if(!is.null(ID.spacetime)){
    if(!all(ID.spacetime %in% colnames(data))){
      stop("space and time indicators of observations must be inlucded in the dataset")
    }
  }
  
  if(is.null(ID.spacetime)){
    ID.spacetime.df <- NULL
  }else{
    ID.spacetime.df <- data[,ID.spacetime]
  }
  
  warning.text.a = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of location-specific random effects",
                          " and the structure of location-specific random effects")
  warning.text1.a = paste0("the structure can only be specified as ",
                           "exch, or CAR")
  warning.text2.a = paste0("the variable name indicating",
                           " the clusters of location-specific random effects",
                           " must be inlcuded in the dataset")
  
  if(is.null(spatial.str)){
    stop(warning.text.a)
  }
  if(!is.null(spatial.str)){
    if(!is.list(spatial.str)){
      stop(warning.text.a)
    }
    ll = length(spatial.str)
    if(ll != 2){
      stop(warning.text.a)
    }
    if(ll == 2){
      if(!(spatial.str[[2]] %in% c("exch", "CAR"))){
        stop(warning.text1.a)
      }
      if(!(spatial.str[[1]] %in% colnames(data))){
        stop(warning.text2.a)
      }
    }
  }
  id.loc.name = spatial.str[[1]]
  nlocs = length(unique(as.matrix(data[,id.loc.name])))
  if(nlocs == num.obs){
    X.design.spatial <- Diagonal( x = rep(1, num.obs))
  }else if(nlocs < num.obs){
    comp.vec <- paste0("(1|", id.loc.name, ")")
    formula.final <- as.formula(paste(c(deparse(formula), comp.vec), 
                                      collapse = "+"))
    X.all <- lme4::lFormula(formula = formula.final, data = data)
    X.design.spatial <- t(X.all$reTrms$Zt)
  }
  
  ###### SCENARIO 1: no temporal random effects and no time-varying coefs ######
  if(rand.int == FALSE & slope.tvarying == FALSE){
    X.design.fix <- model.matrix(formula, data = data)
    re0 <- fit.NB.st.add.s4(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix, 
                            offset.vec = offset.vec,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]], 
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  } # END SCENARIO 1
  
  
  ###### SCENARIO 2: no temporal random effects and time-varying coefs ######
  if(rand.int == FALSE & slope.tvarying == TRUE){
    warning.text.a = paste0("please provide a list containing",
                            " vectors that specify: ",
                            " (1) variable names of time-varying coefs,",
                            " (2) variable names indicating",
                            " the clusters of time-varying coefs",
                            " , and (3) the structure of time-varying coefs")
    warning.text.b = paste0("the variable name indicating",
                            " the clusters of time-varying coefs",
                            " must be inlcuded in the dataset")
    warning.text.c = paste0("variables having",
                            " time-varying coefs",
                            " must be inlcuded in the model")
    
    if(is.null(str.tvarying)){
      stop(warning.text.a)
    }
    if(!is.null(str.tvarying)){
      if(!is.list(str.tvarying)){
        stop(warning.text.a)
      }
      ll.vec = sapply(str.tvarying, length)
      if(!all(ll.vec == 3)){
        stop(warning.text.a)
      }
      check <- sapply(str.tvarying, function(x) x[1])
      check.1 <- sapply(str.tvarying, function(x) x[2])
      check.2 <- sapply(str.tvarying, function(x) x[3])
      if(!all(check %in% pred.names)){
        stop(warning.text.c)
      }
      if(!all(check.1 %in% colnames(data))){
        stop(warning.text.b)
      }
      if(!all(check.2 %in% c("exch", "AR1", "AR2"))){
        stop(warning.text1)
      }
    }
    
    X.design.fix <- model.matrix(formula, data = data)
    var.slop.name = sapply(str.tvarying, function(x) x[1])
    id.slop.name = sapply(str.tvarying, function(x) x[2])
    str.slop.name = sapply(str.tvarying, function(x) x[3])
    
    X.design.rand.slop <- list()
    for(lll in 1:length(var.slop.name)){
      ntimes.slop.lll = length(unique(as.matrix(data[,id.slop.name[lll]])))
      if(ntimes.slop.lll == num.obs){
        X.design.rand.slop[[lll]] <- Diagonal( x = rep(1, num.obs))
      }else if(ntimes.slop.lll < num.obs){
        comp.vec <- paste0("(0+", var.slop.name[lll] ,"|", id.slop.name[lll], ")")
        formula.lll <- as.formula(paste(c(deparse(formula), comp.vec), 
                                        collapse = "+"))
        X.inter <- lme4::lFormula(formula = formula.lll, data = data)
        X.design.rand.slop[[lll]] <- t(X.inter$reTrms$Zt)
      }
    }
    names(X.design.rand.slop) <- var.slop.name
    var.fix.name = colnames(X.design.fix)
    
    re0 <- fit.NB.st.add.s3(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix, 
                            X.design.rand.slop = X.design.rand.slop,
                            offset.vec = offset.vec,
                            rand.slop.str = str.slop.name,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]], 
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  }
  
  
  ## SCENARIO 3: random intercept and no time-varying coefs
  if(rand.int == TRUE & slope.tvarying == FALSE){
    
    warning.text = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of time-specific random effects",
                          " and the structure of time-specific random effects")
    warning.text1 = paste0("the structure can only be specified as ",
                           "exch, AR1, or AR2")
    warning.text2 = paste0("the variable name indicating",
                           " the clusters of time-specific random effects",
                           " must be inlcuded in the dataset")
    
    if(is.null(str.int)){
      stop(warning.text)
    }
    if(!is.null(str.int)){
      if(!is.list(str.int)){
        stop(warning.text)
      }
      ll = length(str.int)
      if(ll != 2){
        stop(warning.text)
      }
      if(ll == 2){
        if(!(str.int[[2]] %in% c("exch", "AR1", "AR2"))){
          stop(warning.text1)
        }
        if(!(str.int[[1]] %in% colnames(data))){
          stop(warning.text2)
        }
      }
    }
    id.int.name = str.int[[1]]
    ntimes.int = length(unique(as.matrix(data[,id.int.name])))
    if(ntimes.int == num.obs){
      X.design.fix <- model.matrix(formula, data = data)
      X.design.rand.int <- Diagonal( x = rep(1, num.obs))
    }else if(ntimes.int < num.obs){
      comp.vec <- paste0("(1|", id.int.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec), 
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.rand.int <- t(X.all$reTrms$Zt)
      X.design.fix <- X.all$X
    }
    
    re0 <- fit.NB.st.add.s2(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix, 
                            X.design.rand.int = X.design.rand.int,
                            offset.vec = offset.vec,
                            rand.int.str = str.int[[2]],
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]], 
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
    
  } # END SCENARIO 2
  
  
  ###### SCENARIO 4: random intercept and random slopes ######
  if(rand.int == TRUE & slope.tvarying == TRUE){
    
    warning.text = paste0("please provide a list containing",
                          " a variable name indicating",
                          " the clusters of time-specific random effects",
                          " and the structure of time-specific random effects")
    warning.text1 = paste0("the structure can only be specified as ",
                           "exch, AR1, or AR2")
    warning.text2 = paste0("the variable name indicating",
                           " the clusters of time-specific random effects",
                           " must be inlcuded in the dataset")
    if(is.null(str.int)){
      stop(warning.text)
    }
    if(!is.null(str.int)){
      if(!is.list(str.int)){
        stop(warning.text)
      }
      ll = length(str.int)
      if(ll != 2){
        stop(warning.text)
      }
      if(ll == 2){
        if(!(str.int[[2]] %in% c("exch", "AR1", "AR2"))){
          stop(warning.text1)
        }
        if(!(str.int[[1]] %in% colnames(data))){
          stop(warning.text2)
        }
      }
    }
    
    warning.text.a = paste0("please provide a list containing",
                            " vectors that specify: ",
                            " (1) variable names of time-varying coefs,",
                            " (2) variable names indicating",
                            " the clusters of time-varying coefs",
                            " , and (3) the structure of time-varying coefs")
    warning.text.b = paste0("the variable name indicating",
                            " the clusters of time-varying coefs",
                            " must be inlcuded in the dataset")
    warning.text.c = paste0("variables having",
                            " time-varying coefs",
                            " must be inlcuded in the model")
    
    if(is.null(str.tvarying)){
      stop(warning.text.a)
    }
    if(!is.null(str.tvarying)){
      if(!is.list(str.tvarying)){
        stop(warning.text.a)
      }
      ll.vec = sapply(str.tvarying, length)
      if(!all(ll.vec == 3)){
        stop(warning.text.a)
      }
      check <- sapply(str.tvarying, function(x) x[1])
      check.1 <- sapply(str.tvarying, function(x) x[2])
      check.2 <- sapply(str.tvarying, function(x) x[3])
      if(!all(check %in% pred.names)){
        stop(warning.text.c)
      }
      if(!all(check.1 %in% colnames(data))){
        stop(warning.text.b)
      }
      if(!all(check.2 %in% c("exch", "AR1", "AR2"))){
        stop(warning.text1)
      }
    }
    
    warning.text.a = paste0("please provide a list containing",
                            " a variable name indicating",
                            " the clusters of location-specific random effects",
                            " and the structure of location-specific random effects")
    warning.text1.a = paste0("the structure can only be specified as ",
                             "exch, or CAR")
    warning.text2.a = paste0("the variable name indicating",
                             " the clusters of location-specific random effects",
                             " must be inlcuded in the dataset")
    
    if(is.null(spatial.str)){
      stop(warning.text.a)
    }
    if(!is.null(spatial.str)){
      if(!is.list(spatial.str)){
        stop(warning.text.a)
      }
      ll = length(spatial.str)
      if(ll != 2){
        stop(warning.text.a)
      }
      if(ll == 2){
        if(!(spatial.str[[2]] %in% c("exch", "CAR"))){
          stop(warning.text1.a)
        }
        if(!(spatial.str[[1]] %in% colnames(data))){
          stop(warning.text2.a)
        }
      }
    }
    
    id.int.name = str.int[[1]]
    ntimes.int = length(unique(as.matrix(data[,id.int.name])))
    if(ntimes.int == num.obs){
      X.design.fix <- model.matrix(formula, data = data)
      X.design.rand.int <- Diagonal( x = rep(1, num.obs))
    }else if(ntimes.int < num.obs){
      comp.vec <- paste0("(1|", id.int.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec), 
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.rand.int <- t(X.all$reTrms$Zt)
      X.design.fix <- X.all$X
    }
    
    var.slop.name = sapply(str.tvarying, function(x) x[1])
    id.slop.name = sapply(str.tvarying, function(x) x[2])
    str.slop.name = sapply(str.tvarying, function(x) x[3])
    
    X.design.rand.slop <- list()
    for(lll in 1:length(var.slop.name)){
      ntimes.slop.lll = length(unique(as.matrix(data[,id.slop.name[lll]])))
      if(ntimes.slop.lll == num.obs){
        X.design.rand.slop[[lll]] <- Diagonal( x = rep(1, num.obs))
      }else if(ntimes.slop.lll < num.obs){
        comp.vec <- paste0("(0+", var.slop.name[lll] ,"|", id.slop.name[lll], ")")
        formula.lll <- as.formula(paste(c(deparse(formula), comp.vec), 
                                        collapse = "+"))
        X.inter <- lme4::lFormula(formula = formula.lll, data = data)
        X.design.rand.slop[[lll]] <- t(X.inter$reTrms$Zt)
      }
    }
    names(X.design.rand.slop) <- var.slop.name
    var.fix.name = colnames(X.design.fix)
    
    id.loc.name = spatial.str[[1]]
    nlocs = length(unique(as.matrix(data[,id.loc.name])))
    if(nlocs == num.obs){
      X.design.spatial <- Diagonal( x = rep(1, num.obs))
    }else if(nlocs < num.obs){
      comp.vec <- paste0("(1|", id.loc.name, ")")
      formula.final <- as.formula(paste(c(deparse(formula), comp.vec), 
                                        collapse = "+"))
      X.all <- lme4::lFormula(formula = formula.final, data = data)
      X.design.spatial <- t(X.all$reTrms$Zt)
    }
    
    re0 <- fit.NB.st.add.s1(outcome.vec = outcome.vec,
                            X.design.fix = X.design.fix, 
                            X.design.rand.int = X.design.rand.int,
                            X.design.rand.slop = X.design.rand.slop,
                            offset.vec = offset.vec,
                            rand.int.str = str.int[[2]],
                            rand.slop.str = str.slop.name,
                            niter = niter, burn_in = burn_in,
                            covariate.AD = covariate.AD,
                            countsonly = countsonly,
                            X.design.spatial = X.design.spatial,
                            spatial.str = spatial.str[[2]], 
                            W.s = adj.mat.s,
                            ID.spacetime.df = ID.spacetime.df)
  } # END SCENARIO 4 (space-time additive)
  
  return(re0)
  ### OUTPUTS: a list contains:
  ## (1) "beta": posterior samples of fixed effects or 
  ###    mean of time-varying coefficients
  ## (2) "parm": posterior samples of parameters controlling temporal dependency
  ###  of time-specific random effects, 
  ###   parameters controlling spatial dependency of location-specific random effects,
  ###   and the over-dispersion parameter
  ## (3) "rand.int.t": posterior samples of time-specific random effects 
  ###     (when countsonly = FALSE)
  ## (4) "rand.int.s":  posterior samples of location-specific effects 
  ###     (when countsonly = FALSE)   
  ## (5) "beta.tvarying": a list contains posterior samples of
  ### time-varying coefs (when slope.tvarying = TRUE and countsonly = FALSE)
  ## (6) "parm.beta.tvarying": a list contains posterior samples of
  ### parameters controlling temporal dependency
  ###  of time-varying coefs (when slope.tvarying = TRUE and countsonly = FALSE)
  ## (7) "pred.counts": posterior samples (burn-in samples are discarded) of
  ### predicted outcome counts
  ## (8) "pred.null.counts": posterior samples (burn-in samples are discarded) of
  ### predicted deaths assuming no certain infections
  ## (9) "AD": posterior samples (burn-in samples are discarded) of
  ### attributable deaths 
  ## (10): "accept.rate": acceptance rate of M-H algorithm applied for updating 
  ### the over-dispersion parameter 
  ## (11): "WAIC": the WAIC measuring how the model fits to the data
}



### space-time additive model
### time-varying coefficients and time-specific random effects ###
fit.NB.st.add.s1 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.int,
                             X.design.rand.slop,
                             offset.vec = NULL,
                             rand.int.str,
                             rand.slop.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL, 
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  beta.fix.names <- colnames(X.design.fix)
  
  outcome.vec <- as.numeric(outcome.vec)
  
  nbeta.rand.vec <- sapply(X.design.rand.slop, ncol)
  beta.rand.names <- names(X.design.rand.slop)
  
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  beta.fix.names.0 <- beta.fix.names[! beta.fix.names %in% beta.rand.names]
  nbeta.fix.0 <- length(beta.fix.names.0)
  if(nbeta.fix.0 == 1){
    c.beta0 = matrix(0, ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix.0), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix.0))
  }
  
  C.beta.bar <- matrix(1/100, ncol = 1)
  c.beta.bar <- matrix(0, ncol = 1)
  
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  colnames(beta.mat) <- beta.fix.names
  beta.mat[1, ] <- rep(0, nbeta.fix)
  
  ntimes.int <- ncol(X.design.rand.int)
  e.mat <- matrix(NA, ncol = ntimes.int, nrow = niter + 1)
  e.mat[1, ] <- rnorm(ntimes.int, 0, 1)
  
  ntimes.slop.vec <- as.numeric(nbeta.rand.vec)
  beta.rand.list <- list()
  parm.beta.list <- list()
  W.t.slop <- Dw.t.slop <- lambda.t.slop <- ll.rho.slop <- list()
  nvar.rand <- length(nbeta.rand.vec)
  for(lll in 1:nvar.rand){
    inter.lll <- matrix(NA, nrow = niter + 1,
                        ncol = nbeta.rand.vec[lll])
    inter.lll[1, ] <- rnorm(nbeta.rand.vec[lll], 0, 1)
    beta.rand.list[[lll]] <- inter.lll
    
    if(rand.slop.str[lll] == "exch"){
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 1, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0("sigma2.", 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- 1
    }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
      
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 2, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0(c("sigma2.", "rho."), 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- c(1, 0.5)
      
      W.t.slop[[lll]] = compute.adjmat(ntimes = ntimes.slop.vec[lll], 
                                       str = rand.slop.str[lll])
      Dw.t.slop[[lll]] <- Diagonal(x = apply(W.t.slop[[lll]], 1, sum))
      lambda.t.slop[[lll]] = eigen(solve(Dw.t.slop[[lll]])%*%W.t.slop[[lll]], 
                                   only.values = TRUE)$values
      rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
      ll.rho.slop[[lll]] = sapply(rho.prior.val, 
                                  function(x) 0.5*sum(log(1-x*lambda.t.slop[[lll]])), 
                                  simplify = TRUE)
      
    }
  }
  names(beta.rand.list) <- beta.rand.names
  names(parm.beta.list) <- beta.rand.names
  # names(W.t.slop) <- names(Dw.t.slop) <- beta.rand.names
  # names(lambda.t.slop) <- names(ll.rho.slop) <- beta.rand.names
  
  if(rand.int.str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "xi")
    parm.mat[1, ] <- c(1, round(mean(outcome.vec)))
  }else if(rand.int.str %in% c("AR1", "AR2")){
    
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "rho.t", "xi")
    parm.mat[1, ] <- c(1, 0.5, round(mean(outcome.vec)))
    
    W.t.int = compute.adjmat(ntimes = ntimes.int, str = rand.int.str)
    Dw.t.int <- Diagonal(x = apply(W.t.int, 1, sum))
    lambda.t.int = eigen(solve(Dw.t.int)%*%W.t.int, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.int.1 = sapply(rho.prior.val, 
                          function(x) 0.5*sum(log(1-x*lambda.t.int)), simplify = TRUE)
    
  }
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  X.design.fix.0 <- matrix(X.design.fix[,beta.fix.names.0], ncol = nbeta.fix.0)
  colnames(X.design.fix.0) <- beta.fix.names.0
  
  pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
  
  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[1, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.int = as.numeric(X.design.rand.int %*%
                               matrix(e.mat[1, ], ncol = 1))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][1, ], ncol = 1)))
  comp.spatial = as.numeric(X.design.spatial %*%
                              matrix(theta.mat[1, ], ncol = 1))
  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + comp.fix + comp.rand.int + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }else{
    pg.parm2 <- comp.fix + comp.rand.int + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }
  pg.parm2 <- as.numeric(pg.parm2)
  ntotal = num.obs
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2.1 update constant beta coefficients 
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix.0*sqrt(omega.vec)) + C.beta0.pre)
    e.vec.i <- as.numeric(X.design.rand.int %*% matrix(e.mat[i, ], ncol = 1))
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    cov.rand.list.i <- get.randcovariate.list(X.design.rand.slop = X.design.rand.slop,
                                              beta.rand.list = beta.rand.list,
                                              curr.index = i)
    cov.rand.vec.i <- Reduce("+", cov.rand.list.i)
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - cov.rand.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - cov.rand.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix.0)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, beta.fix.names.0] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff.0 <- as.numeric(X.design.fix.0 %*% 
                              matrix(beta.mat[i+1, beta.fix.names.0], ncol = 1))
    
    # 2.2 update time-varying beta coefficients  
    # 2.3 update mean of time-varying coefficients 
    # 2.4 update parameters controlling smoothness of random slopes
    update.beta <- function(X.design.beta, C.pre, c.pri, 
                            omega.vec, res.vec.i){
      M.beta <- solve(crossprod(X.design.beta*sqrt(omega.vec)) + C.pre)
      m.beta <- M.beta%*%(C.pre%*%c.pri + 
                            t(sqrt(omega.vec)*X.design.beta)%*%
                            (sqrt(omega.vec)*res.vec.i))
      return(list(M = M.beta, m = m.beta))
    }
    
    if(!is.null(offset.vec)){
      res.beta.i0 <- z.vec - fix.eff.0 - e.vec.i - theta.vec.i - offset.vec
    }else{
      res.beta.i0 <- z.vec - fix.eff.0 - e.vec.i - theta.vec.i
    }
    
    for(lll in 1:nvar.rand){
      
      ## 2.2 update time-specific coefficients
      if(rand.slop.str[lll] == "exch"){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- Diagonal(x = rep(1/parm.beta.lll[1], 
                                           ntimes.slop.vec[lll]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- as.matrix((1/parm.beta.lll[1])*
                                      (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }
      
      if(lll == 1 | lll == nvar.rand){
        
        if(nvar.rand == 1){
          res.beta.i <- res.beta.i0
        }else{
          index.iter <- ifelse(lll==1, i, i+1)
          cov.rand.vec.lll = 
            get.randcovariate.list(X.design.rand.slop = 
                                     X.design.rand.slop[-lll],
                                   beta.rand.list = beta.rand.list[-lll],
                                   curr.index = index.iter)
          res.beta.i <- res.beta.i0 - Reduce("+", cov.rand.vec.lll)
        }
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }else{
        
        cov.rand.vec.lll.0 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[(1:(lll-1))],
                                 beta.rand.list = beta.rand.list[(1:(lll-1))],
                                 curr.index = i+1)
        cov.rand.vec.lll.1 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[((lll+1):nvar.rand)],
                                 beta.rand.list = 
                                   beta.rand.list[((lll+1):nvar.rand)],
                                 curr.index = i)
        res.beta.i <- res.beta.i0 - 
          Reduce("+", cov.rand.vec.lll.0) - 
          Reduce("+", cov.rand.vec.lll.1)
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }
      
      ## 2.3 update mean of time-varying coefficients
      X.design.betabar = matrix(rep(1, as.numeric(nbeta.rand.vec[lll])), 
                                ncol = 1)
      pre.multi.beta = matrix(colSums(C.beta.pre.lll), nrow = 1)
      M <- solve(pre.multi.beta%*%X.design.betabar + C.beta.bar)
      m <- M%*%(C.beta.bar%*%c.beta.bar + 
                  pre.multi.beta%*%matrix(beta.rand.list[[lll]][i+1, ], 
                                          ncol = 1))
      beta.mat[i+1, beta.rand.names[lll]] <- 
        as.numeric( mvrnorm(1, mu = m, Sigma = M) )
      
      ## 2.4 update parameter controlling smoothness of random slopes
      ncluster = as.numeric(nbeta.rand.vec[lll])
      if(rand.slop.str[lll] == "exch"){
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, b1 + as.numeric(sum(res.vec.i^2))/2)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, 
                   b1 + as.numeric(t(matrix(res.vec.i, ncol = 1))%*%
                                     (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]])%*%
                                     matrix(res.vec.i, ncol = 1))/2)
        
        # paste0("rho.", beta.rand.names[lll]) (M-H)
        inter = as.numeric( t(matrix(res.vec.i, ncol = 1))%*%
                              W.t.slop[[lll]]%*%matrix(res.vec.i, ncol = 1) )
        ll.rho = ll.rho.slop[[lll]] + rho.prior.val/(2*parm.beta.list[[lll]][i+1, 1]^2)*inter
        parm.beta.list[[lll]][i+1, 2] <- sample(x = rho.prior.val, size = 1, 
                                                prob = exp(ll.rho - max(ll.rho)))
        
      } 
      ## end 2.4
    } # end loop over time-varying coefficients 
    
    # 3. update e (temporal random effects; random intercept) 
    ## posterior mean and variance of e ##
    rand.eff.i <- get.randcovariate.list(X.design.rand.slop = 
                                           X.design.rand.slop,
                                         beta.rand.list = beta.rand.list,
                                         curr.index = i+1)
    rand.eff.i <- Reduce("+", rand.eff.i)
    
    if(rand.int.str == "exch"){
      pre.t.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.t"]), ntimes.int))
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - theta.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%
        (t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, 
                                          Sigma = Sigma.t.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.t"] <- 1/rgamma(1, ntimes.int/2 + a1, 
                                          b1 + as.numeric(sum(e.mat[i+1, ]^2))/2)
      
    }else if(rand.int.str %in% c("AR1", "AR2")){
      
      pre.t.0 <- (1/parm.mat[i, "tau2.t"])*(Dw.t.int - parm.mat[i, "rho.t"]*W.t.int)
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - theta.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.t"] <- 
        1/rgamma(1, ntimes.int/2 + a1, 
                 b1 + as.numeric(t(matrix(e.mat[i+1, ], ncol = 1))%*%
                                   (Dw.t.int - parm.mat[i, "rho.t"]*W.t.int)%*%
                                   matrix(e.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(e.mat[i+1, ], ncol = 1))%*%
                            W.t.int%*%matrix(e.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.int.1 + rho.prior.val/(2*parm.mat[i+1, "tau2.t"]^2)*inter
      parm.mat[i+1, "rho.t"] <- sample(x = rho.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    e.vec.i <- as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - e.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - e.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + fix.eff.0 + rand.eff.i + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = fix.eff.0 + rand.eff.i + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }
    q.nb = 1/(1+exp(eta.vec)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, 
                            outcome.vec + parm.mat[i+1, "xi"], eta.vec)
    
    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
  } # end MCMC
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                          do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                          e.mat[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix.0, 
                         do.call(cbind, X.design.rand.slop),
                         X.design.rand.int,
                         X.design.spatial)
  X.design.comb <- as.matrix(X.design.comb)
  pred.val.mat <- get_predvals(parm_post = parm.post.comb,
                               X_all = X.design.comb,
                               ntotal = ntotal)
  
  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)), 
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1, 
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  if(all(covariate.AD %in% names(X.design.rand.slop))){
    index.rand <- !(names(X.design.rand.slop) %in% covariate.AD)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list[index.rand])[-c(1:burn_in),],
                            e.mat[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- cbind(X.design.fix.0, 
                           do.call(cbind, X.design.rand.slop[index.rand]),
                           X.design.rand.int,
                           X.design.spatial)
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(any(covariate.AD %in% names(X.design.rand.slop))){
    index.curr <- covariate.AD %in% names(X.design.rand.slop)
    index.fix <- covariate.AD[!index.curr]
    index.rand <- !names(X.design.rand.slop)%in%covariate.AD
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,index.fix] <- 0
    beta.rand.comb <- beta.rand.list[index.rand]
    if(length(beta.rand.comb) == 0){
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              e.mat[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             X.design.rand.int,
                             X.design.spatial)
    }else{
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              do.call(cbind, beta.rand.comb)[-c(1:burn_in), ],
                              e.mat[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             do.call(cbind, 
                                     X.design.rand.slop[index.rand]),
                             X.design.rand.int,
                             X.design.spatial)
    }
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(!any(covariate.AD %in% names(X.design.rand.slop))){
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,covariate.AD] <- 0
    X.design.null <- cbind(X.design.fix.null, 
                           do.call(cbind, X.design.rand.slop),
                           X.design.rand.int,
                           X.design.spatial)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                            e.mat[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }
  
  if(!is.null(offset.vec)){
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }else{
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }
  AD.mat <- pred.counts.mat - pred.counts.null.mat
  
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.mat <- data.frame(ID.spacetime.df, pred.counts.mat)
    pred.counts.null.mat <- data.frame(ID.spacetime.df, pred.counts.null.mat)
  }
  
  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat, 
               parm.beta.tvarying = parm.beta.list, 
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.t = e.mat, 
               rand.int.s = theta.mat,
               beta.tvarying = beta.rand.list,
               parm.beta.tvarying = parm.beta.list,
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  return(re)
  
  
}


### space-time additive model
### random intercept (baseline temporal trend) and constant slopes ###
fit.NB.st.add.s2 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.int,
                             offset.vec = NULL,
                             rand.int.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  outcome.vec <- as.numeric(outcome.vec)
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  if(nbeta.fix == 1){
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix))
  }
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ntimes = ncol(X.design.rand.int)
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  e.mat <- matrix(NA, ncol = ntimes, nrow = niter + 1)
  beta.mat[1, ] <- rep(0, nbeta.fix)
  e.mat[1, ] <- rnorm(ntimes, 0, 1)
  
  if(rand.int.str == "exch"){
    parm.mat <- as.data.frame(matrix(NA, ncol = 2, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "xi")
    parm.mat[1, ] <- c(1, round(mean(outcome.vec)))
  }else if(rand.int.str %in% c("AR1", "AR2")){
    parm.mat <- as.data.frame(matrix(NA, ncol = 3, nrow = niter + 1))
    colnames(parm.mat) <- c("tau2.t", "rho.t", "xi")
    parm.mat[1, ] <- c(1, 0.5, round(mean(outcome.vec)))
  }
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  if(!is.null(offset.vec)){
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.rand.int %*% matrix(e.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }
  
  if(rand.int.str == "AR1"){
    adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))
    
    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.1 = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)), simplify = TRUE)
  }else if(rand.int.str == "AR2"){
    adj.list <- as.matrix(rbind(cbind(c(1:(ntimes-1)), c(2:ntimes)),
                                cbind(c(1:(ntimes-2)), c(3:ntimes))))
    adj.list <- adj.list[order(adj.list[,1]), ]
    W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
    Dw.t <- Diagonal(x = apply(W.t, 1, sum))
    
    lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
    rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.1 = sapply(rho.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)), simplify = TRUE)
  }
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix*sqrt(omega.vec)) + C.beta0.pre)
    e.vec.i <- as.numeric(X.design.rand.int %*% matrix(e.mat[i, ], ncol = 1))
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - e.vec.i - theta.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - e.vec.i - theta.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    # 3. update e (temporal random effects) 
    ## posterior mean and variance of e ##
    fix.eff <- X.design.fix %*% matrix(beta.mat[i+1, ], ncol = 1)
    if(rand.int.str == "exch"){
      pre.t.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.t"]), ntimes))
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - theta.vec.i -  offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.t"] <- 1/rgamma(1, ntimes/2 + a1, 
                                          b1 + as.numeric(sum(e.mat[i+1, ]^2))/2)
      
    }else if(rand.int.str %in% c("AR1", "AR2")){
      
      pre.t.0 <- (1/parm.mat[i, "tau2.t"])*(Dw.t - parm.mat[i, "rho.t"]*W.t)
      pre.t.post <- t(X.design.rand.int)%*%Omega%*%X.design.rand.int + pre.t.0
      Sigma.t.post <- solve(pre.t.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - theta.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - theta.vec.i
      }
      mu.t.post <- Sigma.t.post%*%(t(X.design.rand.int)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      e.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.t.post, Sigma = Sigma.t.post) )
      
      # 4. update tau22 
      parm.mat[i+1, "tau2.t"] <- 
        1/rgamma(1, ntimes/2 + a1, 
                 b1 + as.numeric(t(matrix(e.mat[i+1, ], ncol = 1))%*%
                                   (Dw.t - parm.mat[i, "rho.t"]*W.t)%*%
                                   matrix(e.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(e.mat[i+1, ], ncol = 1))%*%
                            W.t%*%matrix(e.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.1 + rho.prior.val/(2*parm.mat[i+1, "tau2.t"]^2)*inter
      parm.mat[i+1, "rho.t"] <- sample(x = rho.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
      
    }
    
    
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    e.vec.i <- as.numeric(X.design.rand.int%*%e.mat[i+1, ], ncol = 1)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - e.vec.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - e.vec.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff  - e.vec.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 9. update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        e.vec.i + as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }
    q.nb = 1/(1+exp(eta.vec)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, 
                            outcome.vec + parm.mat[i+1, "xi"], eta.vec)
    
    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
    
  } # end MCMC
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), ],
                          e.mat[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix, 
                         X.design.rand.int,
                         X.design.spatial)
  X.design.comb <- as.matrix(X.design.comb)
  pred.val.mat <- get_predvals(parm_post = parm.post.comb,
                               X_all = X.design.comb,
                               ntotal = ntotal)
  
  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)), 
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1, 
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  X.design.fix.null <- X.design.fix
  X.design.fix.null[,covariate.AD] <- 0
  X.design.null <- cbind(X.design.fix.null, 
                         X.design.rand.int,
                         X.design.spatial)
  X.design.null <- as.matrix(X.design.null)
  pred.val.null.mat <- get_predvals(parm_post = parm.post.comb,
                                    X_all = X.design.null,
                                    ntotal = ntotal)
  
  if(!is.null(offset.vec)){
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }else{
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }
  AD.mat <- pred.counts.mat - pred.counts.null.mat
  
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.mat <- data.frame(ID.spacetime.df, pred.counts.mat)
    pred.counts.null.mat <- data.frame(ID.spacetime.df, pred.counts.null.mat)
  }
  
  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat, 
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.t = e.mat,
               rand.int.s = theta.mat,
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  
  return(re)
  
}




### time-varying coefficients and NO time-specific random effects ###
fit.NB.st.add.s3 <- function(outcome.vec,
                             X.design.fix, 
                             X.design.rand.slop,
                             offset.vec = NULL,
                             rand.slop.str,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  beta.fix.names <- colnames(X.design.fix)
  
  outcome.vec <- as.numeric(outcome.vec)
  
  nbeta.rand.vec <- sapply(X.design.rand.slop, ncol)
  beta.rand.names <- names(X.design.rand.slop)
  
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  beta.fix.names.0 <- beta.fix.names[! beta.fix.names %in% beta.rand.names]
  nbeta.fix.0 <- length(beta.fix.names.0)
  if(nbeta.fix.0 == 1){
    c.beta0 = matrix(0, ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix.0), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix.0))
  }
  
  C.beta.bar <- matrix(1/100, ncol = 1)
  c.beta.bar <- matrix(0, ncol = 1)
  
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  colnames(beta.mat) <- beta.fix.names
  beta.mat[1, ] <- rep(0, nbeta.fix)
  
  ntimes.slop.vec <- as.numeric(nbeta.rand.vec)
  beta.rand.list <- list()
  parm.beta.list <- list()
  W.t.slop <- Dw.t.slop <- lambda.t.slop <- ll.rho.slop <- list()
  nvar.rand <- length(nbeta.rand.vec)
  for(lll in 1:nvar.rand){
    inter.lll <- matrix(NA, nrow = niter + 1,
                        ncol = nbeta.rand.vec[lll])
    inter.lll[1, ] <- rnorm(nbeta.rand.vec[lll], 0, 1)
    beta.rand.list[[lll]] <- inter.lll
    
    if(rand.slop.str[lll] == "exch"){
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 1, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0("sigma2.", 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- 1
    }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
      
      parm.beta.list[[lll]] <- as.data.frame(matrix(NA, ncol = 2, 
                                                    nrow = niter + 1))
      colnames(parm.beta.list[[lll]]) <- c(paste0(c("sigma2.", "rho."), 
                                                  beta.rand.names[lll]))
      parm.beta.list[[lll]][1, ] <- c(1, 0.5)
      
      W.t.slop[[lll]] = compute.adjmat(ntimes = ntimes.slop.vec[lll], 
                                       str = rand.slop.str[lll])
      Dw.t.slop[[lll]] <- Diagonal(x = apply(W.t.slop[[lll]], 1, sum))
      lambda.t.slop[[lll]] = eigen(solve(Dw.t.slop[[lll]])%*%W.t.slop[[lll]], 
                                   only.values = TRUE)$values
      rho.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
      ll.rho.slop[[lll]] = sapply(rho.prior.val, 
                                  function(x) 0.5*sum(log(1-x*lambda.t.slop[[lll]])), 
                                  simplify = TRUE)
      
    }
  }
  names(beta.rand.list) <- beta.rand.names
  names(parm.beta.list) <- beta.rand.names
  # names(W.t.slop) <- names(Dw.t.slop) <- beta.rand.names
  # names(lambda.t.slop) <- names(ll.rho.slop) <- beta.rand.names
  
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  X.design.fix.0 <- matrix(X.design.fix[,beta.fix.names.0], ncol = nbeta.fix.0)
  colnames(X.design.fix.0) <- beta.fix.names.0
  
  pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
  
  comp.fix =  as.numeric(X.design.fix.0 %*%
                           matrix(beta.mat[1, beta.fix.names.0],
                                  ncol = nbeta.fix.0))
  comp.rand.slop.list = lapply(1:nvar.rand,
                               function(x)
                                 as.numeric(X.design.rand.slop[[x]] %*%
                                              matrix(beta.rand.list[[x]][1, ], ncol = 1)))
  comp.spatial = as.numeric(X.design.spatial %*%
                              matrix(theta.mat[1, ], ncol = 1))
  if(!is.null(offset.vec)){
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + comp.fix + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }else{
    pg.parm2 <- comp.fix + 
      Reduce("+", comp.rand.slop.list) + comp.spatial
  }
  pg.parm2 <- as.numeric(pg.parm2)
  ntotal = num.obs
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2.1 update constant beta coefficients 
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix.0*sqrt(omega.vec)) + C.beta0.pre)
    theta.vec.i <- as.numeric(X.design.spatial %*% matrix(theta.mat[i, ], ncol = 1))
    cov.rand.list.i <- get.randcovariate.list(X.design.rand.slop = X.design.rand.slop,
                                              beta.rand.list = beta.rand.list,
                                              curr.index = i)
    cov.rand.vec.i <- Reduce("+", cov.rand.list.i)
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec - theta.vec.i - cov.rand.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec - theta.vec.i - cov.rand.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix.0)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, beta.fix.names.0] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff.0 <- as.numeric(X.design.fix.0 %*% 
                              matrix(beta.mat[i+1, beta.fix.names.0], ncol = 1))
    
    # 2.2 update time-varying beta coefficients  
    # 2.3 update mean of time-varying coefficients 
    # 2.4 update parameters controlling smoothness of random slopes
    update.beta <- function(X.design.beta, C.pre, c.pri, 
                            omega.vec, res.vec.i){
      M.beta <- solve(crossprod(X.design.beta*sqrt(omega.vec)) + C.pre)
      m.beta <- M.beta%*%(C.pre%*%c.pri + 
                            t(sqrt(omega.vec)*X.design.beta)%*%
                            (sqrt(omega.vec)*res.vec.i))
      return(list(M = M.beta, m = m.beta))
    }
    
    if(!is.null(offset.vec)){
      res.beta.i0 <- z.vec - fix.eff.0 - theta.vec.i - offset.vec
    }else{
      res.beta.i0 <- z.vec - fix.eff.0 - theta.vec.i
    }
    
    for(lll in 1:nvar.rand){
      
      ## 2.2 update time-specific coefficients
      if(rand.slop.str[lll] == "exch"){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- Diagonal(x = rep(1/parm.beta.lll[1], 
                                           ntimes.slop.vec[lll]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        C.beta.pre.lll <- as.matrix((1/parm.beta.lll[1])*
                                      (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]]))
        c.beta = matrix(rep(beta.mat[i, beta.rand.names[lll]], 
                            nbeta.rand.vec[lll]), 
                        ncol = 1)
      }
      
      if(lll == 1 | lll == nvar.rand){
        
        if(nvar.rand == 1){
          res.beta.i <- res.beta.i0
        }else{
          index.iter <- ifelse(lll==1, i, i+1)
          cov.rand.vec.lll = 
            get.randcovariate.list(X.design.rand.slop = 
                                     X.design.rand.slop[-lll],
                                   beta.rand.list = beta.rand.list[-lll],
                                   curr.index = index.iter)
          res.beta.i <- res.beta.i0 - Reduce("+", cov.rand.vec.lll)
        }
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }else{
        
        cov.rand.vec.lll.0 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[(1:(lll-1))],
                                 beta.rand.list = beta.rand.list[(1:(lll-1))],
                                 curr.index = i+1)
        cov.rand.vec.lll.1 = 
          get.randcovariate.list(X.design.rand.slop = 
                                   X.design.rand.slop[((lll+1):nvar.rand)],
                                 beta.rand.list = 
                                   beta.rand.list[((lll+1):nvar.rand)],
                                 curr.index = i)
        res.beta.i <- res.beta.i0 - 
          Reduce("+", cov.rand.vec.lll.0) - 
          Reduce("+", cov.rand.vec.lll.1)
        res.beta.i <- matrix(res.beta.i, ncol = 1)
        
        mM.list <- update.beta(X.design.beta = 
                                 X.design.rand.slop[[lll]], 
                               C.pre = C.beta.pre.lll, 
                               c.pri = c.beta, 
                               omega.vec = omega.vec, 
                               res.vec.i = res.beta.i)
        
        beta.rand.list[[lll]][i+1, ] <- as.numeric( 
          mvrnorm(1, mu = mM.list$m, Sigma = mM.list$M) )
        
      }
      
      ## 2.3 update mean of time-varying coefficients
      X.design.betabar = matrix(rep(1, as.numeric(nbeta.rand.vec[lll])), 
                                ncol = 1)
      pre.multi.beta = matrix(colSums(C.beta.pre.lll), nrow = 1)
      M <- solve(pre.multi.beta%*%X.design.betabar + C.beta.bar)
      m <- M%*%(C.beta.bar%*%c.beta.bar + 
                  pre.multi.beta%*%matrix(beta.rand.list[[lll]][i+1, ], 
                                          ncol = 1))
      beta.mat[i+1, beta.rand.names[lll]] <- 
        as.numeric( mvrnorm(1, mu = m, Sigma = M) )
      
      ## 2.4 update parameter controlling smoothness of random slopes
      ncluster = as.numeric(nbeta.rand.vec[lll])
      if(rand.slop.str[lll] == "exch"){
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, b1 + as.numeric(sum(res.vec.i^2))/2)
      }else if(rand.slop.str[lll] %in% c("AR1", "AR2")){
        
        parm.beta.lll = as.numeric(parm.beta.list[[lll]][i, ])
        res.vec.i <- as.numeric(beta.rand.list[[lll]][i+1, ] - 
                                  beta.mat[i+1, beta.rand.names[lll]])
        
        # paste0("sigma2.", beta.rand.names[lll])
        parm.beta.list[[lll]][i+1, 1] <- 
          1/rgamma(1, ncluster/2 + a1, 
                   b1 + as.numeric(t(matrix(res.vec.i, ncol = 1))%*%
                                     (Dw.t.slop[[lll]] - parm.beta.lll[2]*W.t.slop[[lll]])%*%
                                     matrix(res.vec.i, ncol = 1))/2)
        
        # paste0("rho.", beta.rand.names[lll]) (M-H)
        inter = as.numeric( t(matrix(res.vec.i, ncol = 1))%*%
                              W.t.slop[[lll]]%*%matrix(res.vec.i, ncol = 1) )
        ll.rho = ll.rho.slop[[lll]] + rho.prior.val/(2*parm.beta.list[[lll]][i+1, 1]^2)*inter
        parm.beta.list[[lll]][i+1, 2] <- sample(x = rho.prior.val, size = 1, 
                                                prob = exp(ll.rho - max(ll.rho)))
        
      } 
      ## end 2.4
    } # end loop over time-varying coefficients 
    
    
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    rand.eff.i <- get.randcovariate.list(X.design.rand.slop = 
                                           X.design.rand.slop,
                                         beta.rand.list = beta.rand.list,
                                         curr.index = i+1)
    rand.eff.i <- Reduce("+", rand.eff.i)
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff.0 - rand.eff.i
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff.0  - rand.eff.i - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff.0  - rand.eff.i
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + fix.eff.0 + rand.eff.i + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = fix.eff.0 + rand.eff.i + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }
    q.nb = 1/(1+exp(eta.vec)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, 
                            outcome.vec + parm.mat[i+1, "xi"], eta.vec)
    
    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
  } # end MCMC
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                          do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix.0, 
                         do.call(cbind, X.design.rand.slop),
                         X.design.spatial)
  X.design.comb <- as.matrix(X.design.comb)
  pred.val.mat <- get_predvals(parm_post = parm.post.comb,
                               X_all = X.design.comb,
                               ntotal = ntotal)
  
  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)), 
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1, 
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  if(all(covariate.AD %in% names(X.design.rand.slop))){
    index.rand <- !(names(X.design.rand.slop) %in% covariate.AD)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list[index.rand])[-c(1:burn_in),],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- cbind(X.design.fix.0, 
                           do.call(cbind, X.design.rand.slop[index.rand]),
                           X.design.spatial)
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(any(covariate.AD %in% names(X.design.rand.slop))){
    index.curr <- covariate.AD %in% names(X.design.rand.slop)
    index.fix <- covariate.AD[!index.curr]
    index.rand <- !names(X.design.rand.slop)%in%covariate.AD
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,index.fix] <- 0
    beta.rand.comb <- beta.rand.list[index.rand]
    if(length(beta.rand.comb) == 0){
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             X.design.spatial)
    }else{
      parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                              do.call(cbind, beta.rand.comb)[-c(1:burn_in), ],
                              theta.mat[-c(1:burn_in), ])
      X.design.null <- cbind(X.design.fix.0,
                             do.call(cbind, 
                                     X.design.rand.slop[index.rand]),
                             X.design.spatial)
    }
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }else if(!any(covariate.AD %in% names(X.design.rand.slop))){
    X.design.fix.null <- X.design.fix.0
    X.design.fix.null[,covariate.AD] <- 0
    X.design.null <- cbind(X.design.fix.null, 
                           do.call(cbind, X.design.rand.slop),
                           X.design.spatial)
    parm.post.null <- cbind(beta.mat[-c(1:burn_in), beta.fix.names.0],
                            do.call(cbind, beta.rand.list)[-c(1:burn_in), ],
                            theta.mat[-c(1:burn_in), ])
    X.design.null <- as.matrix(X.design.null)
    pred.val.null.mat <- get_predvals(parm_post = parm.post.null,
                                      X_all = X.design.null,
                                      ntotal = ntotal)
  }
  
  if(!is.null(offset.vec)){
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }else{
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }
  AD.mat <- pred.counts.mat - pred.counts.null.mat
  
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.mat <- data.frame(ID.spacetime.df, pred.counts.mat)
    pred.counts.null.mat <- data.frame(ID.spacetime.df, pred.counts.null.mat)
  }
  
  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat, 
               parm.beta.tvarying = parm.beta.list, 
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.s = theta.mat,
               beta.tvarying = beta.rand.list,
               parm.beta.tvarying = parm.beta.list,
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat,
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  return(re)
  
  
}

### space-time additive model
### no random intercept (baseline temporal trend) and constant slopes ###
fit.NB.st.add.s4 <- function(outcome.vec,
                             X.design.fix, 
                             offset.vec = NULL,
                             niter, burn_in,
                             covariate.AD,
                             countsonly,
                             X.design.spatial,
                             spatial.str, W.s = NULL,
                             ID.spacetime.df){
  ## INPUTs:
  ## niter (# MCMC iterations), nstates (# locations), nweeks (# time points)
  
  num.obs <- length(outcome.vec)
  nbeta.fix <- ncol(X.design.fix)
  outcome.vec <- as.numeric(outcome.vec)
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  if(nbeta.fix == 1){
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = matrix(1/100, ncol = 1)
  }else{
    c.beta0 = matrix(rep(0, nbeta.fix), ncol = 1)
    C.beta0.pre = diag(rep(1/100, nbeta.fix))
  }
  ## tau2.t, inverse-Gamma 
  a1 = 0.1; b1 = 0.1
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = num.obs, nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.fix, nrow = niter + 1)
  beta.mat[1, ] <- rep(0, nbeta.fix)
  
  parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
  colnames(parm.mat) <- c("xi")
  parm.mat[1, ] <- round(mean(outcome.vec))
  
  nlocs = ncol(X.design.spatial)
  theta.mat <- matrix(NA, ncol = nlocs, nrow = niter + 1)
  theta.mat[1, ] <- rnorm(nlocs, 0, 1)
  if(spatial.str == "exch"){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
  }else if(spatial.str != "exch" & !is.null(W.s)){
    parm.mat[["tau2.s"]] <- c(1, rep(NA, niter))
    parm.mat[["rho.s"]] <- c(0.5, rep(NA, niter))
    Dw.s <- Diagonal(x = apply(W.s, 1, sum))
    
    lambda.s = eigen(solve(Dw.s)%*%W.s, only.values = TRUE)$values
    rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
    ll.rho.s.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.s)), simplify = TRUE)
    
  }
  
  ## initial values ##
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  
  if(!is.null(offset.vec)){
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    offset.vec <- as.numeric(offset.vec)
    pg.parm2 <- offset.vec + X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) + 
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }else{
    pg.parm1 <- outcome.vec + parm.mat[1, "xi"]
    pg.parm2 <- X.design.fix%*%matrix(beta.mat[1, ], ncol = 1) +
      X.design.spatial %*% matrix(theta.mat[1, ], ncol = 1)
    pg.parm2 <- as.numeric(pg.parm2)
    ntotal = num.obs
    omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  }
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    # 2. update beta (conjugate)
    ## compute latent normal variables ##
    z.vec <- (outcome.vec - parm.mat[i, "xi"])/(2*omega.vec)
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.fix*sqrt(omega.vec)) + C.beta0.pre)
    theta.vec.i <- as.numeric(X.design.spatial %*% 
                                matrix(theta.mat[i, ], ncol = 1))
    if(!is.null(offset.vec)){
      res.vec.i <- z.vec  - theta.vec.i - offset.vec
    }else{
      res.vec.i <- z.vec  - theta.vec.i
    }
    m <- M%*%(C.beta0.pre%*%c.beta0 + 
                t(sqrt(omega.vec)*X.design.fix)%*%
                (sqrt(omega.vec)*matrix(res.vec.i, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    fix.eff <- X.design.fix %*% matrix(beta.mat[i+1, ], ncol = 1)
    # 4. update theta (spatial random effects) 
    ## posterior mean and variance of e ##
    if(spatial.str == "exch"){
      pre.s.0 <- Diagonal(x = rep((1/parm.mat[i, "tau2.s"]), nlocs))
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff  - offset.vec
      }else{
        res.vec.i <- z.vec - fix.eff 
      }
      mu.s.post <- Sigma.s.post%*%
        (t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 1/rgamma(1, nlocs/2 + a1, 
                                          b1 + as.numeric(sum(theta.mat[i+1, ]^2))/2)
      
    }else if(spatial.str == "CAR"){
      
      pre.s.0 <- (1/parm.mat[i, "tau2.s"])*(Dw.s - parm.mat[i, "rho.s"]*W.s)
      pre.s.post <- t(X.design.spatial)%*%Omega%*%X.design.spatial + pre.s.0
      Sigma.s.post <- solve(pre.s.post)
      if(!is.null(offset.vec)){
        res.vec.i <- z.vec - fix.eff - offset.vec
      }else{
        res.vec.i <-  z.vec - fix.eff 
      }
      mu.s.post <- Sigma.s.post%*%(t(X.design.spatial)%*%Omega%*%matrix(res.vec.i, ncol = 1))
      theta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = mu.s.post, 
                                              Sigma = Sigma.s.post) )
      
      # 4. update tau2.t 
      parm.mat[i+1, "tau2.s"] <- 
        1/rgamma(1, nlocs/2 + a1, 
                 b1 + as.numeric(t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                                   (Dw.s - parm.mat[i, "rho.s"]*W.s)%*%
                                   matrix(theta.mat[i+1, ], ncol = 1))/2)
      
      
      # 5. update rho.t (M-H)
      inter = as.numeric( t(matrix(theta.mat[i+1, ], ncol = 1))%*%
                            W.s%*%matrix(theta.mat[i+1, ], ncol = 1) )
      ll.rho = ll.rho.s.1 + rho1.prior.val/(2*parm.mat[i+1, "tau2.s"]^2)*inter
      parm.mat[i+1, "rho.s"] <- sample(x = rho1.prior.val, size = 1, 
                                       prob = exp(ll.rho - max(ll.rho)))
      
    }
    
    
    # 9. update xi (over dispersion) 
    if(!is.null(offset.vec)){
      eta.vec = offset.vec + X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }else{
      eta.vec = X.design.fix%*%matrix(beta.mat[i+1, ], ncol = 1) + 
        as.numeric(X.design.spatial%*%theta.mat[i+1, ], ncol = 1)
      eta.vec <- as.numeric(eta.vec)
    }
    q.nb = 1/(1+exp(eta.vec)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(outcome.vec, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(outcome.vec, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, 
                            outcome.vec + parm.mat[i+1, "xi"], eta.vec)
    
    ## tuning parameters ##
    if(i <= 4000 & i%%100 == 0){
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
    }
    
  } # end MCMC
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  parm.post.comb <- cbind(beta.mat[-c(1:burn_in), ],
                          theta.mat[-c(1:burn_in), ])
  X.design.comb <- cbind(X.design.fix, 
                         X.design.spatial)
  X.design.comb <- as.matrix(X.design.comb)
  pred.val.mat <- get_predvals(parm_post = parm.post.comb,
                               X_all = X.design.comb,
                               ntotal = ntotal)
  
  compq <- function(pred.val.vec, offset.vec){
    if(!is.null(offset.vec)){
      q = 1/(1+exp(offset.vec + pred.val.vec))
    }else{
      q = 1/(1+exp(pred.val.vec))
    }
    return(q)
  }
  q.list <- lapply(1:(ncol(pred.val.mat)), 
                   function(y) compq(pred.val.vec = pred.val.mat[,y],
                                     offset.vec = offset.vec))
  q.vec <- do.call(cbind, q.list)
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(outcome.vec, q.vec), 1, 
                    compll, xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  X.design.fix.null <- X.design.fix
  X.design.fix.null[,covariate.AD] <- 0
  X.design.null <- cbind(X.design.fix.null, 
                         X.design.spatial)
  X.design.null <- as.matrix(X.design.null)
  pred.val.null.mat <- get_predvals(parm_post = parm.post.comb,
                                    X_all = X.design.null,
                                    ntotal = ntotal)
  
  if(!is.null(offset.vec)){
    pred.counts.mat <- get_count(offset = offset.vec,
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = offset.vec,
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }else{
    pred.counts.mat <- get_count(offset = rep(0, ntotal),
                                 pred_mat = pred.val.mat,
                                 overdisp_post = 
                                   as.numeric(parm.mat[-c(1:burn_in), "xi"]))
    pred.counts.null.mat <- get_count(offset = rep(0, ntotal),
                                      pred_mat = pred.val.null.mat,
                                      overdisp_post = 
                                        as.numeric(parm.mat[-c(1:burn_in), "xi"]))
  }
  AD.mat <- pred.counts.mat - pred.counts.null.mat
  
  if(!is.null(ID.spacetime.df)){
    AD.mat <- data.frame(ID.spacetime.df, AD.mat)
    pred.counts.mat <- data.frame(ID.spacetime.df, pred.counts.mat)
    pred.counts.null.mat <- data.frame(ID.spacetime.df, pred.counts.null.mat)
  }
  
  if(countsonly){
    re <- list(beta = beta.mat,
               parm = parm.mat, 
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }else{
    re <- list(beta = beta.mat, 
               parm = parm.mat, 
               rand.int.s = theta.mat,
               pred.counts = pred.counts.mat,
               pred.null.counts = pred.counts.null.mat,
               AD = AD.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), 1]), 
               WAIC = -2*(lppd - pWAIC), 
               lppd = lppd, pWAIC = pWAIC)
  }
  
  return(re)
  
}



# Functions to summarize simulation results
get_sumres <- function(x){
  re.tab = data.frame(est = colMeans(x),
                      sd = apply(x, 2, sd),
                      lci = apply(x, 2, quantile, 0.025),
                      uci = apply(x, 2, quantile, 0.975))
  return(re.tab)
}

get.counts.total <- function(x, dat.sim){
  tab.list <- list()
  jcount = 1
  names.vec = c("y.obs", "y0.true", "y.CV19.true",
                "y.obs", "y0.true", "y.CV19.true")
  for(j in c("pred.counts", "pred.null.counts", "AD",
             "pred.counts.randdraw",
             "pred.null.counts.randdraw",
             "AD.randdraw")){
    re.curr.post <- colSums(x[[j]])
    re.tab <- data.frame(true = sum(dat.sim[,names.vec[jcount]]),
                         est = mean(re.curr.post),
                         sd = sd(re.curr.post),
                         lci = quantile(re.curr.post, 0.025),
                         uci = quantile(re.curr.post, 0.975))
    rownames(re.tab) <- NULL
    tab.list[[jcount]] <- re.tab
    jcount = jcount + 1
  }
  names(tab.list) <- c("pred.counts", "pred.null.counts", "AD",
                       "pred.counts.randdraw",
                       "pred.null.counts.randdraw",
                       "AD.randdraw")
  return(tab.list)
}

get.counts.timeseries <- function(x, dat.sim){
  dat.true <- dat.sim %>% group_by(week.index) %>%
    summarise(y.obs = sum(y.obs),
              y0.true = sum(y0.true),
              y.CV19.true = sum(y.CV19.true))
  
  tab.list <- list()
  jcount = 1
  names.vec = c("y.obs", "y0.true", "y.CV19.true",
                "y.obs", "y0.true", "y.CV19.true")
  for(j in c("pred.counts", "pred.null.counts", "AD",
             "pred.counts.randdraw",
             "pred.null.counts.randdraw",
             "AD.randdraw")){
    inter <- sum.counts(counts.post = cbind.data.frame(dat.sim[,c("week.index", 
                                                                  "states")], 
                                                       x[[j]]),
                        ID.aggre = "week.index",
                        timeperiod = unique(dat.sim$week.index),
                        locs = unique(dat.sim$states),
                        ID.time = "week.index", ID.loc = "states")
    inter <- inter %>% left_join(dat.true[,c("week.index", names.vec[jcount])],
                                 by = c("week.index" = "week.index"))
    colnames(inter)[colnames(inter) == names.vec[jcount]] <- "true"
    tab.list[[jcount]] <- inter
    jcount = jcount + 1
  }
  names(tab.list) <- c("pred.counts", "pred.null.counts", "AD",
                       "pred.counts.randdraw",
                       "pred.null.counts.randdraw",
                       "AD.randdraw")
  return(tab.list)
  
}

get.counts.map <- function(x, dat.sim){
  dat.true <- dat.sim %>% group_by(states) %>%
    summarise(y.obs = sum(y.obs),
              y0.true = sum(y0.true),
              y.CV19.true = sum(y.CV19.true))
  
  tab.list <- list()
  jcount = 1
  names.vec = c("y.obs", "y0.true", "y.CV19.true",
                "y.obs", "y0.true", "y.CV19.true")
  for(j in c("pred.counts", "pred.null.counts", "AD",
             "pred.counts.randdraw",
             "pred.null.counts.randdraw",
             "AD.randdraw")){
    inter <- sum.counts(counts.post = cbind.data.frame(dat.sim[,c("week.index", 
                                                                  "states")], 
                                                       x[[j]]),
                        ID.aggre = "states",
                        timeperiod = unique(dat.sim$week.index),
                        locs = unique(dat.sim$states),
                        ID.time = "week.index", ID.loc = "states")
    inter <- inter %>% left_join(dat.true[,c("states", names.vec[jcount])],
                                 by = c("states" = "states"))
    colnames(inter)[colnames(inter) == names.vec[jcount]] <- "true"
    tab.list[[jcount]] <- inter
    jcount = jcount + 1
  }
  names(tab.list) <- c("pred.counts", "pred.null.counts", "AD",
                       "pred.counts.randdraw",
                       "pred.null.counts.randdraw",
                       "AD.randdraw")
  return(tab.list)
}


get.counts.st <- function(x, dat.sim){
  tab.list <- list()
  jcount = 1
  names.vec = c("y.obs", "y0.true", "y.CV19.true",
                "y.obs", "y0.true", "y.CV19.true")
  for(j in c("pred.counts", "pred.null.counts", "AD", 
             "pred.counts.randdraw",
             "pred.null.counts.randdraw",
             "AD.randdraw")){
    inter <- data.frame(week.index = dat.sim[,"week.index"],
                        states = dat.sim[,"states"],
                        true = dat.sim[,names.vec[jcount]],
                        est = rowMeans(x[[j]]),
                        sd = apply(x[[j]], 1, sd),
                        lci = apply(x[[j]], 1, quantile, 0.025),
                        uci = apply(x[[j]], 1, quantile, 0.975))
    tab.list[[jcount]] <- inter
    jcount = jcount + 1
  }
  names(tab.list) <- c("pred.counts", "pred.null.counts", "AD",
                       "pred.counts.randdraw",
                       "pred.null.counts.randdraw",
                       "AD.randdraw")
  return(tab.list)
}

get_sumall <- function(re1, re2, dat.curr, burn_in = 5000){
  # re1: results from fitting model with constant coefs
  # re2: results from fitting model with time-varying coefs
  
  # burn_in = 5000
  beta.hat.c <- get_sumres(x = cbind(re1$beta[-c(1:burn_in), ],
                                     rowSums(re1$beta[-c(1:burn_in), -1])))
  rownames(beta.hat.c) <- c(paste0("beta", 0:2), "beta1+2")
  
  beta.hat.t <- get_sumres(x = cbind(re2$beta[-c(1:burn_in), ],
                                     rowSums(re2$beta[-c(1:burn_in), -1])))
  rownames(beta.hat.t) <- c(paste0("beta", 0:2), "beta1+2")
  
  
  parm.hat.c <- get_sumres(x = re1$parm[-c(1:burn_in), ])
  parm.hat.t <- get_sumres(x = re2$parm[-c(1:burn_in), ])
  
  
  beta.rand.mat <- lapply(re2$beta.tvarying, 
                          function(y) get_sumres(y[-c(1:burn_in), ]))
  names(beta.rand.mat) <- c("beta1", "beta2")
  beta.rand.mat[["beta1+2"]] <- get_sumres(Reduce("+", 
                                                  re2$beta.tvarying)[-c(1:burn_in), ])
  parm.beta.t <- do.call(rbind.data.frame, 
                         lapply(re2$parm.beta.tvarying,
                                function(y) get_sumres(y[-c(1:burn_in), ])))
  
  
  rand.list <- list(time = list(cons = get_sumres(re1$rand.int.t[-c(1:burn_in), ]),
                                tvary = get_sumres(re2$rand.int.t[-c(1:burn_in), ])),
                    state = list(cons = get_sumres(re1$rand.int.s[-c(1:burn_in), ]),
                                 tvary = get_sumres(re2$rand.int.s[-c(1:burn_in), ])))
  
  
  counts.list <- list()
  counts.list[["total"]] <- list(cons = 
                                   get.counts.total(x = re1, 
                                                    dat.sim = dat.curr),
                                 tvary = 
                                   get.counts.total(x = re2, 
                                                    dat.sim = dat.curr))
  
  counts.list[["timeseries"]] <- list(cons = get.counts.timeseries(x = re1,
                                                                   dat.sim = dat.curr),
                                      tvary = get.counts.timeseries(x = re2,
                                                                    dat.sim = dat.curr))
  
  counts.list[["map"]] <- list(cons = get.counts.map(x = re1,
                                                     dat.sim = dat.curr),
                               tvary = get.counts.map(x = re2,
                                                      dat.sim = dat.curr))
  
  counts.list[["st"]] <- list(cons = get.counts.st(x = re1,
                                                   dat.sim = dat.curr),
                              tvary = get.counts.st(x = re2,
                                                    dat.sim = dat.curr))
  
  waic.list <- list(waic = c(re1$WAIC, re2$WAIC),
                    pwaic = c(re1$pWAIC, re2$pWAIC),
                    lppd = c(re1$lppd, re2$lppd))
  
  re.all <- list(beta = list(cons = beta.hat.c,
                             tvary = beta.hat.t),
                 parm = list(cons = parm.hat.c,
                             tvary = parm.hat.t),
                 rand.effct = rand.list,
                 beta.rand = beta.rand.mat,
                 beta.rand.parm = parm.beta.t,
                 counts.list = counts.list,
                 waic = waic.list)
  return(re.all)
  
}

### A function to provide summarized results of predicted counts ###
sum.counts <- function(counts.post, ID.time, ID.loc,
                       ID.aggre, timeperiod, locs){
  ## INPUTs:
  # counts.post: posterior samples (burn-in samples are discarded) of 
  ## predicted counts obtained from applying functions fit.NB.st/fit.NB.timeseries
  # ID.aggre: a string indicates the variable that posterior samples are aggregated by
  # timeperiod: time-period of interest
  # locs: locations of interest
  # ID.time: a string indicates the time indicator of posterior samples
  # ID.loc: a string indicates the space indicator of posterior samples
  
  ID.all <- c(ID.time, ID.loc)
  col.names.post <- colnames(counts.post)
  if(!(ID.aggre %in% col.names.post)){
    stop("aggregation index must be included in the dataframe containing posterior samples")
  }
  if(!all(ID.all %in% col.names.post)){
    stop("variables indicating space and time of posterior samples must be included in the dataframe")
  }
  if(!(ID.aggre %in% ID.all)){
    stop("aggregation index must be one of space/time indicators")
  }
  index1 = counts.post[,ID.time] %in% timeperiod
  index2 = counts.post[,ID.loc] %in% locs
  post.sub <- subset(counts.post, index1 & index2)
  
  ID.notaggre <- ID.all[!(ID.all %in% ID.aggre)]
  post.sub1 <- post.sub %>% dplyr::select(-ID.notaggre) %>% 
    group_by(across(ID.aggre)) %>% summarise_all(sum)
  
  get.sumres <- function(df){
    re = data.frame(est = rowMeans(df),
                    sd = apply(df, 1, sd),
                    lci = apply(df, 1, quantile, 0.025),
                    uci = apply(df, 1, quantile, 0.975))
    return(re)
  }
  
  post.sub2 <- post.sub1[,!colnames(post.sub1) %in% ID.aggre]
  re0 <- get.sumres(post.sub2)
  return(data.frame(post.sub1[,ID.aggre], re0))
}

get_sumparms <- function(res.all, true.vec, s1){
  re.tab.list <- NULL
  count = 1
  for(s1.1 in c("cons", "tvary")){
    beta.est.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"est"]))
    beta.sd.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"sd"]))
    beta.lci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"lci"]))
    beta.uci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"uci"]))
    
    re.tab <- data.frame(true = true.vec,
                         avg.est = rowMeans(beta.est.all),
                         avg.bias = rowMeans(beta.est.all) - true.vec,
                         emp.sd = apply(beta.est.all, 1, sd),
                         avg.sd = rowMeans(beta.sd.all),
                         cover = sapply(1:length(true.vec),
                                        function(x) mean(beta.lci.all[x,] <= true.vec[x] & 
                                                           true.vec[x] <= beta.uci.all[x,])))
    rownames(re.tab) <- rownames(res.all[[1]][[s1]][[s1.1]])
    re.tab.list[[count]] <- re.tab
    count = count + 1
  }
  names(re.tab.list) <- c("cons", "tvary")
  re.tab <- left_join(re.tab.list[[1]], re.tab.list[[2]],
                      by = c("true" = "true"))
  rownames(re.tab) <- rownames(res.all[[1]][[s1]][[s1.1]])
  return(re.tab)
}

get_sumbetarands <- function(res.all, beta.rand.true.mat){
  s1 = "beta.rand"
  re.tab.list <- list()
  re.sum.tab <- NULL
  for(s1.1 in 1:3){
    beta.est.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"est"]))
    beta.sd.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"sd"]))
    beta.lci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"lci"]))
    beta.uci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][,"uci"]))
    true.vec <- beta.rand.true.mat[,s1.1]
    
    re.tab <- data.frame(true = true.vec,
                         avg.est = rowMeans(beta.est.all),
                         avg.bias = rowMeans(beta.est.all) - true.vec,
                         emp.sd = apply(beta.est.all, 1, sd),
                         avg.sd = rowMeans(beta.sd.all),
                         cover = sapply(1:length(true.vec),
                                        function(x) mean(beta.lci.all[x,] <= true.vec[x] & 
                                                           true.vec[x] <= beta.uci.all[x,])))
    
    re.sum.tab <- rbind(re.sum.tab, c(avg.true = mean(true.vec),
                                      avg.bias = mean(re.tab$avg.bias),
                                      avg.cover = mean(re.tab$cover),
                                      mse = get.mse(re.tab$avg.est, re.tab$true),
                                      rse = get.rse(re.tab$avg.est, 
                                                    re.tab$true, mean(re.tab$true)),
                                      mae = get.mae(re.tab$avg.est, re.tab$true),
                                      rae = get.rae(re.tab$avg.est, 
                                                    re.tab$true, mean(re.tab$true)),
                                      mape = get.mape(re.tab$avg.est, re.tab$true),
                                      rrmse = get.rrmse(re.tab$avg.est, re.tab$true)))
    
    re.tab.list[[s1.1]] <- re.tab
  }
  names(re.tab.list) <- c("beta1_t", "beta2_t", "beta1+2_t")
  rownames(re.sum.tab) <- c("beta1_t", "beta2_t", "beta1+2_t")
  return(list(pointwise = re.tab.list,
              sum = re.sum.tab))
  
}

get.mse <- function(x.pred, x.obs){
  re = mean((x.pred - x.obs)^2)
  return(re)
}
get.rse <- function(x.pred, x.obs, x.obs.mean){
  re = sum((x.pred - x.obs)^2)/sum((x.obs - x.obs.mean)^2)
  return(re)
}
get.mae <- function(x.pred, x.obs){
  re = mean(abs(x.pred - x.obs))
  return(re)
}
get.rae <- function(x.pred, x.obs, x.obs.mean){
  re = sum(abs(x.pred - x.obs))/sum(abs(x.obs - x.obs.mean))
  return(re)
}
get.mape <- function(x.pred, x.obs){
  re = mean(abs(x.pred - x.obs)/x.obs)
  return(re)
}
get.rrmse <- function(x.pred, x.obs){
  re = sqrt(mean((x.pred - x.obs)^2)/sum(x.pred^2))
  return(re)
}
get.rse.nested <- function(x.pred, x.obs){
  nsce = nrow(x.pred)
  de <- do.call(c, lapply(1:nsce, function(x) (x.pred[x, ] - x.obs[x, ])^2))
  nu <- do.call(c, lapply(1:nsce, function(x) (x.obs[x, ] - mean(x.obs[x, ]))^2))
  return(sum(de)/sum(nu))
}

get.rae.nested <- function(x.pred, x.obs){
  nsce = nrow(x.pred)
  de <- do.call(c, lapply(1:nsce, function(x) abs(x.pred[x, ] - x.obs[x, ])))
  nu <- do.call(c, lapply(1:nsce, function(x) abs(x.obs[x, ] - mean(x.obs[x, ]))))
  return(sum(de)/sum(nu))
}

get_sumcounts <- function(res.all, s1.1, s1.2, s1.3){
  
  # s1.1.vec <- c("total", "timeseries", "map", "st")
  # s1.2.vec <- c("cons", "tvary")
  # s1.3.vec <- c("pred.counts", "pred.null.counts", "AD")
  # s1 = "counts.list"
  # s1.1 = "total"
  # s1.2 = "cons"
  # s1.3 = "AD"
  
  f.total <- function(res.all, s1 = "counts.list", s1.1 = "total", s1.2, s1.3){
    est.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"est"]))
    sd.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"sd"]))
    lci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"lci"]))
    uci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"uci"]))
    true.vec <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"true"]))
    
    re.tab <- data.frame(avg.true = mean(true.vec),
                         avg.est = rowMeans(est.all),
                         avg.bias = rowMeans(est.all - true.vec),
                         avg.rel.bias = mean((est.all - true.vec)/true.vec),
                         emp.sd = apply(est.all, 1, sd),
                         avg.sd = rowMeans(sd.all),
                         cover = mean(lci.all[1,] <= true.vec & 
                                        true.vec <= uci.all[1,]),
                         mse = get.mse(as.numeric(est.all), 
                                       as.numeric(true.vec)),
                         rse = get.rse(est.all, true.vec, mean(true.vec)),
                         mae = get.mae(as.numeric(est.all), 
                                       as.numeric(true.vec)),
                         rae = get.rae(est.all, true.vec, mean(true.vec)),
                         mape = get.mape(as.numeric(est.all), 
                                         as.numeric(true.vec)),
                         rrmse = get.rrmse(as.numeric(est.all), 
                                           as.numeric(true.vec)))
    
    return(re.tab)
  } # end f.total
  
  f <- function(res.all, s1 = "counts.list", s1.1, s1.2, s1.3){
    
    est.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"est"]))
    sd.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"sd"]))
    lci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"lci"]))
    uci.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"uci"]))
    true.all <- do.call(cbind, lapply(1:100, function(x) 
      res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"true"]))
    
    num.sces <- nrow(est.all)
    
    re.tab <- data.frame(avg.true = rowMeans(true.all),
                         avg.est = rowMeans(est.all),
                         avg.bias = rowMeans(est.all - true.all),
                         avg.rel.bias = rowMeans((est.all - true.all)/true.all),
                         emp.sd = apply(est.all, 1, sd),
                         avg.sd = rowMeans(sd.all),
                         cover = sapply(1:num.sces, function(x) 
                           mean(lci.all[x,] <= true.all[x,] & 
                                  true.all[x,] <= uci.all[x,])),
                         mse = sapply(1:num.sces, function(x) 
                           get.mse(est.all[x, ], true.all[x, ])),
                         rse = sapply(1:num.sces, function(x) 
                           get.rse(est.all[x, ], true.all[x, ],
                                   mean(true.all[x, ]))),
                         mae = sapply(1:num.sces, 
                                      function(x) get.mae(est.all[x, ],
                                                          true.all[x, ])),
                         rae = sapply(1:num.sces, 
                                      function(x) get.rae(est.all[x, ],
                                                          true.all[x, ],
                                                          mean(true.all[x, ]))),
                         mape = sapply(1:num.sces, 
                                       function(x) get.mape(est.all[x, ],
                                                            true.all[x, ])),
                         rrmse = sapply(1:num.sces, 
                                        function(x) get.rrmse(est.all[x, ],
                                                              true.all[x, ])))
    
    if(s1.1 == "st"){
      est.all.vec <- do.call(c, lapply(1:100, function(x) 
        res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"est"]))
      true.all.vec <- do.call(c, lapply(1:100, function(x) 
        res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"true"]))
      re.tab.sum <- data.frame(avg.true = mean(re.tab$avg.true),
                               avg.est = mean(re.tab$avg.est),
                               avg.bias = mean(re.tab$avg.bias),
                               avg.rel.bias = mean(re.tab$avg.rel.bias[is.finite(re.tab$avg.rel.bias)]),
                               per.upbias = mean(re.tab$avg.bias > 0),
                               avg.cover = mean(re.tab$cover),
                               mse = get.mse(est.all.vec, true.all.vec),
                               rse = get.rse(est.all.vec, 
                                             true.all.vec, mean(true.all.vec)),
                               mae = get.mae(est.all.vec, true.all.vec),
                               rae = get.rae(est.all.vec, 
                                             true.all.vec, mean(true.all.vec)),
                               mape = get.mape(est.all.vec, true.all.vec),
                               rrmse = get.rrmse(est.all.vec, true.all.vec))
      return(list(pointwise = re.tab, sum = re.tab.sum))
    }else{
      est.all.vec <- do.call(c, lapply(1:100, function(x) 
        res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"est"]))
      true.all.vec <- do.call(c, lapply(1:100, function(x) 
        res.all[[x]][[s1]][[s1.1]][[s1.2]][[s1.3]][,"true"]))
      re.tab.sum <- data.frame(avg.true = mean(re.tab$avg.true),
                               avg.est = mean(re.tab$avg.est),
                               avg.bias = mean(re.tab$avg.bias),
                               avg.rel.bias = mean(re.tab$avg.rel.bias[is.finite(re.tab$avg.rel.bias)]),
                               per.upbias = mean(re.tab$avg.bias > 0),
                               avg.cover = mean(re.tab$cover),
                               mse = get.mse(est.all.vec, true.all.vec),
                               rse = get.rse.nested(est.all, true.all),
                               mae = get.mae(est.all.vec, true.all.vec),
                               rae = get.rae.nested(est.all, true.all),
                               mape = get.mape(est.all.vec, true.all.vec),
                               rrmse = get.rrmse(est.all.vec, true.all.vec))
      return(list(pointwise = re.tab, sum = re.tab.sum))
    }
  } # end f
  
  if(s1.1 == "total"){
    re.final <- f.total(res.all = res.all,
                        s1.2 = s1.2, s1.3 = s1.3)
    return(re.final)
  }else if(s1.1 != "total"){
    re.final <- f(res.all = res.all,
                  s1.1 = s1.1,
                  s1.2 = s1.2, s1.3 = s1.3)
    return(re.final)
  }
  
  
  # s1 = "counts.list"
  # # s1.1 = "timeseries" 
  # s1.1 = "st"
  # s1.2 = "cons"
  # s1.3 = "AD"
  
} # end get_sumcounts 