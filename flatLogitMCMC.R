MHflatLogit = function(nIter, nBurnIn, y, X, delta, plot = TRUE) {
  
  ## PERFORMS A METROPOLIS-HASTINGS ALGORITHM FOR A PROBIT MODEL
  
  ## PARAMETERS
  
  # nIter -> number of iterations
  # y     -> binary response data
  # X     -> design matrix of covariates
  # delta -> witdh of the proposal distribution
  # plot  -> LOGICAL, plot traceplots and density distributions?
  
  ## RETURNS
  
  # array with samples from coefficient posterior distributions
  
  if (nIter <= nBurnIn) stop("Number of warmp-up iterations needs to be smaller than total iterations")
  
  logitLL = function(beta, y, X) {
    
    ## COMPUTES THE LOG-LIKELIHOOD FOR A LOGIT MODEL
    
    ## PARAMETERS
    
    # beta -> model coefficients, n x p matrix
    # y    -> binary response data
    # X    -> design matrix of covariates
    
    ## RETURNS
    
    # log-likelihood of parameter values given the data
    
    if (is.matrix(beta)==F) beta = as.matrix(t(beta))
    n   = dim(beta)[1]
    pll = rep(0,n)
    for (i in 1:n) {
      lf1    = plogis(X%*%beta[i,], log.p = T)
      lf2    = plogis(-X%*%beta[i,], log.p = T)
      pll[i] = sum(y*lf1+(1-y)*lf2) 
    }
    pll
  }
  
  p        = dim(X)[2] # get number of coefficients
  freqMod  = summary(glm(y~1+X,family = binomial(link = "logit")))
  beta     = matrix(0, nIter, p) # empty matrix to store betas, could also be NAs?
  beta[1,] = as.vector(freqMod$coeff[,1]) # set initial values
  variance = as.matrix(freqMod$cov.unscaled) # get initial variance
  
  pb       = txtProgressBar(min = 0,
                            max = nIter,
                            initial = 0,
                            style = 3) # terminal progress bar
  
  for(i in 2:nIter) {
    
    setTxtProgressBar(pb,i)
    
    # THE PROPOSAL DISTRIBUTION
    betaHat = rmnorm(1, beta[i-1,], delta*variance)
    
    # COMPUTE LOG ACCEPTANCE RATIO
    logR    = logitLL(betaHat, y, X) - logitLL(beta[i-1,], y, X)
    
    # ACCEPT WITH PROBABILITY MIN(1, R)
    if (runif(1) <= exp(logR)) beta[i,] = betaHat
    else beta[i,] = beta[i-1,]
    
    close(pb)
  }
  
  # EXCLUDE WARMP-UP ITERATIONS
  beta = beta[-(1:nBurnIn),]
  
  if (plot == TRUE) {
    
    par(mfrow = c(2,3))
    # PLOTTING TRACEPLOTS
    plot(beta[,1], type = "l", main = expression(alpha))
    plot(beta[,2], type = "l", main = paste(expression(beta),"1"))
    plot(beta[,3], type = "l", main = paste(expression(beta),"2"))
    
    # PLOTTING PROBABILITY DENSITIES
    plot(density(beta[,1]), main = expression(alpha))
    plot(density(beta[,2]), main = paste(expression(beta),"1"))
    plot(density(beta[,3]), main = paste(expression(beta),"2"))
    
  }
  
  # POSTERIOR SUMMARIES
  
  require(dplyr)
  
  summary = data.frame(
    c("alpha","beta1","beta2"),
    c(median(beta[,1]), median(beta[,2]), median(beta[,3])),
    c(sd(beta[,1]), sd(beta[,2]), sd(beta[,3]))
  )
  
  names(summary) = c("Parameter", "Mean", "SD")
  colnames(beta) = c("alpha", "beta1", "beta2")
  
  output = list(summary, beta)
  names(output) = c("posteriorSummary", "posteriorSamples")
  
  output
  
}

# RUNNING MCMC
a = MHflatLogit(nIter = 50000, nBurnIn = 1000,
                y = y, X = X, delta = 0.5, plot = TRUE)
