# SIMULATING DATA WITH KNOWN COEFFICIENTS

y = rep(0,1000)  
y[1] = 0.25
for(i in 2:1000) {
  y[i] = 0.15 + y[i-1]*0.4 + rnorm(1, 0, 0.05) 
}

graphics.off()
plot(y, type = "l")

data = y

AR_LL = function(theta, y = data, nIter = 1) {
  # RETURNS THE LOG-LIKELIHOOD OF AN AR(1) MODEL
  # ARGUMENTS
  # @ theta: two-dimensional vector with alpha and beta values
  # @ y: sample
  # @ nIter: number of iterations for MCMC
  N     = nIter
  LL    = rep(0, N)
  alpha = theta[1]
  beta  = theta[2]
  for (i in 1:N) {
    for (k in 2:length(y)) {
      logpdf = 0.5*dnorm(y[k] - alpha - beta*y[k-1], log = TRUE)
      LL[i] = -0.5*dnorm(y[1] - alpha, log = TRUE)^2 - sum(logpdf^2)
    }
  }
  -LL
}


# EVALUATING LOG-LIKELIHOOD FOR ALPHA AND BETA
X = rep(NA, 1000)
Y = rep(NA, 1000)
alpha = seq(0.01, 1, along.with = X)
beta  = seq(0, 1.5, along.with = X)

for (i in 1:length(X)) {
  X[i] = AR_LL(theta = c(alpha[i], 0.8), y = y, nIter = 1)
  Y[i] = AR_LL(theta = c(0.1, beta[i]), y = y, nIter = 1)
}

par(mfrow = c(1,2))
plot(X ~ alpha, type = "l")
plot(Y ~ beta, type = "l")



logL <- function(par) {
  
  alpha = par[1]
  beta  = par[2]
  
  -sum(dnorm(y[-1],alpha+beta*y[1:length(y)-1],log=TRUE))
  
}

# TESTING OPTIMIZATION
res = optim(par = list(alpha = 0.3, beta = 0.1), fn = logL, hessian = T)

AR_LL = function(par, y) {
  
  alpha = as.numeric(par[1])
  beta  = as.numeric(par[2])
  y = as.numeric(y)
  
  -sum(dnorm(y[-1],alpha+beta*y[1:length(y)-1],log=TRUE))
  
  
}

res = optim(par = list(alpha = 0.3, beta = 0.1), fn = AR_LL, hessian = T)


AR_LL(par = list(alpha = 0.3, beta = 0.1), y = y)

-sum(dnorm(y[-1],0.1+0.7*y[1:length(y)-1],log=TRUE))

# MCMC -------------------------------------------------------------------

MH_AR1 = function(nIter, nBurnIn, y, delta, inits, plot = TRUE) {
  
  if (nIter <= nBurnIn) stop("Number of warmp-up iterations needs to be smaller than total iterations")
  
  AR_LL = function(par, y) {
    
    # RETURNS LOG-LIKELIHOOD OF AR(1) MODEL
    
    alpha = as.numeric(par[1])
    beta  = as.numeric(par[2])
    y     = as.numeric(y)
    
    sum(dnorm(y[-1],alpha+beta*y[1:length(y)-1],log=TRUE))
    
  }
  
  # INITIALIZE PARAMETERS
  alpha     = matrix(0, nIter, 1)
  alpha[1,] = inits[1]
  beta      = matrix(0, nIter, 1)
  beta[1,]  = inits[2]
  
  pb        = txtProgressBar(min = 0,
                             max = nIter,
                             initial = 0,
                             style = 3)
  
  # SAMPLING -------------------------------------------------------------------
  for (i in 2:nIter) {
    
    setTxtProgressBar(pb, i)
    
    # PROPOSAL DISTRIBUTIONS
    betaHat  = rnorm(1, beta[i-1,],  delta)
    alphaHat = rnorm(1, alpha[i-1,], delta)
    
    # COMPUTING LOG ACCEPTANCE RATIO
    logR     = AR_LL(par = list(alpha = alphaHat, beta = betaHat), y = y) - 
      AR_LL(par = list(alpha = alpha[i-1,], beta = beta[i-1,]), y = y)
    
    # DECIDE WHETHER TO ACCEPT OR REJECT
    # ACCEPTS WITH PROBABILITY P = min(U(1), R)
    if (runif(1) < exp(logR)-0.95) alpha[i,] = alphaHat
    else alpha[i,] = alpha[i-1,]
    
    if (runif(1) < exp(logR)-0.95) beta[i,] = betaHat
    else beta[i,] = beta[i-1,]
    
    close(pb)
    
  }
  # END OF SAMPLING STAGE ------------------------------------------------------
  
  # EXCLUDE BURN-IN ITERATIONS
  alphaSamples = alpha[(nBurnIn:nIter),]
  betaSamples  = beta[(nBurnIn:nIter),]
  
  
  # PLOTTING RESULTS
  if (plot == TRUE) {
    
    par(mfrow = c(2,2))
    
    plot(alpha, type = "l", main = expression(alpha))
    plot(beta, type = "l", main = expression(beta))
    
    plot(density(alphaSamples), main = expression(alpha))
    plot(density(betaSamples), main = expression(beta))
    
  } 
  
  # RETURNING OUTPUT
  MCMC_samples        = data.frame(as.vector(alphaSamples),
                                   as.vector(betaSamples))
  names(MCMC_samples) = c("alpha", "beta") 
  MCMC_summary        = data.frame(c(mean(MCMC_samples$alpha),
                                     mean(MCMC_samples$beta)),
                                   c(sd(MCMC_samples$alpha),
                                     sd(MCMC_samples$beta)))
  colnames(MCMC_summary) = c("Mean", "SD")
  rownames(MCMC_summary) = c("alpha", "beta")
  
  output = list(MCMC_samples, MCMC_summary)
  names(output) = c("MCMC_samples", "MCMC_summary")
  
  return(output)
}

inits = c(0, 0)

a = MH_AR1(1000000, 10000, y = y,
           delta = .01, inits = inits)
