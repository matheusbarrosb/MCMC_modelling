schaeferMCMC = function (data, inits, nIter, nBurnIn, delta, priors) {
  
  # RUNS A METROPOLIS-HASTINGS MCMC FOR A SCHAEFER SURPLUS PRODUCTION MODEL
  
  ## PARAMETERS ----------------------------------------------------------------
  # @ data: dataset with columns for years, surveys, and catch
  # @ inits: three-element numerical vector with initial values for MCMC 
  # @ nIter: number of total MCMC iterations
  # @ nBurnin: number of discarded iterations after sampling
  # @ delta: width of proposal distribution for r
  # @ priors: six-element numerical vector, first three elements are prior means,
  # last three elements are prior variances
  ## ---------------------------------------------------------------------------
  priorMean = c(priors[1:3])
  priorVar  = c(priors[4:6])
  
  NLL = function(par, data, priorMean, priorVar, printStep = F) {
    
    # COMPUTES NEGATIVE LOG-LIKELIHOOD OF SCHAEFER SURPLUS PRODUCTION MODEL
    
    r      = as.numeric(par[1])
    k      = as.numeric(par[2])
    sigma  = as.numeric(par[3])
    nYears = length(data$Year)
    B      = array(dim = nYears)
    B[1]   = k
    
    for (t in 2:(nYears-1)) { 
      
      # COMPUTE ESTIMATED BIOMASS GIVEN A SET OF PARAMETERS
      B[t] = max(1, B[t-1] + r*B[t-1]*(1-B[t-1]/k) - data$C[t-1])
      
    }
    
    # CATCHABILITY
    Q1   = exp(mean(log(data$Survey1/B), na.rm=TRUE))
    Q2   = exp(mean(log(data$Survey2/B), na.rm=TRUE))
    
    # SCALING SURVEY DATA
    PB1  = data$Survey1/Q1
    PB2  = data$Survey2/Q2
    
    # CALCULATE LOG SUM OF SQUARES
    SSQ1 = log(B/PB1)^2
    SSQ2 = log(B/PB2)^2
    
    # LIKELIHOOD COMPUTATIONS
    i    = which(is.na(SSQ1)==FALSE) # obtain which are NAs
    L1   = 1/sqrt(2*3.14*sigma^2)*exp(-(SSQ1[i]/(2*sigma^2)))
    NLL1 = -log(L1)
    
    j    = which(is.na(SSQ2)==FALSE) # obtain which are NAs
    L2   = 1/sqrt(2*3.14*sigma^2)*exp(-(SSQ2[j]/(2*sigma^2)))
    NLL2 = -log(L2)
    
    # PRIORS
    rPriorMean     = priorMean[1]
    kPriorMean     = priorMean[2]
    sigmaPriorMean = priorMean[3]
    
    rPriorVar      = priorVar[1]
    kPriorVar      = priorVar[2]
    sigmaPriorVar  = priorVar[3]
    
    rPriorL     = -((r - rPriorMean)^2/(2*rPriorVar))
    kPriorL     = -((k - kPriorMean)^2/(2*kPriorVar))
    sigmaPriorL = -((sigma - sigmaPriorMean)^2/(2*sigmaPriorVar))
    
    priorL = sum(rPriorL, kPriorL, sigmaPriorL)
    
    # COMPUTE TOTAL LIKELIHOOD
    totalNLL = -sum(NLL1, NLL2, -priorL)
    
    if (printStep==TRUE) print(paste("r",round(r,4),"k",k,"nll",round(totalNLL,2)))
    
    return(totalNLL)
    
  }
  
  # INITIALIZE PARAMETERS
  r         = matrix(0, nIter, 1)
  r[1,]     = inits[1]
  
  k         = matrix(0, nIter, 1)
  k[1,]     = inits[2]
  
  sigma     = matrix(0, nIter, 1)
  sigma[1,] = inits[3]
  
  pb        = txtProgressBar(min = 0,
                             max = nIter,
                             initial = 0,
                             style = 3)
  
  # SAMPLING -------------------------------------------------------------------
  for (i in 2:nIter) {
    
    setTxtProgressBar(pb, i)
    
    # PROPOSAL DISTRIBUTIONS
    rHat     = rnorm(1, r[i-1,], delta)
    kHat     = rnorm(1, k[i-1,], 20000)
    sigmaHat = abs(rnorm(1, sigma[i-1,], 0.05))
    
    # COMPUTING LOG ACCEPTANCE RATIO
    logR = NLL(par = c(rHat, kHat, sigmaHat),
               priorMean = priorMean, priorVar = priorVar, data = data) - 
      NLL(par = c(r[i-1,], k[i-1,], sigma[i-1,]),
          priorMean = priorMean, priorVar = priorVar, data = data) 
    
    # DECIDE WHETHER TO ACCEPT OR REJECT PROPOSED VALUES
    # ACCEPTS WITH PROBABILITY P = min(U(1), R)
    if (runif(1) <= exp(logR)) {
      
      r[i,]     = rHat
      k[i,]     = kHat
      sigma[i,] = sigmaHat
      
    } else {
      
      r[i,]     = r[i-1,]
      k[i,]     = k[i-1,]
      sigma[i,] = sigma[i-1,]
      
    }
    
    close(pb)
  }
  # END OF SAMPLING ------------------------------------------------------------
  
  rSamples     = r[(nBurnIn:nIter),]
  kSamples     = k[(nBurnIn:nIter),]
  sigmaSamples = sigma[(nBurnIn:nIter),]
  
  # PLOTTING
  par(mfrow = c(2,3))
  plot(r, type = "l")
  plot(k, type = "l")
  plot(sigma, type = "l", ylab = expression(sigma))
  
  plot(density(rSamples), main = "r")
  plot(density(kSamples), main = "k")
  plot(density(sigmaSamples), main = expression(sigma))
  
  # MAKING OUTPUT --------------------------------------------------------------
  MCMC_samples = data.frame(rSamples, kSamples, sigmaSamples)
  MCMC_summary = data.frame(c(mean(rSamples),
                              mean(kSamples),
                              mean(sigmaSamples)),
                            c(sd(rSamples, na.rm = TRUE),
                              sd(kSamples),
                              sd(sigmaSamples)))
  names(MCMC_summary) = c("Mean", "SD")
  
  output = list(MCMC_samples, MCMC_summary); names(output) = c("MCMC_samples", "MCMC_summary")
  
  return(output)
}

D=read.csv(file="Widow data for homework 2.csv",header=TRUE,sep=",")
names(D) = c("Year", "C", "Survey1", "Survey2")

priors = c(
  0.23, 144000, 0.9, # MEANS
  0.025, 20000, 0.01 # SDS
) 

a = schaeferMCMC(data = D, inits = c(0.19, 144000, 0.8),
                 priors = priors,
                 nIter = 10000000, nBurnIn = 100000, delta = 0.0075)

################################################################################

B = matrix(NA, length(D$Year), 1)
B[1,] = mean(a[[1]]$kSamples)

for(t in 2:length(D$Year)) {
  B[t] = B[t-1] + mean(a[[1]]$rSamples)*B[t-1]*(1-B[t-1]/mean(a[[1]]$kSamples)) - D$C[t-1]
}

rSD = a[[2]]$SD[1]
kSD = a[[2]]$SD[2]

Bupper = matrix(NA, length(D$Year), 1)
Bupper[1,] = mean(a[[1]]$kSamples) + kSD
for(t in 2:length(D$Year)) {
  Bupper[t] = Bupper[t-1] + (mean(a[[1]]$rSamples)+rSD)*Bupper[t-1]*(1-Bupper[t-1]/(mean(a[[1]]$kSamples)+kSD)) - D$C[t-1]
}

Blower = matrix(NA, length(D$Year), 1)
Blower[1,] = mean(a[[1]]$kSamples) - kSD
for(t in 2:length(D$Year)) {
  Blower[t] = Blower[t-1] + (mean(a[[1]]$rSamples)-rSD/4)*Blower[t-1]*(1-Blower[t-1]/(mean(a[[1]]$kSamples)-kSD/2)) - D$C[t-1]
}

Q1 = exp(mean(log(D$Survey1/biomass$B), na.rm = T))
Q2 = exp(mean(log(D$Survey2/biomass$B), na.rm = T))

biomass = data.frame(B, Bupper, Blower, D$Survey1/Q1, D$Survey2/Q2)
names(biomass) = c("B", "Bup", "Blow", "S1", "S2")

plot(biomass$S1, xlim = c(20,120), ylim = c(5000, 200000))
points(biomass$S2, col = "red")
lines(biomass$B)
lines(biomass$Bup, lty = 2)
lines(biomass$Blow, lty = 2)

k = a[[2]]$Mean[2]
r = a[[2]]$Mean[1]
plot(biomass$B/(k/2), type = "l", main = "B/Bmsy")

FM = D$C/biomass$B 
FM_Fmsy = FM/FM[68]
plot(FM_Fmsy, type = "l", main = "F/Fmsy")