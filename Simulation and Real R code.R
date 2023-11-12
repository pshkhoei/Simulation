# -------------------------------------------------------------------------
       #1 Run Bayesian Approach with Weibull distribution.
# -------------------------------------------------------------------------  

{
    rm(list = ls())
  # Install packages:survival & R2openBUGS.
    library(survival)
    library(R2OpenBUGS)
  # Set working directory and modelfile.
    getwd()
    bugswd = paste0(getwd(),"/bugswd"); bugswd
    modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Generate Data
    set.seed(12345)
    n = 200  # n=100; 200; 300
    x = rep(0:1, c(0.50*n, 0.50*n))  # Weibull scale parameter related to x.
    table(x)
    shape = 2  # Shape: 0.5; 1; 2
    b = c(-3, 0.3) # set b1 and b2 with table 2 in paper.
    lambda = exp(b[1] + b[2]*x)   # Link the parameter to covariate x .
    summary(lambda)
    scale = lambda^(-1/shape)   # Since weibull formulla in winbugs is different to R, we need to convert
                              # formula to get similar results.
    summary(scale)   # Mean scale parameter is near to 4.
  #Generate Observed time
    y = rweibull(n,shape, scale )    
    summary(y)
    range(y)
  # Generate censored time
    delta1 = rep(1,n)   # to make censored data
    cen = rexp(n,0.06)             # Censored time
    delta = as.numeric(y < cen)
    cenper = 1 - mean(delta); cenper   # Get percent of censoring
  # Merge observed and censored time.
    z = pmin(y,cen)  # to select observed time or censored time. Every one that is lesser than other.
  # make variable "t" as observed time and variable "c" as censored time to use in BUGS.
    t <- ifelse(delta == 1, z, NA)
    c <- ifelse(delta == 1, 0, z)
  # Run model in BUGS.
      modeltext = "model {
      for(i in 1:n){
        t[i] ~ dweib(shape,lambda[i])C(c[i], )
      	log(lambda[i]) <- b[1]+b[2]*x[i]
        cim[i] <- step(c[i]-1.0E-5)*pow(log(2)/lambda[i]+pow(c[i],shape), 1/shape)
        }
      	# priors
      	shape ~ dgamma(0.01,0.01)  # Non-informative prior
      	for(j in 1:2) {b[j]~dnorm(0,0.01)}		
      }
      "
      # write BUGS output into file.
      cat(modeltext, file = modelfile) #file.show(modelfile)
      modeldata = list(n = n, x = x, t = t, c = c)
      modelinit = list(list(b = rep(0,length(b)), shape = shape))
      param = c("shape","b","cim")
      # bugs ----------------------------------------
      bugsOut <- bugs(
        working.directory = bugswd,
        model.file = modelfile,
        data = modeldata,
        inits = modelinit,
        #inits = NULL,
        parameters.to.save = param,
        n.chains = 1,
        n.iter = 11000,
        n.burnin = 1000,
        n.thin = 10
        #, debug = TRUE
        #, codaPkg = TRUE
      )
  # output ----------------------------------------
    bugsOut$DIC
    # Which records is censored:
    ic = which(delta==0); ic; length(ic)
    # Dimension of output:
    dim(bugsOut$sims.array)
  # Describe censored simulations.
    bugsOut$summary[c(1:3,3+ic),c(1,2)] 
    # Describe parameter simulations:
    parsim = bugsOut$sims.array[,1,1:3]   #parameter simulation
    parsim[1:5,]  # Only five rows of 10.000 simulation for parameters.
  # print median of simulations for every censor that replaced.
    bugsOut$median$cim[ic]  
  #
}

# Convergence: Geweke
  library(coda)
  geweke.diag(parsim, frac=0.10, frac2 = 0.50)
  matparsin <- mcmc(as.matrix(parsim))
  geweke.plot(matparsin)
  acfplot(matparsin)
  autocorr.diag(matparsin)
  traceplot(matparsin)

#-------------------------------------------------------------------------------
#          Figures 7 in paper, figures 1 and 2 in supplementary material.

{
  # Kaplan-Meier Curve:
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("t~Weibull(2,4), c~Exp(0.06), p=0.20, n=200") ) 
  
  # Curve with Median of Simulated Times
    # output ----------------------------------------
    # imputation      h=hat
    bh = bugsOut$mean$b; bh
    shapeh = bugsOut$mean$shape; shapeh
    lambdah = exp(bh[1] + bh[2]*x); lambdah #every person has specific lambda because it has specific X.
    scaleh = lambdah^(-1/shapeh); scaleh
    # Compute median of Simulations.
    zmed = qweibull(.5*pweibull(cen,shapeh,scaleh, lower.tail = FALSE),shapeh, scaleh, lower.tail = FALSE)
  
    zimp = rep(NA,n)
    zimp[ic] = zmed[ic]
    zimp[-ic] = z[-ic]  # zimp = failure times+imputed censored times
  
    curve2 = survfit(Surv(zimp,delta1) ~ x); curve2     # Bayesian Imputation
    lines(curve2, mark.time = TRUE, col = "Blue", lty = 1)
  
  #Curve without Censored Times 
    tOC = z[delta==1]   #time omitting censored
    deltaOC = rep(1, length(tOC))
    curve3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); curve3       # Omitting_Censored
    lines(curve3, mark.time = TRUE, col = "Red", lty = 1)
    
    legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
}

#-------------------------------------------------------------------------------
  #           Figure 4 in paper. 

{
  # Kaplan-Meier Curve:  
  km1 = survfit(Surv(z,delta) ~ x); km1
  plot(km1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
       main = paste("t~Weibull(2,4), c~Exp(0.06), p=0.20, n=200"))  #KM_Estimation
  
  # Curve for 10,000 Times Imputation. 
  timp=t
  impsim = bugsOut$sims.array[,1,3+ic]
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    
  }
  # Curve for Imputations Mean
  timp[ic] <- colMeans(impsim)
  kmmean = survfit(Surv(timp,delta1) ~ x)
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
  
}



# -------------------------------------------------------------------------
  #2 Run Bayesian Approach with Birnbaum-Saunders (BS) distribution.
# -------------------------------------------------------------------------  

{
    rm(list = ls())
    # Set working directory and modelfile.
    getwd()
    bugswd = paste0(getwd(),"/bugswd"); bugswd
    modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
    
    # generate Data
    set.seed(12345)
    n = 200 # n=100; 200; 300
    x = rep(0:1, c(0.50*n, 0.50*n))  # BS scale parameter related to x.
    table(x)
    shape = 2    # Shape: 0.5; 1; 2                   
    b = c(1.37, 0.15)  # set b1 and b2 with table 4 in paper.
    lambda = exp(b[1] + b[2]*x)
    summary(lambda)
    scale = lambda             
    # Define rbn to generate numbers from BS distribution.
    rbn <- function(n, shape, scale){    # shape = a, scale = b
      x <- rnorm(n, 0, shape/2)
      t <- scale * (1 + 2 * x^2 + 2 * x * sqrt(1 + x^2))
      return(t)
    }
    #Generate Observed time
    y <- rbn(n, shape, lambda)
    # Generate censored time
    delta1 = rep(1,n)
    cen = rexp(n,0.021)             # Censored time
    delta = as.numeric(y < cen)
    cenper = 1 - mean(delta); cenper   # % censoring
    # Merge observed and censored time.
    z = pmin(y,cen)
    # make variable "t" as observed time and variable "c" as censored time to use in BUGS.
    t <- ifelse(delta == 1, z, NA)
    c <- ifelse(delta == 1, 0, z)
    # Run model in BUGS.
    modeltext = "model {
  for(i in 1:n){
  	t[i] ~ dbs(shape, lambda[i])C(c[i], )
  	log(lambda[i]) <- b[1]+b[2]*x[i]
    cim[i] <- step(c[i]-1.0E-5)*lambda[i]
    
    }
  	# priors
  	shape ~ dgamma(0.01,0.01)
  	for(j in 1:2) {b[j]~dnorm(0,0.01)}		
  }
  "
    # write BUGS output into file.
    cat(modeltext, file = modelfile) #file.show(modelfile)
    modeldata = list(n = n, x = x, t = t, c = c)
    modelinit = list(list(b = rep(0,length(b)), shape = shape))
    param = c("shape","b","cim")
    # bugs ----------------------------------------
    bugsOut <- bugs(
      working.directory = bugswd,
      model.file = modelfile,
      data = modeldata,
      inits = modelinit,
      #inits = NULL,
      parameters.to.save = param,
      n.chains = 1,
      n.iter = 11000,
      n.burnin = 1000,
      n.thin = 10
      #, debug = TRUE
      #, codaPkg = TRUE
    )
    # output ----------------------------------------
    bugsOut$DIC
    # Which records is censored:
    ic = which(delta==0); ic
    # Dimension of output:
    dim(bugsOut$sims.array)
    # Describe censored simulations.
    bugsOut$summary[c(1:3,3+ic),c(1,2)]  
    # Describe parameter simulations:
    parsim = bugsOut$sims.array[,1,1:3]    #parameter simulation
    parsim[1:5,]  # Only five rows of 10.000 simulation for parameters.
  }  
  
# Convergence: Geweke
  library(coda)
  geweke.diag(parsim, frac=0.10, frac2 = 0.50)
  matparsin <- mcmc(as.matrix(parsim))
  geweke.plot(matparsin)
  acfplot(matparsin)
  autocorr.diag(matparsin)
  traceplot(matparsin)


  #------------------------------------------------------
  # Figures 5 in paper, figures 3 and 4 in supplementary material.

  {
    # Kaplan-Meier Curve:
    curve1 = survfit(Surv(z,delta) ~ x); curve1
    plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
         main = paste("t~BS(2,4), c~Exp(0.02), p=0.20, n=200") )  #KM_Estimation
    # Curve with Median of Simulated Times
    # output ----------------------------------------
    # imputation      h=hat
    bh = bugsOut$mean$b; bh
    shapeh = bugsOut$mean$shape; shapeh
    lambdah = exp(bh[1] + bh[2]*x); lambdah #every person has specific lambda because it has specific X.
    scaleh = lambdah; scaleh
    
    #install.packages("extraDistr")
    library(extraDistr)
    # Compute median of Simulations.
    zmed = qfatigue(.5*pfatigue(cen,shapeh,scaleh, mu = 0, lower.tail = FALSE),shapeh, scaleh,mu = 0, lower.tail = FALSE)
    zmed
    # Make a variable include median of simulations.
    zimp = rep(NA,n)
    zimp[ic] = zmed[ic]
    zimp[-ic] = z[-ic]  # zimp = failure times+imputed censored times
    
    
    #
    curve2 = survfit(Surv(zimp,delta1) ~ x); curve2
    lines(curve2, mark.time = TRUE, col = "Blue", lty = 1)
    
    #Curve without Censored Times 
    tOC = z[delta==1]
    deltaOC = rep(1, length(tOC))
    curve3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); curve3
    lines(curve3, mark.time = TRUE, col = "Red", lty = 1)
    
    legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"),
           lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
    
  }
  
  
#-------------------------------------------------------------------------------
#          Figure 6 in paper.

{
    # Kaplan-Meier Curve:  
    curve1 = survfit(Surv(z,delta) ~ x); curve1
    plot(curve1, mark.time = TRUE,lty = 1, lwd = 2, col = "black", 
         main = paste("t~BS(2,4), c~Exp(0.02), p=0.20, n=200"))  #KM_Estimation
    
    # Curve for 10,000 Times Imputation
    timp=t
    impsim = bugsOut$sims.array[,1,3+ic]  
    for (i in 1:nrow(impsim)) {
      timp[ic] <- impsim[i,]
      kmi = survfit(Surv(timp,delta1) ~ x)
      lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
      
    }
    # Curve for Imputations Mean
    timp[ic] <- colMeans(impsim)
    kmmean = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
    
    legend("topright", 
           c("Kaplan-Meier Curve", "Curve for 10,000 times Imputation", "Curve for Imputations Mean"),
           lty = 1, col = c("Black", "gray","blue"), cex = .7)
    
    
  }


  
  
# -------------------------------------------------------------------------
   #3 Run Bayesian Approach on Breast Cancer Data distributed as Weibull.
# -------------------------------------------------------------------------  
{
  rm(list = ls())
  # Install packages:survival & R2openBUGS.
  library(survival)
  library(R2OpenBUGS)
  # Set working directory and modelfile.
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Import and define variables in Data.
  breast <- read.table("Data_Paper1.txt", header = TRUE)
  t <- breast$t 
  c <- breast$c
  x <- breast$AgeC
  length(t[t == "NA"])/length(t)   # Percent of Censoring, 88 Censor, 40% 
  length(c[c == "0"])/length(c)   # Percent of Observed
  n = length(t); n
  z = breast$z   # Composed from Observed and Censored data
  delta = breast$delta   # delta=0 means Censoring
  ic = which(delta == "0")   # indicator censor
  age <- breast$AgeC
  # Run model in BUGS.
  modeltext = "model {
    for(i in 1:n){
    t[i] ~ dweib(shape,lambda)C(c[i], )
    cim[i]<-step(c[i]-1.0E-5)*pow(log(2)/lambda+pow(c[i],shape),1/shape)  
    }
  	# priors
  	shape ~ dgamma(0.01,0.01)
  	lambda ~ dgamma(0.01, 0.01)
	}
  "
  # write BUGS output into file.
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list(n = n, t = t, c = c)
  modelinit = list(list(shape = 1, lambda = 1 ))
  param = c("shape","lambda", "cim")
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 10
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  
# output ----------------------------------------
  bugsOut$DIC
  # Dimension of output:
  dim(bugsOut$sims.array)   #composed: alpha, lambda, 88 simulation,deviance = 91 columns.
  # Describe censored simulations.
  bugsOut$sims.array[1:5,1,3:90]         # Head
  bugsOut$sims.array[9996:10000,1,3:90]  # Tail
  bugsOut$summary[1:2, c(1:2)]    # mean & sd parameters: alpha & lambda
  # Describe parameter simulations:
  parsim = bugsOut$sims.array[,1,1:2]   #parameter simulation
  impsim = bugsOut$sims.array[,1,3:90]  # imputation simulation
  timp = t
}

# Convergence: Geweke
  library(coda)
  geweke.diag(parsim, frac=0.10, frac2 = 0.50)
  matparsin <- mcmc(as.matrix(parsim))
  geweke.plot(matparsin)
  acfplot(matparsin)
  autocorr.diag(matparsin)
  traceplot(matparsin)
  
#------------------------------------------------------
                # Fiqure 7 in Paper.

{
  # Kaplan-Meier Curve:
  curve1 = survfit(Surv(z,delta) ~ age); curve1
  plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("Posterior Estimate: Shape=1.24,Scale=0.001,DIC=1698"))  #KM_Estimation
  # Curve with Median of Simulated Times
  # output ----------------------------------------
  # imputation      h=hat
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = bugsOut$mean$lambda; lambdah
  cen=c
  # Compute median of Simulations.
  library(miscTools)
  zmed = colMedians(impsim)
  # 
  ic = which(delta==0); ic      #index censor to count number of censored case.
  length(ic)
  zimp = rep(NA,n)
  zimp[ic] = zmed[ic]
  zimp[-ic] = z[-ic]  # zimp = failure times+imputed censored times
  delta1 = rep(1,n)  # after impute, all of times are observed then we made delta1.
  #
  km2 = survfit(Surv(zimp,delta1) ~ x); km2     # Bayesian Imputation
  lines(km2, mark.time = TRUE, col = "Blue", lty = 1)
  
  # Curve without Censored Times
  tOC = z[delta==1]  # number of observed times
  deltaOC = rep(1, length(tOC))
  length(deltaOC)
  km3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); km3       # Omitting_Censored
  lines(km3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
}


#-------------------------------------------------------------------------------
#           Fiqure 7 in Paper.
{

  # Kaplan-Meier Curve
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1, lwd =2, col = "black",
       main = paste("t~Weibull, p=0.40, n=220"))  #KM_Estimation
  
  # Curve with Median of Simulated Times
  # simulation 
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    #Sys.sleep(.5)
  }
  # Curve for Imputations Mean
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
  }  


# -----------------------------------------------------------------------------------
    #4 Run Bayesian Approach on Breast Cancer Data distributed as Bernbaum-Saunders.
# -----------------------------------------------------------------------------------
{
  rm(list = ls())
  # Install packages:survival & R2openBUGS.
  library(survival)
  library(R2OpenBUGS)
  # Set working directory and modelfile.
  getwd()
  bugswd = paste0(getwd(),"/bugswd"); bugswd
  modelfile = paste0(bugswd,"/modelfile.txt"); modelfile
  # Import and define variables in Data.
  breast <- read.table("Data_Paper1.txt", header = TRUE)
  x <- breast$AgeC  
  t <- breast$t  #time based on month
  c <- breast$c
  length(t[t == "NA"])/length(t)   # Percent of Censoring, 88 Censor, 40% 
  length(c[c == "0"])/length(c)   # Percent of Observed
  n = length(x); n
  z = breast$z   # Composed from Observed and Censored data
  delta = breast$delta   # delta=0 means Censoring
  ic = which(delta == 0)   # indicator censor
  length(ic)
  age <- breast$AgeC
  # Run model in BUGS.
  modeltext = "model {
  for(i in 1:n){
  t[i] ~ dbs(shape,lambda)C(c[i], )
  cim[i] <- step(c[i]-1.0E-5)*lambda    #tmed
  }
  # priors
  shape ~ dgamma(0.01,0.01)
  lambda ~ dgamma(0.01, 0.01)
  }
  "
  # write BUGS output into file.
  cat(modeltext, file = modelfile) #file.show(modelfile)
  modeldata = list(n = n, t = t, c = c)
  modelinit = list(list(shape = 4, lambda = 4))
  param = c("shape","lambda", "cim")
  # bugs ----------------------------------------
  bugsOut <- bugs(
    working.directory = bugswd,
    model.file = modelfile,
    data = modeldata,
    inits = modelinit,
    #inits = NULL,
    parameters.to.save = param,
    n.chains = 1,
    n.iter = 11000,
    n.burnin = 1000,
    n.thin = 10
    #, debug = TRUE
    #, codaPkg = TRUE
  )
  
  # output ----------------------------------------
  bugsOut$DIC
  # Dimension of output:
  dim(bugsOut$sims.array)   #composed: alpha, lambda, 88 simulation,deviance = 91 columns.
  # Describe censored simulations.
  bugsOut$sims.array[1:5,1,3:90]  # report 1 till 5 from 100 times censored times simulations.
  bugsOut$summary[1:2, c(1:2)]    # mean & sd parameters: alpha & lambda
  # Describe parameter simulations:
  parsim = bugsOut$sims.array[,1,1:2]   #parameter simulation: 10000*2
  impsim = bugsOut$sims.array[,1,3:90]  # imputation simulation: 10000*88
  timp = t
  
}

# Convergence: Geweke
  library(coda)
  geweke.diag(parsim, frac=0.10, frac2 = 0.50)
  matparsin <- mcmc(as.matrix(parsim))
  geweke.plot(matparsin)
  acfplot(matparsin)
  autocorr.diag(matparsin)
  traceplot(matparsin)
  
  
#------------------------------------------------------
#                Fiqure 8

{ 
  # Kaplan-Meier Curve:
  curve1 = survfit(Surv(z,delta) ~ age); curve1
  plot(curve1, mark.time = TRUE,lty = 1,conf.int = FALSE,  col = "black",
       main = paste("Posterior Estimate: Shape=1.22, Scale=145.21, DIC=1510"))  #KM_Estimation
  
  # Curve with Median of Simulated Times
  # output ----------------------------------------
  # imputation      h=hat
  shapeh = bugsOut$mean$shape; shapeh
  lambdah = bugsOut$mean$lambda; lambdah
  scaleh = lambdah; scaleh
  cen=c
  # Compute median of Simulations.
  #install.packages("extraDistr")
  library(extraDistr)
  # How calculate median times in Birnbaum-Saunders distribution:
  zmed = qfatigue(.5*pfatigue(cen,shapeh,scaleh, mu = 0, lower.tail = FALSE),shapeh, scaleh,mu = 0, lower.tail = FALSE)
  #
  ic = which(delta==0); ic      #index censor to count number of censored case.
  zimp <- rep(NA, n)
  zimp[ic] <- zmed[ic]
  zimp[-ic] <- z[-ic]  # zimp = failure times+imputed censored times
  delta1 = rep(1,n)  # after impute, all of times are observed then we made delta1.
  #
  curve2 = survfit(Surv(zimp,delta1) ~ x); curve2     # Bayesian Imputation
  lines(curve2, mark.time = TRUE, col = "Blue", lty = 1)
  
  # Curve without Censored Times
  tOC = z[delta==1]  # number of observed times
  deltaOC = rep(1, length(tOC))
  length(deltaOC)
  curve3 = survfit(Surv(tOC, deltaOC) ~ x[delta==1]); curve3       # Omitting_Censored
  lines(curve3, mark.time = TRUE, col = "Red", lty = 1)
  
  legend("topright", c("Kaplan-Meier Curve", "Curve with Median of Simulated Times", "Curve without Censored Times"), lty= 1, col = c("black", "Blue", "Red"), cex = 0.7)
  
}  


#-------------------------------------------------------------------------------
                     # Fiqure 8

{  
  # Kaplan-Meier Curve    
  curve1 = survfit(Surv(z,delta) ~ x); curve1
  plot(curve1, mark.time = TRUE,lty = 1, lwd=2, col = "black",
       main = paste("t~Birnbaum-Saunders, p=0.40, n=220"))  #KM_Estimation
  # Curve with Median of Simulated Times  
  # simulation 
  for (i in 1:nrow(impsim)) {
    timp[ic] <- impsim[i,]
    kmi = survfit(Surv(timp,delta1) ~ x)
    lines(kmi, mark.time = TRUE, col = "gray", lty = 1)    # n time Imputation
    #Sys.sleep(.5)
  }
  # Curve for Imputations Mean
  timp[ic] <- colMeans(impsim)
  kmmean = survfit(Surv(timp,delta1) ~ x)
  lines(kmi, mark.time = TRUE, col = "blue", lty = 2, lwd = 2)  # Mean of n times Imputation
  
  legend("topright", 
         c("Kaplan-Meier Curve", "Curve for 10,000 Times Imputation", "Curve for Imputations Mean"),
         lty = 1, col = c("Black", "gray","blue"), cex = .7)
}  

