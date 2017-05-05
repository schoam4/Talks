
## Compare mediation power analysis with MC conf int and Bootstrap conf int

## Generate data from population covariance matrix
## Use example from Hayes (2013) study on econmic stress by Pollack et al. (2012)

## X = estress (economic stress)
## M = affect (depressed affect)
## Y = withdraw (withdrawl intentions)

##Make one big function for MClapply
N = 20
simfunc <- function(N){

popCor <- matrix(c(1, 0.340, 0.064, 0.340, 1, 0.417, 0.064, 0.417, 1), nrow = 3, byrow = TRUE)
popSD <- diag(c(1.42361, 0.72372, 1.24687))
popCov <- popSD%*%popCor%*%popSD


rep = 1000
bootRep = 1000
MCrep = 20000

## Generate data

library(MASS)

#function to generate 1 data set

metaRep <- 500
MCtime <- rep(NA, metaRep)
MCpow <- rep(NA, metaRep)
MCpowSE <- rep(NA, metaRep)
boottime <- rep(NA, metaRep)
bootpow <- rep(NA, metaRep)
bootpowSE <- rep(NA, metaRep)

for(i in 1:metaRep){

genMed <- function(seed = 1234, pCov = popCov, pMeans = c(0,0,0), Ns = N){
   require(MASS)
   mvrnorm(Ns, mu = pMeans, Sigma = pCov)
  
}

##Generate R data sets using apply

dats <- lapply(sample(1:50000, rep), genMed)


##Analyze results with MC conf int

powMC <- function(dat, mcmcReps = MCrep, conf = 95){
  
  #require(MASS)
  
  # Run regressions
  m1 <- lm(dat[,2] ~ dat[,1])
  m2 <- lm(dat[,3] ~ dat[,2] + dat[,1])
  
  # Output parameter estimates and standard errors
#   pest <- c(, coef(m2)[2])
#   covmat <- diag(c((diag(vcov(m1)))[2],
#                    (diag(vcov(m2)))[2]))
  
  # Simulate draws of a, b from multivariate normal distribution
  a <- rnorm(mcmcReps, coef(m1)[2], sqrt(vcov(m1)[2,2]))
  b <- rnorm(mcmcReps, coef(m2)[2], sqrt(vcov(m2)[2,2]))
  ab <- a*b
  
  # Calculate confidence intervals
  low <- (1 - (conf / 100)) / 2
  upp <- ((1 - conf / 100) / 2) + (conf / 100)
  LL <- quantile(ab, low)
  UL <- quantile(ab, upp)
  
  # Is rep significant?
  LL*UL >= 0
  
}

MCtime[i] <- system.time(pow <- lapply(dats, powMC))
MCpow[i] <- sum(unlist(pow)) / rep
MCpowSE[i] <- sqrt((MCpow*(1-MCpow))/rep)

## Lets try with lavaan


powboot <- function(dat, nBoot = bootRep, conf = 95){
  
  require(lavaan)
  
  # Run mediation
  dat <- data.frame(dat)
  
  mod <- 'X2 ~ a*X1
          X3 ~ X1 + b*X2
          ab := a*b
         '
  fit <- sem(mod, dat, boot= nBoot)
  
  ci <-  parameterestimates(fit)[7,]
  
  # Is rep significant?
  ci[9]*ci[10] >= 0
  
}

boottime[i] <- system.time(pow <- lapply(dats, powboot))
bootpow[i] <- sum(unlist(pow)) / rep
bootpowSE[i] <- sqrt((bootpow*(1-bootpow))/rep)

}

return(cbind(MCtime,MCpow, MCpowSE, boottime, bootpow, bootpowSE))
}

#now for a vector of sample sizes
N <- seq(10, 300, 10)

#mclapply hack for windows

source('http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R')

res <- mclapply(N, simfunc)

save(res, "results.RData")

# 
# powboot <- function(dat, nBoot = bootrep, conf = 95){
#   
#   require(boot)
#   
#   # Run mediation
#   
#   bootMed <- function(dat, nBoot) {
#     require(boot)
#     
#     #make sure variable names are characters
#     # dv <- 'advance'
#     # iv <- 'complaints'
#     # med <- 'learning'
#     # numBoot=50
#     
#     #Create function to pass to boot
#     medReg <- function(dat, i){
#       
#       d <- dat[i,]
#       
#       m1 <- lm(d[,2] ~ d[,1])
#       m2 <- lm(d[,3] ~ d[,2] + d[,1])
#       
#       ab <- (coef(m1)[2])*(coef(m2)[2])
#       
#       return(as.numeric(ab))
#     }
#     #Bootstrap
#     medBoot <- boot(dat, statistic=medReg, R=nBoot, sim = "ordinary", stype = "i")
#     return(medBoot)
#   }
#   
#   ci <-  boot.ci(test, type = "bca")$bca
#   
#   # Is rep significant?
#   ci[4]*ci[5] >= 0
#   
# }
# 
# Boottime <- system.time(pow <- lapply(dats, powboot))
# Bootpow <- sum(unlist(pow)) / rep
# BootpowSE <- sqrt((Bootpow*(1-Bootpow))/rep)
# Boottime
# Bootpow
# BootpowSE


powMC <- function(dat, mcmcReps = MCrep, conf = 95){
  
  require(MASS)
  require(lavaan)
  
  # Run mediation
  dat <- data.frame(dat)
  
  mod <- 'X2 ~ a*X1
  X3 ~ X1 + b*X2
  ab := a*b
  '
  fit <- sem(mod, dat, boot= nBoot)
  
  est <-  parameterestimates(fit)
  # Output parameter estimates and standard errors
  pest <- c(est[1,5], est[3,5])
  
  covmat <- diag(c(vcov(fit)[1,1], vcov(fit)[3,3]))
  
 
  ## Simulate draws of a, b from multivariate normal distribution
  mcmc <- mvrnorm(mcmcReps, pest, covmat, empirical = FALSE)
  ab <- mcmc[, 1] * mcmc[, 2]
  
  # Calculate confidence intervals
  low <- (1 - (conf / 100)) / 2
  upp <- ((1 - conf / 100) / 2) + (conf / 100)
  LL <- quantile(ab, low)
  UL <- quantile(ab, upp)
  
  # Is rep significant?
  LL*UL >= 0
  
}

MCtime <- system.time(pow <- lapply(dats, powMC))
MCpow <- sum(unlist(pow)) / rep
MCpowSE <- sqrt((MCpow*(1-MCpow))/rep)


### Check out results

res2 <- lapply(res, colMeans)

res3 <- do.call(rbind.data.frame, res2)
names(res3) <- names(res2[[1]])
summary(res3)
res3
res3 <- cbind(res3, N)
res3
cor(res3)

write.csv(res3, 'simres.csv')

