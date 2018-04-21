## estimating

library(MCMCpack)
library(mvtnorm)
#options(digits = 10)
set.seed(1)
# Importing data

# iStart <- 1996
# iEnd   <- 2017.917
# iFreq  <- 12
# mY     <- read.csv("data-raw/brmacrom.csv", dec=",")
# mY     <- ts(mY, start = iStart, frequency = iFreq)



iP  <- 12  #defasagens

iStart <- 1996
iEnd   <- 2017.75
iFreq  <- 4
mY     <- read.csv("data-raw/brmacroq.csv", dec=",")
mY     <- ts(mY, start = iStart, frequency = iFreq)

# library(seasonal)
# 
# aux <- seas(x = mY[,1])
# mY[,1] <- (aux$data[,"trend"])
# 
# aux <- seas(x = mY[,4])
# mY[,4] <- (aux$data[,"trend"])


mY[,1] <- log(mY[,1]/1000)
mY[,4] <- log(mY[,4])

mdY    <- diff(mY)

plot(mY)
plot(mdY)

mTP <- read.csv("data-raw/brrecess.csv", dec=",")


mdY <- scale(mdY, scale = F)
cN  <- ncol(mdY)
cT  <- nrow(mdY)


# OLS

mx <- NULL
for (i in 1:iP) {
  mx <- ts.union(mx, lag(mdY, -i))
}
ms_p       <- window(mdY, iStart + 1/iFreq, iStart + iP / iFreq)
vs_p       <- c(t(ms_p))
my         <- window(mdY, iStart + (1/iFreq) + (iP / iFreq), iEnd)
my         <- matrix(my, cT - iP, cN)
mx         <- window(mx, iStart + (1/iFreq) + (iP / iFreq), iEnd)
mx         <- matrix(mx, cT - iP, iP * cN)
lOLS       <- ApplyOLS(my, mx)
mPhi_OLS   <- lOLS$mB
mSigma_OLS <- lOLS$mS

# Empirical Bayes

lhyper  <- optim(c(1, cN + 2), lnML, method = "L-BFGS-B", lower = c(.0001, cN - 1))
dLambda <- lhyper$par[1]
dN_0    <- lhyper$par[2]
#lhyper  <- optimize(lnML, c(0, 100))
#dLambda <- lhyper$minimum
#dN_0    <- cN + 2

# Prior

lPrior  <- SpecifyPrior(mY, iP, dLambda, dN_0)
mM_0    <- lPrior$mM_0
mD_0    <- lPrior$mD_0
dN_0    <- lPrior$dN_0
mS_0    <- lPrior$mS_0
mD_0Inv <- solve(mD_0)

# Posterior

lPosterior <- SpecifyPosterior(mM_0, mD_0, dN_0, mS_0, my, mx)
mM_1       <- lPosterior$mM_1
mD_1       <- lPosterior$mD_1
dN_1       <- lPosterior$dN_1
mS_1       <- lPosterior$mS_1
mD_1Inv    <- solve(mD_1)

# Log marginal likelihood

dlnML <- ComputeLogMarginalLikelihood(mD_0, dN_0, mS_0, mD_1, dN_1, mS_1)

# Bayes factor of AR vs VAR (S-D density ratio)

dbf <- ComputeSDRatio(mD_0Inv, dN_0, mS_0, mM_1, mD_1Inv, dN_1, mS_1)

# Initial value

msigma <- riwish(dN_1, mS_1)
mphi   <- DrawPhi(msigma, mM_1, mD_1Inv)
#mgamma <- ComputeInitStateVariance(mphi, msigma)

# M-H algorithm

cBurn    <- 0
cR       <- 10000
mphi_s   <- mcmc(matrix(, cR, iP * cN ^ 2))
msigma_s <- mcmc(matrix(, cR, cN * (cN + 1) / 2))
for (r in (1 - cBurn):cR) {
  
  # Draw parameters
  
  msigma_star <- riwish(dN_1, mS_1)
  mphi_star   <- DrawPhi(msigma_star, mM_1, mD_1Inv)
  #	mgamma_star <- ComputeInitStateVariance(mphi_star, msigma_star)
  #	dalpha      <- ComputePmove(mgamma, mgamma_star, vs_p)
  #	if (runif(1) <= dalpha) {
  mphi   <- mphi_star
  msigma <- msigma_star
  #		mgamma <- mgamma_star
  #	}
  
  # Save draws
  
  if (r >= 1) {
    mphi_s[r,]   <- c(t(mphi))		
    msigma_s[r,] <- vech(msigma)
  }
  print(r)
}

# B-N Decomposition

amphi   <- array(dim = c(cN, iP * cN, cR))
amsigma <- array(dim = c(cN, cN, cR))
for (r in 1:cR) {
  amphi[,, r]   <- t(matrix(mphi_s[r,], iP * cN, cN))
  amsigma[,, r] <- xpnd(msigma_s[r,], cN)
}

ms <- NULL
for (i in 0:(iP - 1)) {
  ms <- ts.union(ms, lag(mdY, -i))
}
ms      <- window(ms, iStart + iP / iFreq, iEnd)
ms      <- matrix(ms, cT - iP + 1, iP * cN)
amgap   <- array(dim = c(cN, cT - iP + 1, cR))
amdgap  <- array(dim = c(cN, cT - iP, cR))
amd2gap <- array(dim = c(cN, cT - iP - 1, cR))
amcorr  <- array(dim = c(cN, cN, cR))
mC      <- cbind(diag(cN), matrix(0, cN, (iP - 1) * cN))
for (r in 1:cR) {
  ma            <- GetCompanionMatrix(amphi[,, r])
  mw            <- -mC %*% solve(diag(iP * cN) - ma) %*% ma
  amgap[,, r]   <- mw %*% t(ms)
  amdgap[,, r]  <- t(diff(t(amgap[,, r])))
  amd2gap[,, r] <- t(diff(t(amgap[,, r]), 2))
  amcorr[,, r]  <- cor(t(amgap[,, r]))
}
mgap_med   <- t(apply(amgap, 1:2, median))
mgap_lower <- t(apply(amgap, 1:2, quantile, prob = .025))
mgap_upper <- t(apply(amgap, 1:2, quantile, prob = .975))
mgap_med   <- ts(mgap_med, start = iStart + iP / iFreq, frequency = iFreq)
mgap_lower <- ts(mgap_lower, start = iStart + iP / iFreq, frequency = iFreq)
mgap_upper <- ts(mgap_upper, start = iStart + iP / iFreq, frequency = iFreq)
mY         <- window(mY, iStart + iP / iFreq, iEnd)
mnr_med    <- mY - mgap_med
mnr_lower  <- mY - mgap_upper
mnr_upper  <- mY - mgap_lower
amgapDI    <- (amgap > 0)
mgapDI     <- t(apply(amgapDI, 1:2, mean))
mgapDI     <- ts(mgapDI, start = iStart + iP / iFreq, frequency = iFreq)
amets      <- (amd2gap[, 1:(cT - iP - 3),] > 0) * (amdgap[, 2:(cT - iP - 2),] > 0) * (amgap[, 3:(cT - iP - 1),] > 0) * (amdgap[, 3:(cT - iP - 1),] < 0) * (amd2gap[, 3:(cT - iP - 1),] < 0)
amcts      <- (amd2gap[, 1:(cT - iP - 3),] < 0) * (amdgap[, 2:(cT - iP - 2),] < 0) * (amgap[, 3:(cT - iP - 1),] < 0) * (amdgap[, 3:(cT - iP - 1),] > 0) * (amd2gap[, 3:(cT - iP - 1),] > 0)
mets       <- t(apply(amets, 1:2, mean))
mcts       <- t(apply(amcts, 1:2, mean))
mets       <- ts(mets, start = iStart + iP / iFreq + .5, frequency = iFreq)
mcts       <- ts(mcts, start = iStart + iP / iFreq + .5, frequency = iFreq)

# Stats

mY[, 4]         <- exp(mY[, 4])
mnr_med[, 4]    <- exp(mnr_med[, 4])
mnr_lower[, 4]  <- exp(mnr_lower[, 4])
mnr_upper[, 4]  <- exp(mnr_upper[, 4])
mgap_med[, 4]   <- mY[, 4] - mnr_med[, 4];
mgap_lower[, 4] <- mY[, 4] - mnr_upper[, 4];
mgap_upper[, 4] <- mY[, 4] - mnr_lower[, 4];


save(mTP, mY, mnr_med, mgap_med, mgap_lower, mgap_upper,
     amcorr, mgapDI, mets, mcts,
     file = "resultados.RData")
