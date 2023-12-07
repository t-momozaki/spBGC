###-----------------------------------------------------------###
###               Oneshot simulation study for                ###
###             spatial Bayesian Gaussian Copula              ###
###-----------------------------------------------------------###
rm(list=ls())


# Load packages -----------------------------------------------------------
library(mvtnorm)
library(sbgcop)
library(coda)
library(bayesplot)
source("spBGC-function.R")


# Setting -----------------------------------------------------------------
set.seed(3)
n <- 300           # sample size 
p <- 4             # dimension of mixed outcomes
phi <- 0.5         # spatial range parameter 
iter <- 3000       # length of MCMC
warmup <- 1000     # length of burn-in 


# Data generation ---------------------------------------------------------
## correlation matrix -----------------------------------------------------
R <- matrix(NA, p, p)  # 
for(j in 1:p){
  for(k in 1:p){
    R[j, k] <- (0.6)^(abs(k-j))
  }
}
comb = expand.grid(1:p,1:p)
names.par = apply( subset(comb,Var1<Var2), 1, function(x) paste0("R",paste0(x,collapse="")) )
names(names.par) = NULL
## Location information ---------------------------------------------------
SP <- matrix(runif(n*2, -2,2), n, 2)
## Data generation --------------------------------------------------------
CF <- function(dist.mat, phi) exp(-dist.mat/phi)
H <- CF(as.matrix(dist(SP)), phi=phi)
Z <- drop( mvtnorm::rmvnorm(1, mean=rep(0,p*n), sigma=H%x%R, method="chol") )
mat_Z <- t(matrix(Z, p, n))
data <- cbind(qbinom(pnorm(mat_Z[,1]), 1, 0.2),  # Binomial(0.2)
              qpois(pnorm(mat_Z[,2]), 5),        # Poisson(5)
              qpois(pnorm(mat_Z[,3]), 10),       # Poisson(10)
              mat_Z[,-(1:3)])                    # Multivariate Normal


# MCMC (spBGC; proposed) --------------------------------------------------
L <- 30
Phi.range <- c(0, quantile(dist(SP),0.50))
Phi.set <- seq(Phi.range[1], Phi.range[2], length=L+1)[-1]
post_spBGC <- spBGC(data=data, SP=SP, iter=iter, warmup=warmup, Phi.set=Phi.set)


# MCMC (BGC; Hoff, 2007) --------------------------------------------------
post_BGC <- sbgcop::sbgcop.mcmc(data, nsamp=iter, odens=1, plugin.threshold = 9999)


# Summary -----------------------------------------------------------------
pos_spBGC <- post_spBGC$post.R
post_R_spBGC <- c()
MC <- (iter-warmup)
for (r in 1:MC) {
  post_R_spBGC <- rbind(post_R_spBGC, pos_spBGC[,,r][upper.tri(pos_spBGC[,,r])])
}
post_R_spBGC <- as.data.frame(post_R_spBGC)
colnames(post_R_spBGC) <- names.par

pos_BGC <- post_BGC$C.psamp
post_R_BGC <- c()
for (r in (warmup+1):iter) {
  post_R_BGC <- rbind(post_R_BGC, pos_BGC[,,r][upper.tri(pos_BGC[,,r])])
}
post_R_BGC <- as.data.frame(post_R_BGC)
colnames(post_R_BGC) <- names.par

## Posterior mean ---------------------------------------------------------
(hR_spBGC <- colMeans(post_R_spBGC))  # proposed
(hR_BGC <- colMeans(post_R_BGC))      # Hoff (2007)
### MSE of correlation matrix ---------------------------------------------
mean( (hR_spBGC-R[upper.tri(R)])^2 )  # proposed
mean( (hR_BGC-R[upper.tri(R)])^2 )    # Hoff (2007)
## 95% credible interval --------------------------------------------------
(CI_R_spBGC <- t(apply(post_R_spBGC, 2, quantile, probs=c(0.025,0.975))))  # proposed
(CI_R_BGC <- t(apply(post_R_BGC, 2, quantile, probs=c(0.025,0.975))))      # Hoff (2007)
## Traceplot --------------------------------------------------------------
bayesplot::mcmc_trace(post_R_spBGC)
## Auto-correlation function ----------------------------------------------
bayesplot::mcmc_acf(post_R_spBGC)
## ESS of correlation matrix ----------------------------------------------
coda::effectiveSize(post_R_spBGC)



