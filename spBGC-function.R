library(TruncatedNormal)
library(tmvtnorm)

###-----------------------------------------------------------------###
###-----------------------------------------------------------------###
###             spBGC: spatial Bayesian Gaussian Copula             ###
###-----------------------------------------------------------------###
###-----------------------------------------------------------------###


# Input -------------------------------------------------------------------
# data: (n,p)-data matrix 
# SP: (n,2)-matrix of location information
# CF: correlation function; default is "exponential correlation function"
# Phi.set: set of spatial range parameters 
# m: number of neighbors in NNGP; default is 5
# iter: length of MCMC 
# warmup: length of burn-in 
# R0: (p,p)-positive definite matrix of inverse-Wishart prior; default is identity matrix
# r0: a positive scalar degrees of freedom parameter; default is p+2

# Output (list object) ----------------------------------------------------
# post.R: MCMC samples of R (correlation matrix of spatial Gaussian Copula)
# post.phi: MCMC samples of phi (spatial range parameter)
# acceptance prob.: acceptance probability of phi in MH step


# Main function -----------------------------------------------------------
spBGC <- function(data, SP, CF=function(dist.mat, phi){ exp(-dist.mat/phi) },
                  Phi.set=NULL, m=5,
                  iter=2500, warmup=500, 
                  R0=diag(dim(data)[2]), r0=dim(data)[2]+2) {
  # Check prior -------------------------------------------------------------
  ok_R0 <- all(eigen(R0)$val>0) & dim(R0)[1]==dim(data)[2] & dim(R0)[2]==dim(data)[2]
  ok_r0 <- (r0>=0)
  if(!ok_R0) { stop("Error: R0 must be a positive definite p x p matrix \n") }
  if(!ok_r0) { stop("Error: r0 must be positive \n") }
  
  # Preparation -------------------------------------------------------------
  vnames <- colnames(data)
  data <- as.matrix(data)
  colnames(data) <- vnames
  n <- dim(data)[1]
  p <- dim(data)[2]
  Rank_data <- apply(data, 2, function(x){ match(x, sort(unique(x))) })
  Rlevels <- apply(Rank_data, 2, max, na.rm=TRUE)
  Ranks <- apply(data, 2, rank, ties.method="max",na.last="keep")
  N <- apply(!is.na(Ranks), 2, sum)
  U <- t( t(Ranks)/(N+1) )
  
  ## initial values ----------
  Z <- qnorm(U) # latent variable 
  Zfill <- matrix(rnorm(n*p),n,p)
  Z[is.na(data)] <- Zfill[is.na(data)] # imputation of latent variable for missing data 
  V <- cov(Z) # Gaussian covariance matrix 
  L <- length(Phi.set)
  if(is.null(Phi.set)){
    L <- 10
    Phi.range <- c(0, median(dist(SP)))
    Phi.set <- seq(Phi.range[1], Phi.range[2], length=L+1)[-1]
  }
  Phi.index <- ifelse(L==1,1,round(L/2))  # spatial range parameter \phi
  
  ## missing index ----------
  na_index <- which(is.na(data), arr.ind=TRUE)
  na_id <- na_index[,1] + n*{na_index[,2]-1}
  
  ## NNGP parameter ----------
  dd <- as.matrix( dist(SP) )
  NN <- matrix(NA, n, m) # neighbor children location 
  for(i in 2:n){
    if(i<=m){  NN[i, 1:(i-1)] <- sort(order(dd[i,1:(i-1)]))  }
    if(i>m){   NN[i,] <- sort(order(dd[i,1:(i-1)])[1:m])  }
  }
  In <- function(x, i){ i %in% x }
  Parent <- lapply(1:n, function(i){ (1:n)[apply(NN, 1, In, i=i)] }) # parent location
  Dependence <- lapply(1:n, function(i){ if(i==1) c(1, Parent[[i]]) else c(NN[i,][!is.na(NN[i,])], i, Parent[[i]]) })  # dependent location
  BB.set <- list()    # coefficient of latent variables
  FF.set <- list()    # conditional variance of latent variables
  for(l in 1:L){
    BB <- list()
    FF <- c()
    for(i in 1:n){
      mat <- CF(as.matrix( dist(SP[Dependence[[i]],])), Phi.set[l])
      C1 <- solve(mat[which(i!=Dependence[[i]]), which(i!=Dependence[[i]])])
      C2 <- as.matrix(mat[which(i==Dependence[[i]]), which(i!=Dependence[[i]])])
      FF[i] <- drop( mat[which(i==Dependence[[i]]), which(i==Dependence[[i]])] - emulator::quad.form(C1,C2) )
      BB[[i]] <- drop( crossprod(C2,C1) )
    }
    BB.set[[l]] <- BB
    FF.set[[l]] <- FF
  }
  BB.set_Phi <- list()    
  FF.set_Phi <- list()    
  for(l in 1:L){
    BB <- matrix(NA, n, m)
    FF <- c()
    FF[1] <- 1
    for(i in 2:n){
      if(i<=m){
        mat <- CF(as.matrix( dist(SP[c(NN[i,1:(i-1)], i),])), Phi.set[l])
        C1 <- solve(mat[1:(i-1), 1:(i-1)])
        C2 <- mat[i, 1:(i-1)]
        FF[i] <- mat[i, i] - emulator::quad.form(C1,C2) 
        BB[i,1:(i-1)] <- as.vector(C2%*%C1) 
      }
      if(i>m){
        mat <- CF(as.matrix(dist(SP[c(NN[i,], i),])), Phi.set[l])
        C1 <- solve(mat[1:m, 1:m])
        C2 <- mat[m+1, 1:m]
        FF[i] <- mat[m+1, m+1] - emulator::quad.form(C1,C2) 
        BB[i,] <- drop( crossprod(C2,C1) ) 
      }
    }
    BB.set_Phi[[l]] <- BB
    FF.set_Phi[[l]] <- FF
  }
  
  
  # MCMC iterations ---------------------------------------------------------
  c <- 1
  post_R <- array(dim=c(p,p,(iter-warmup)))
  post_phi <- numeric((iter-warmup))
  prob_acp <- numeric(iter)
  pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                   total = iter,
                                   complete = "=",   # Completion bar character
                                   incomplete = "-", # Incomplete bar character
                                   current = ">",    # Current bar character
                                   clear = FALSE,    # If TRUE, clears the bar when finish
                                   width = 100)      # Width of the progress bar
  for (r in 1:iter) {
    pb$tick() # progress bar
    
    ## update latent variable ----------
    for (i in 1:n) {
      bound <- (sapply(1:p, function(j){
        lb <- suppressWarnings(max( Z[ which((Rank_data[i,j]-1)==Rank_data[,j]),j],na.rm=TRUE))
        ub <- suppressWarnings(min( Z[ which((Rank_data[i,j]+1)==Rank_data[,j]),j],na.rm=TRUE))
        matrix(c(lb,ub),nrow=2)
      }))
      Di <- Dependence[[i]][i!=Dependence[[i]]]
      mu_Zi <- colSums(BB.set[[Phi.index]][[i]]*Z[Di,])
      Sigma_Zi <- FF.set[[Phi.index]][i]*V
      if(i %in% na_id) {
        js <- na_index[na_index[,1]==i,2]
        Z[i,-js] <- TruncatedNormal::rtmvnorm(1, mu=mu_Zi, sigma=Sigma_Zi, lb=bound[1,], ub=bound[2,])[-js]
        for (j in js) {
          sig_na <- drop( Sigma_Zi[j,j] - emulator::quad.form.inv(Sigma_Zi[-j,-j], Sigma_Zi[j,-j]) )
          mu_na <- mu_Zi[j] + drop( emulator::quad.3form.inv(Sigma_Zi[-j,-j],Sigma_Zi[j,-j],(Z[i,-j]-mu_Zi[-j])) )
          Z[i,j] <- rnorm(1,mu_na,sqrt(sig_na))
        }
      }else Z[i,] <- TruncatedNormal::rtmvnorm(1, mu=mu_Zi, sigma=Sigma_Zi, lb=bound[1,], ub=bound[2,])
    }
    
    ## update correlation matrix ----------
    ZHZ <- matrix(0,p,p)
    for (i in 1:n) {
      if(i==1) {
        ZHZ <- ZHZ + tcrossprod(Z[i,])
      }else if(i==2) {
        bar_Z_i <- BB.set_Phi[[Phi.index]][i,1:(i-1)]*Z[NN[i,1:(i-1)],]
        ZHZ <- ZHZ + FF.set_Phi[[Phi.index]][i]^{-1} * tcrossprod(Z[i,]-bar_Z_i)
      }else if(i<=m) {
        bar_Z_i <- colSums(BB.set_Phi[[Phi.index]][i,1:(i-1)]*Z[NN[i,1:(i-1)],])
        ZHZ <- ZHZ + FF.set_Phi[[Phi.index]][i]^{-1} * tcrossprod(Z[i,]-bar_Z_i)
      }else {
        bar_Z_i <- colSums(BB.set_Phi[[Phi.index]][i,]*Z[NN[i,],])
        ZHZ <- ZHZ + FF.set_Phi[[Phi.index]][i]^{-1} * tcrossprod(Z[i,]-bar_Z_i)
      }
    }
    V <- round(MCMCpack::riwish( r0+n, r0*R0 + ZHZ ), 7)
    R <- round( V/tcrossprod(sqrt(diag(V)),sqrt(diag(V))), 7 )
    
    ## update spatial range parameter ----------
    new.index <- Phi.index + sample(c(1, -1), 1) 
    if(new.index<1){ new.index <- 1 }
    if(new.index>L){ new.index <- L }
    new.BB_Phi <- BB.set_Phi[[new.index]]
    new.FF_Phi <- FF.set_Phi[[new.index]]
    BB_Phi <- BB.set_Phi[[Phi.index]]
    FF_Phi <- FF.set_Phi[[Phi.index]]
    L_old <- 0
    L_new <- 0
    for (i in 1:n) {
      if(i==1) {
        L_old <- L_old + drop(TruncatedNormal::.dmvnorm_arma(t(Z[i,]), rep(0,p), FF_Phi[i]*R, TRUE))
        L_new <- L_new + drop(TruncatedNormal::.dmvnorm_arma(t(Z[i,]), rep(0,p), new.FF_Phi[i]*R, TRUE))
      }else if(i<=m) {
        L_old <- L_old + drop(
          TruncatedNormal::.dmvnorm_arma(t(Z[i,]),
                                         drop( tcrossprod(t(BB_Phi[i,1:(i-1)])%x%diag(p), t(as.vector(t(Z[NN[i,1:(i-1)],]))) ) ),
                                         FF_Phi[i]*R, TRUE)
        )
        L_new <- L_new + drop(
          TruncatedNormal::.dmvnorm_arma(t(Z[i,]),
                                         drop( tcrossprod(t(new.BB_Phi[i,1:(i-1)])%x%diag(p), t(as.vector(t(Z[NN[i,1:(i-1)],]))) ) ),
                                         new.FF_Phi[i]*R, TRUE)
        )
      }else{
        L_old <- L_old + drop(
          TruncatedNormal::.dmvnorm_arma(t(Z[i,]),
                                         drop( tcrossprod(t(BB_Phi[i,])%x%diag(p), t(as.vector(t(Z[NN[i,],]))) ) ),
                                         FF_Phi[i]*R, TRUE)
        )
        L_new <- L_new + drop(
          TruncatedNormal::.dmvnorm_arma(t(Z[i,]),
                                         drop( tcrossprod(t(new.BB_Phi[i,])%x%diag(p), t(as.vector(t(Z[NN[i,],]))) ) ),
                                         new.FF_Phi[i]*R, TRUE)
        )
      }
    }
    pp <- min(1, exp(L_new-L_old))
    if(runif(1)<pp){ Phi.index <- new.index }
    
    ## store the posterior samplings ----------
    prob_acp[r] <- pp
    if(r>warmup){
      post_R[,,c] <- R
      post_phi[c] <- Phi.set[Phi.index]
      c=c+1
    }
  }
  
  # Summary posterior sampling ----------------------------------------------
  post = list(post_R, post_phi, prob_acp)
  names(post) = c("post.R", "post.phi", "acceptance prob.")
  return(post)
}




