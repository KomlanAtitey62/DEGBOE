DEGBOE_atitey <- function(n,M,H,dataset){
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% resampling function
  
  resample <- function(x, size, replace = TRUE, prob = NULL)
  {
    if(length(x)<1)
      if(!missing(size) && size>0)
        stop("Requested sample of size ", size, " from list of length 0")
    else
      x[FALSE]
    else if(length(x)==1)
    {
      if(missing(size) || size==1)
        x
      else if(size>=1 && replace==TRUE)
        rep(x, size)
      else if(size < 1)
        x[FALSE]
      else
        stop("Cannot cannot take a sample larger than the population",
             " when 'replace = FALSE'")
    }
    else
      sample(x, size, replace, prob)
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computation of the correlation matrix using the Browmian Law 
  
  N <- 1000
  n <- 100
  H = 0.58
  coef = 1.1 #0.61 #the inverse of the summation of the AR and MA orders 1/(a+b) which represents the mean
  #TimeStep =101
  nu.square = 1.78
  mu = rep(0,N)
  mu <- data.frame(mu)
  mu <- t(t(mu))
  par=c(nu.square,coef)
  
  library(MASS)
  
  cor <- c()
  for(j in 1:N){
    cor[j] = 0.5*(((abs(j+1))^(2.*H))-2*((abs(j))^(2.*H))+((abs(j-1))^(2.*H)))
  }
  cor.matrix <- toeplitz(cor)
  sigma <- par[1]*cor.matrix # variance of the mulvariate Normal distribution
  noise.dist <- rnorm(mu,sigma) # prior non stationnary states
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% arma coefficient
  dataset <- 100*as.numeric(dataset[2:length(dataset)])
  
  obs.vec <- c()
  obs.mat <- c()
  gam.fun <- c()
  arma.coef.var <- c()
  for(i in 1:n){
    obs.vec[i] <- dataset[i] # observation vector
    obs.mat[i] <- obs.vec[i]*t(obs.vec[i])
    gam.fun[i] <- (1/par[2])*obs.mat[i] # coef.^-1.*ObsMat(1,j); % R^-1(y*y')
    arma.coef <- exp(sum(diag(gam.fun[i]))) # initial value of the distribution of ARMA coef P = exp(Tr(R^-1(y*y')))
    arma.coef.var[i] <- exp(sum(diag(gam.fun[i]))) # prior distribution P = exp(Tr(R^-1(y*y')))
  }
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% smc
  
  sir.est <- c()
  variance <- matrix(0,nrow=n,ncol=N)
  inv.cov.mat <- matrix(0,nrow=n,ncol=N)
  particles = matrix(0,nrow=n,ncol=N)
  
  ## Calculate the prior density %%%%%%%%%%%%%%%
  for(t in 1:n){
    for(k in 1:N){
      particles[,k] <- arma.coef*noise.dist[k] # importance density is the transition density 
    }
    ## Calculate importance weights %%%%%%%%%%%%%%% 
    obs.vec[t] <- dataset[t] # exp(P53(1,t)./2)*vEps(t,1); % observation vector
    variance[t,] = exp(particles[t,])
    inv.cov.mat[t,] = 1/abs(variance[t,]) # inverse of the covariance matrix
    
    weights.pred = rep(1,N)/N
    weights.pred <- data.frame(weights.pred)
    weights_pred <- t(t(weights.pred))
    log.weight.0 = log(weights.pred)
    log.weight.1 = 0.0015; 
    log.weight.2 = -0.5*n*log(2*pi);
    log.weight.3= -0.5*log(abs(variance[t,]))
    log.weight.4 = -0.5*(obs.vec[t])*inv.cov.mat[t,]*obs.vec[t]
    log.weight.5 = sum(t(obs.vec[t] - particles[t,])*arma.coef.var[t]*(obs.vec[t] - particles[t,]))
    log.weight = log.weight.0 + log.weight.1 + log.weight.2 + log.weight.3 + log.weight.4 + log.weight.5
    
    ## Stabilize the importance weight %%%%%%%%%%%%%%% 
    d.max.weight <- max(log.weight)
    v.weight <- exp(log.weight - d.max.weight)
    v.weight <- v.weight/sum(v.weight)
    sir.est[t] = sum(v.weight*particles[t,]);
    
    ## Compute the ESS %%%%%%%%%%%%%%%
    n.thr = 0.25*N
    n.eff = 1/(sum(v.weight^2))
    
    ## Resample and reset the weights %%%%%%%%%%%%%%% 
    if(n.eff < n.thr) { 
      index <- resample(v.weight,1:N)  
      index <- as.numeric(index$weights.pred)
      #particles <- particles[,index]
      particles <- matrix(index,             # Duplicate vector in matrix rows
                          nrow = n,
                          ncol = length(index),
                          byrow = TRUE)
      v.weight = rep(1,N)/N
    } 
    
  }
  
  for(l in 1:n){
    if(sir.est[l] < 0) {
      sir.est[l] <- 0
    } else {
      sir.est[l] <- sir.est[l]
    }
  }
  return(list(sir = sir.est, obs = obs.vec))
}