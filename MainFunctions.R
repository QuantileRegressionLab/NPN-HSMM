sd.trim <- function(x, trim=0, na.rm=FALSE, ...)
{
  if(!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    warning("argument is not numeric or logical: returning NA")
    return(NA_real_)
  }
  if(na.rm) x <- x[!is.na(x)]
  if(!is.numeric(trim) || length(trim) != 1)
    stop("'trim' must be numeric of length one")
  n <- length(x)
  if(trim > 0 && n > 0) {
    if(is.complex(x)) stop("trimmed sd are not defined for complex data")
    if(trim >= 0.5) return(0)
    lo <- floor(n * trim) + 1
    hi <- n + 1 - lo
    x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
  }
  sd(x)
}

cov2cor <- function(V)
{
  ## Purpose: Covariance matrix |--> Correlation matrix -- efficiently
  ## ----------------------------------------------------------------------
  ## Arguments: V: a covariance matrix (i.e. symmetric and positive definite)
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 12 Jun 2003, 11L:50
  p <- (d <- dim(V))[1L]
  if(!is.numeric(V) || length(d) != 2L || p != d[2L])
    stop("'V' is not a square numeric matrix")
  Is <- sqrt(1/diag(V)) # diag( 1/sigma_i )
  if(any(!is.finite(Is)))
    warning("diag(.) had 0 or NA entries; non-finite result is doubtful")
  r <- V # keep dimnames
  r[] <- Is * V * rep(Is, each = p)
  ##	== D %*% V %*% D  where D = diag(Is)
  r[cbind(1L:p,1L:p)] <- 1 # exact in diagonal
  r
}

##########################################################################################

#Sojourn distributions
shift.poi <- function(x,lambda,zeta){
  return(dshift.poi <- exp(-lambda)*lambda^(x-zeta)/factorial(x-zeta))
}

.dpois.hsmm.sojourn <- function(x=NULL,lambda,shift,log=FALSE)   {
  if(shift<0) stop(".dpois.hsmm.sojourn: shift must be > 0")
  if(log) dpois(x-shift,lambda,log=TRUE)
  else dpois(x-shift,lambda)
}

shift.nb <- function(x,k,p){
  dshift.nb <- gamma(x+k-1)/(factorial(x-1)*gamma(k))*p^k*(1-p)^(x-1)
  return(dshift.nb)
}

.fitnbinom <- function(eta) {
  shiftthresh=1e-20
  maxshift =  match(TRUE,eta>shiftthresh)
  Mtmp = tail(which(eta>shiftthresh),1)
  
  fun1 <- function(shift) {
    m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
    v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
    size <- if (v > m) m^2/(v - m) else 100
    densfun <- function(par) sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])
    optim(c(size,m),densfun,control=list(fnscale=-1))$value
  }
  
  shift = which.max(sapply(1:maxshift,fun1))
  m <- weighted.mean((maxshift:Mtmp)-shift,eta[maxshift:Mtmp])
  v <- as.numeric(cov.wt(data.frame((maxshift:Mtmp)-shift),wt=eta[maxshift:Mtmp])$cov)
  size <- if (v > m) m^2/(v - m) else 100
  densfun <- function(par) sum(dnbinom((maxshift:Mtmp)-shift,size=par[1],mu=par[2],log=TRUE)*eta[maxshift:Mtmp])
  tmp = optim(c(size,m),densfun,control=list(fnscale=-1))$par
  c(shift = shift,size=tmp[1],mu=tmp[2],prob=tmp[1]/(sum(tmp)))
}

getInitNPNParams  <- function(x, K) {
  D <- ncol(x)
  kMeansRes <- kmeans(x=x, centers=K)
  assignation <- kMeansRes$cluster
  sigma <- array(NA, dim=c(D, D, K))
  for(k in 1:K) {
    indices <- which(assignation==k)
    sigma[,,k] <- cov(x[indices,])
  }
  phi <- sigma
  return(phi)
}

F_tilde = function(x){
  N = nrow(x)
  D = ncol(x)
  F = apply(x,2,ecdf)
  delta_N = 1/(4*N^(1/4)*sqrt(pi)*log(N))
  rip = matrix(NA,N,D)
  for(i in 1:D){
    rip[,i] = F[[i]](x[,i])
  }
  rip[which(rip < delta_N)] = delta_N
  rip[which(rip > (1 - delta_N))] = 1 - delta_N
  return(rip)
}

NPNdensity = function(x, mu, sigma, g, gp){
  D = length(x)
  f = 1 / ((2*pi)^(D/2) * det(sigma)^(1/2)) *
    exp(-1/2 * t(g - mu) %*% solve(sigma) %*% (g - mu)) * prod(abs(gp))
  return(f)
}

getNPNProbabilities <- function(x,phi) {
  N <- ncol(x)
  D <- nrow(x)
  K <- dim(phi$sigma)[3]
  rip <- F_tilde(x)
  rip <- rip*N/(N+1)
  # rip[which(rip == 0.5)] = ...
  g_tilde = qnorm(rip)
  p <- matrix(NA, nrow=N, ncol=K)
  for(n in 1:N) {
    for(k in 1:K) {
      p[n,k] <- NPNdensity(x=x[,n], sigma=phi$sigma[,,k], g=g_tilde[,n])
    }
  }
  return(p)
}

getUpdatedNPNParams <- function(x, mu, u, g_tilde, K, penalty) {
  N = nrow(x)
  D = ncol(x)
  sigma = array(0, dim = c(D, D, K))
  omega = array(0, dim = c(D, D, K))
  S_tilde = array(0, dim = c(D, D, K))
  loglik = rep(0,K)
  for(k in 1:K){
    num = matrix(0, D, D)
    for(i in 1:N){
      g = g_tilde[i,]
      temp = u[i,k] * (g - mu) %*% t(g - mu)
      num = num + temp
    }
    S_tilde[,,k] = num / sum(u[,k])
  }
  for(k in 1:K){
    obj = glasso(S_tilde[,,k], penalty, nobs = N, penalize.diagonal = F)
    sigma[,,k] = obj$w
    omega[,,k] = obj$wi
  }
  phi <- list(sigma=sigma, omega=omega)
  return(phi)
}

# Simulate data from a multiple HSMM regression model
hsmm.multi.gen <- function(ns, P, K, m, delta, gamma, mu, sigma, rho, d, error){
  # ns: sample size (scalar)
  # P: # of response variables (scalar)
  # K: 1 (scalar)
  # m: # of states (scalar)
  # delta: vector of prior probabilities
  # gamma: matrix of state transition probabilities
  # mu: list of state-specific parameters
  # sigma: list of state-specific marginal variances
  # rho: list of state-dependent correlation matrices
  # d: list of state duration probabilities
  # error: error distribution (Gaussian, contaminated Gaussian, exp. transformed Gaussian, Student-t with 5 d.f. or Pearson)
  ld <- sapply(d,length)
  mvect <- 1:m
  state <- numeric(ns)
  x <- matrix(0, ns, P)
  alpha = 0.1
  beta_n = 0.05
  
  state[1] <- sample(mvect, 1, prob=delta)
  dur <- sample(1:ld[state[1]],1,prob=d[[state[1]]])
  
  if (error == "n") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      e[t,]=mvrnorm(1, mu = mu[[state[t]]], Sigma = rho[[state[t]]])
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        e[total+t,]=mvrnorm(1, mu = mu[[state[total+t]]], Sigma = rho[[state[total+t]]])
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
  } else if (error == "outliers") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      e[t,]=mvrnorm(1, mu = mu[[state[t]]], Sigma = rho[[state[t]]])
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        e[total+t,]=mvrnorm(1, mu = mu[[state[total+t]]], Sigma = rho[[state[total+t]]])
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
    n_substitutes = round(beta_n*ns)
    uniforms = matrix(runif(n_substitutes*P,-5,5),n_substitutes,P)
    for(iter in 1:P){
      to_substitute = sample(1:ns,n_substitutes)
      x[to_substitute,iter] = uniforms[,iter]
    }
  } else if (error == "Non-Gaussian") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      e[t,]=mvrnorm(1, mu = mu[[state[t]]], Sigma = rho[[state[t]]])
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        e[total+t,]=mvrnorm(1, mu = mu[[state[total+t]]], Sigma = rho[[state[total+t]]])
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
    n_substitutes = round(alpha*P)
    x[,1:n_substitutes] = exp(x[,1:n_substitutes])
  } else if (error == "t") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      e[t,]=mvtnorm::rmvt(1, df = 5, delta = mu[[state[t]]], sigma = rho[[state[t]]])
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        e[total+t,]=mvtnorm::rmvt(1, df = 5, delta = mu[[state[total+t]]], sigma = rho[[state[total+t]]])
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
  } else if (error == "norm+t") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      normal_data = mvrnorm(1, mu = mu[[state[t]]], Sigma = rho[[state[t]]])
      cdf = array(NA,length(normal_data))
      ppf = array(NA,length(normal_data))
      for(p in 1:length(normal_data)){
        cdf[p] = pnorm(normal_data[p], mean = mu[[state[t]]][p], sd = sqrt(rho[[state[t]]][p,p]))
        ppf[p] = qt(cdf[p], df = 5)
      }
      e[t,] = ppf
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        normal_data = mvrnorm(1, mu = mu[[state[total+t]]], Sigma = rho[[state[total+t]]])
        cdf = array(NA,length(normal_data))
        ppf = array(NA,length(normal_data))
        for(p in 1:length(normal_data)){
          cdf[p] = pnorm(normal_data[p], mean = mu[[state[total+t]]][p], sd = sqrt(rho[[state[total+t]]][p,p]))
          ppf[p] = qt(cdf[p], df = 5)
        }
        e[total+t,] = ppf
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
  } else if (error == "norm+pears") {
    
    e=matrix(NA, ns, P)
    for(t in 1:dur){
      state[t] <- state[1]
      normal_data = mvrnorm(1, mu = mu[[state[t]]], Sigma = rho[[state[t]]])
      cdf = array(NA,length(normal_data))
      ppf = array(NA,length(normal_data))
      for(p in 1:length(normal_data)){
        cdf[p] = pnorm(normal_data[p], mean = mu[[state[t]]][p], sd = sqrt(rho[[state[t]]][p,p]))
        ppf[p] = qpearsonIV(cdf[p], 3, 0.3, 0.003, 0.04)
      }
      e[t,] = ppf
      
      x[t,] <- e[t,]
    }
    
    total <- dur
    
    while(total<ns){
      state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
      dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
      for(t in 1:dur){
        if(total+t>ns) break
        state[total+t] <- state[total+1]
        normal_data = mvrnorm(1, mu = mu[[state[total+t]]], Sigma = rho[[state[total+t]]])
        cdf = array(NA,length(normal_data))
        ppf = array(NA,length(normal_data))
        for(p in 1:length(normal_data)){
          cdf[p] = pnorm(normal_data[p], mean = mu[[state[total+t]]][p], sd = sqrt(rho[[state[total+t]]][p,p]))
          ppf[p] = qpearsonIV(cdf[p], 3, 0.3, 0.003, 0.04)
        }
        e[total+t,] = ppf
        
        x[total+t,] <- e[total+t,]
      }
      total <- total + dur
    }
  }
  
  return(list(series=x, state=state, error=e))
}


# EM algorithm for Hidden Semi Markov Model
EM_HSMM = function(Y, S, Par = NULL, M, sojourn.distribution, lambda, pen_EBIC, seed, err, iterMax){
  # Y: P-variate response variable (matrix)
  # S: # of states
  # Par: initial model parameters (list)
  # M: maximum time spent in each state (scalar)
  # err: convergence threshold
  # iterMax: maximum # of iterations
  
  Y = apply(Y, 2, scale)  # normal scores
  
  N <- nrow(Y)
  P <- ncol(Y)
  M.max = max(M)
  
  if(is.null(Par)) {
    set.seed(seed)
    #mu = replicate(S, colMeans(Y))
    sigma = replicate(S, cor(Y))
    gamma = matrix(c(0,1,1,0),S,S,byrow=TRUE)
    init = c(1, rep(0, S - 1))
    if(S == 1) {
      d = 1
    } else {
      d = matrix(dunif(1:M, 1, M), M, S)
    }
  } else {
    #mu = Par$mu
    sigma = Par$sigma
    gamma = Par$gamma
    init = Par$init
    d = Par$d
  }
  
  dif <- Inf
  t.iter <- 0
  llk.old <- -10^250
  
  rip <- F_tilde(Y)
  rip <- rip*N/(N+1)
  #mean.emp = apply(Y, 2, mean)
  #sd.emp = apply(Y, 2, sd)
  #g_tilde = qnorm(rip, mean = mean.emp, sd = sd.emp)
  g_tilde = qnorm(rip)
  # g_tilde = Y
  dY <- sapply(1:P, function(j) density(Y[,j]))
  dYout <- sapply(1:P, function(j) approx(dY[,j]$x, dY[,j]$y, xout = Y[,j])$y)
  gprime_tilde = dYout/dnorm(qnorm(rip))
  abs_gprime_tilde = apply(abs(gprime_tilde), 1, prod)
  
  t0 = Sys.time()
  while (dif > err && t.iter < iterMax || t.iter < 5) {
    t.iter=t.iter+1
    
    fden <- matrix(NA, nrow = N, ncol = S)
    for (s in 1:S) {
      fden[,s] <- dmvnorm(data = g_tilde, sigma = sigma[,,s]) * abs_gprime_tilde
    }
    fden[is.na(fden) | is.infinite(fden)] = 1e-10
    
    d=head(d,M.max)
    d=matrix(d, M.max, S)
    dD = apply(d,2,function(x) rev(cumsum(rev(x))))
    estep_variables  = .C("backward",
                          transition=as.double(gamma),
                          init=as.double(init),
                          p=as.double(fden),
                          d=as.double(d),
                          D=as.double(dD),
                          timelength=as.integer(N),
                          J=as.integer(S),
                          M=as.integer(rep(M.max, S)),
                          L1=double(NROW(Y)*S),N=double(NROW(Y)),
                          eta=double(M.max*S),
                          F1=double(S*NROW(Y)),
                          si=double(S*NROW(Y)),
                          gamma=double(S*NROW(Y)),
                          nsequences=as.integer(length(N)),
                          totallength=NROW(Y),
                          G=double(S*NROW(Y)),
                          PACKAGE='mhsmm')
    
    if(S == 1) {
      estep_variables$gamma = rep(1, N)
      estep_variables$eta = N
      estep_variables$init = 1
      estep_variables$transition = 1
    } else {
      if(any(estep_variables$gamma<0)) estep_variables$gamma = zapsmall(estep_variables$gamma)
      if(any(estep_variables$eta<0)) estep_variables$eta = zapsmall(estep_variables$eta)
      if(any(estep_variables$N<0))  estep_variables$N = zapsmall(estep_variables$N)
    }
    
    u = matrix(estep_variables$gamma,ncol=S)

    init.new = estep_variables$init
    init.new[init.new < 0] = 0
    gamma.new = matrix(estep_variables$transition,ncol=S)
    gamma.new[gamma.new < 0] = 0
    
    # mu.new = sapply(1:S, function(s) colMeans(g_tilde * u[,s]))
    #mu.new = matrix(0, P, S)
    
    S_tilde = sigma.new = omega.new = array(NA, dim = c(P, P, S))
    for (s in 1:S) {
      S_tilde[,,s] = wcrossprod(x = g_tilde, y = g_tilde, w = u[,s]) / sum(u[,s])
      penalty = 2*lambda / sum(u[,s])
      obj = glasso(S_tilde[,,s], rho = penalty, nobs = N, penalize.diagonal = F, maxit = 1e5)
      sigma.new[,,s] = obj$w
      omega.new[,,s] = obj$wi
    }
    
    d.new = apply(matrix(estep_variables$eta,ncol=S),2,function(x) x/sum(x))
    d.new = zapsmall(d.new)
    
    # ksmoothed-nonparametric
    # ksmooth.thresh = 1e-20
    # d.new = apply(matrix(estep_variables$eta+1e-100,ncol=S),2,function(x) x/sum(x))
    # for(s in 1:S) {
    #   d.new[,s] = approx(density(which(d.new[,s]>ksmooth.thresh),weights=d.new[which(d.new[,s]>ksmooth.thresh),s],from=1,n=M),xout=1:M)$y
    #   d.new[is.na(d.new[,s]),s] = 0
    #   d.new[,s] = (d.new[,s]+1e-300)/sum(d.new[,s])
    # }
    
    llk = sum(log(estep_variables$N), na.rm = T)
    pen_llk = sum(log(estep_variables$N), na.rm = T) - lambda * sum(abs(omega.new))
    dif <- abs(pen_llk-llk.old)

    # print(round(c(t.iter, dif, llk), 3))
    
    gamma=gamma.new
    init=init.new
    #mu=mu.new
    sigma=sigma.new
    omega=omega.new
    d=d.new
    llk.old=pen_llk
  }
  c.time = Sys.time() - t0
  
  omegaT = sapply(1:S, function(s) 2 * omega[,,s] - omega[,,s]%*%S_tilde[,,s]%*%omega[,,s], simplify = "array")
  
  predicted_state = apply(u, 1, which.max)
  post.hmm = u
  E = (sum(omega != 0) - P * S) / 2
  n_parameters = E + P*S + S*(S-2) + sum(d != 0) - S
  aic.crit = -2*llk + 2*n_parameters
  bic.crit = -2*llk + log(N)*n_parameters
  EBIC = bic.crit + 4*pen_EBIC*n_parameters*log(P)
  icl.crit = bic.crit - 2 * sum(post.hmm * ifelse(post.hmm > 0, log(post.hmm), 0))
  
  shiftthresh = 1e-20
  sojourn = list()
  if(sojourn.distribution=="poisson") {
    sojourn$lambda_poi = S
    sojourn$shift = S
    for(i in 1:S){
      eta = d[,i]
      maxshift =  match(TRUE,eta>shiftthresh)
      Mtmp = tail(which(eta>shiftthresh),1)
      sojourn$shift[i] = which.max(sapply(1:maxshift, function(shift) .dpois.hsmm.sojourn(x = maxshift:Mtmp,lambda=((maxshift:Mtmp)-shift)%*%eta[maxshift:Mtmp],shift=shift,log=TRUE)%*%eta[maxshift:Mtmp]))
      sojourn$lambda_poi[i] = ((sojourn$shift[i]:Mtmp)-sojourn$shift[i])%*%eta[sojourn$shift[i]:Mtmp]
    }
  }
  
  else if(sojourn.distribution=="nbinom") {
    sojourn$size = numeric(S)
    sojourn$shift = integer(S)
    sojourn$mu = numeric(S)
    sojourn$prob = numeric(S)
    eta = matrix(estep_variables$eta,ncol=S)
    for(i in 1:S) {
      tryCatch({
        tmp = .fitnbinom(eta[,i])
        sojourn$shift[i] = tmp[1]
        sojourn$size[i] =  tmp[2]
        sojourn$mu[i] =  tmp[3]
        sojourn$prob[i] =  tmp[4]
      }, error = function(e){
        sojourn$shift[i] = NA
        sojourn$size[i] =  NA
        sojourn$mu[i] =  NA
        sojourn$prob[i] =  NA
      })
    }
  }
    
  
  lout = list()
  lout$S = S
  lout$M = M
  lout$iter = t.iter
  lout$dif = dif
  lout$loglik = llk
  lout$penloglik = pen_llk
  lout$gamma = gamma
  #lout$mu = mu
  lout$S_tilde = S_tilde 
  lout$sigma = sigma
  lout$omega = omega
  lout$omegaT = omegaT
  lout$d = d
  lout$u = u
  lout$init = init
  lout$time = c.time
  lout$emission_prob = fden
  lout$predicted_state = predicted_state
  lout$E = E
  lout$n_parameters = n_parameters
  lout$AIC = aic.crit
  lout$BIC = bic.crit
  lout$EBIC = EBIC
  lout$ICL = icl.crit
  lout$sojourn = sojourn
  
  return(lout)
  
}


# EM algorithm for Hidden Markov Model
EM_HMM = function(Y, S, Par = NULL, lambda, pen_EBIC, seed, err, iterMax){
  # Y: P-variate response variable (matrix)
  # S: # of states
  # Par: initial model parameters (list)
  # err: convergence threshold
  # iterMax: maximum # of iterations
  
  Y = apply(Y, 2, scale)  # normal scores
  
  N <- nrow(Y)
  P <- ncol(Y)
  
  if(is.null(Par)) {
    set.seed(seed)
    #mu = replicate(S, colMeans(Y))
    sigma = replicate(S, cor(Y))
  } else {
    #mu = Par$mu
    sigma = Par$sigma
    gamma = Par$gamma
    init = Par$init
  }
  
  dif <- Inf
  t.iter <- 0
  llk.old <- -10^250
  u = double(S*sum(N))
  
  rip <- F_tilde(Y)
  rip <- rip*N/(N+1)
  #mean.emp = apply(Y, 2, mean)
  #sd.emp = apply(Y, 2, sd)
  #g_tilde = qnorm(rip, mean = mean.emp, sd = sd.emp)
  g_tilde = qnorm(rip)
  # g_tilde = Y
  dY <- sapply(1:P, function(j) density(Y[,j]))
  dYout <- sapply(1:P, function(j) approx(dY[,j]$x, dY[,j]$y, xout = Y[,j])$y)
  gprime_tilde = dYout/dnorm(qnorm(rip))
  abs_gprime_tilde = apply(abs(gprime_tilde), 1, prod)
  
  t0 = Sys.time()
  while (dif > err && t.iter < iterMax || t.iter < 5) {
    t.iter=t.iter+1
    
    # E Step
    fden <- matrix(NA, nrow = N, ncol = S)
    for (s in 1:S) {
      fden[,s] <- dmvnorm(data = g_tilde, sigma = sigma[,,s]) * abs_gprime_tilde
    }
    fden[is.na(fden) | is.infinite(fden)] = 1e-10
    
    estep_out = .C("mo_estep_hmm",a=as.double(t(gamma)),pi=as.double(t(init)),p=as.double(t(fden)),
                   N=as.integer(N),nsequences=as.integer(length(N)), K=as.integer(S),
                   alpha=double((S+1)*sum(N)) ,beta=double(S*sum(N)),gam=u,ll=double(1),PACKAGE='mhsmm')
    
    u = matrix(estep_out$gam,ncol=S)
    
    # M Step
    init.new = estep_out$pi
    init.new[init.new < 0] = 0
    gamma.new = matrix(estep_out$a,nrow=S,byrow=T)
    gamma.new[gamma.new < 0] = 0
    
    # mu.new = sapply(1:S, function(s) colMeans(g_tilde * u[,s]))
    #mu.new = matrix(0, P, S)
    
    S_tilde = sigma.new = omega.new = array(NA, dim = c(P, P, S))
    for (s in 1:S) {
      S_tilde[,,s] = wcrossprod(x = g_tilde, y = g_tilde, w = u[,s]) / sum(u[,s])
      penalty = 2*lambda / sum(u[,s])
      obj = glasso(S_tilde[,,s], rho = penalty, nobs = N, penalize.diagonal = F, maxit = 1e5)
      sigma.new[,,s] = obj$w
      omega.new[,,s] = obj$wi
    }
    
    llk = estep_out$ll
    pen_llk = llk - lambda * sum(abs(omega.new))
    dif <- abs(pen_llk-llk.old)
    

    gamma=gamma.new
    init=init.new
    #mu=mu.new
    sigma=sigma.new
    omega=omega.new
    llk.old=pen_llk
  }
  c.time = Sys.time() - t0
  
  predicted_state = apply(u, 1, which.max)
  post.hmm = u
  E = (sum(omega != 0) - P * S) / 2
  n_parameters = E + P*S + S*(S-1)
  aic.crit = -2*llk + 2*n_parameters
  bic.crit = -2*llk + log(N)*n_parameters
  EBIC = bic.crit + 4*pen_EBIC*n_parameters*log(P)
  icl.crit = bic.crit - 2 * sum(post.hmm * ifelse(post.hmm > 0, log(post.hmm), 0))
  
  lout = list()
  lout$S = S
  lout$iter = t.iter
  lout$dif = dif
  lout$loglik = llk
  lout$penloglik = pen_llk
  lout$gamma = gamma
  #lout$mu = mu
  lout$sigma = sigma
  lout$omega = omega
  lout$u = u
  lout$time = c.time
  lout$predicted_state = predicted_state
  lout$E = E
  lout$n_parameters = n_parameters
  lout$AIC = aic.crit
  lout$BIC = bic.crit
  lout$EBIC = EBIC
  lout$ICL = icl.crit
  
  return(lout)
  
}

Viterbi=function(y,transProbs,emissionProbs,initial_distribution) {
  
  T = length(y)
  M = nrow(transProbs)
  prev = matrix(0, T-1, M)
  omega = matrix(0, M, T)
  
  omega[, 1] = log(initial_distribution * emissionProbs[1, ])
  for(t in 2:T){
    for(s in 1:M) {
      probs = omega[, t - 1] + log(transProbs[, s]) + log(emissionProbs[t,s])
      prev[t - 1, s] = which.max(probs)
      omega[s, t] = max(probs)
    }
  }
  
  S = rep(0, T)
  last_state=which.max(omega[,ncol(omega)])
  S[1]=last_state
  
  j=2
  for(i in (T-1):1){
    S[j]=prev[i,last_state] 
    last_state=prev[i,last_state] 
    j=j+1
  }
  
  
  S=rev(S)
  
  return(S)
  
}

hsmm2hmm<-function(omega,dm,eps=1e-10){
  mv<-sapply(dm,length)
  m<-length(mv)
  G<-matrix(0,0,sum(mv))
  for (i in 1:m){
    mi<-mv[[i]]
    F<-cumsum(c(0,dm[[i]][-mi]))
    ci<-ifelse(abs(1-F)>eps,dm[[i]]/(1-F),1)
    cim<-ifelse(1-ci>0,1-ci,0)
    Gi<-matrix(0,mi,0)
    for (j in 1:m){
      if(i==j) { if(mi==1)
      { Gi<-cbind(Gi,c(rep(0,mv[[j]]-1),cim))} else
      { Gi<-cbind(Gi,rbind(cbind(rep(0,mi-1),diag(cim[-mi],mi-1,mi-1)),
                           c(rep(0,mi-1),cim[[mi]])))}
      } else   { if(mi==1)
      { Gi<-cbind(Gi,matrix(c(omega[[i,j]]*ci,rep(0,mv[[j]]-1)),1))} else
      { Gi<-cbind(Gi,cbind(omega[[i,j]]*ci,matrix(0,mv[[i]],mv[[j]]-1)))}
      }
    }
    G<-rbind(G,Gi)
  }
  G }


# invF <- lapply(seq_len(ncol(Y)), function(j)
#   function(u) quantile(Y[, j], probs = u, type = 8))

rNPN <- function(n, rho){
  rho = cov2cor(rho)
  Z <- rmvnorm(n, sigma = rho)
  U <- pnorm(Z)
  X <- qnorm(U) 
  # sapply(seq_along(invF), function(j) invF[[j]](U[, j]))
}

# Resample from a fitted HSMM nonparanormal (parametric bootstrap)
boot.hsmm.multi.gen <- function(ns, P, K, m, delta, gamma, rho, d){
  # ns: sample size (scalar)
  # P: # of response variables (scalar)
  # K: 1 (scalar)
  # m: # of states (scalar)
  # delta: vector of prior probabilities
  # gamma: matrix of state transition probabilities
  # mu: list of state-specific parameters
  # sigma: list of state-specific marginal variances
  # rho: list of state-dependent correlation matrices
  # d: list of state duration probabilities
  # error: error distribution
  ld <- sapply(d,length)
  mvect <- 1:m
  state <- numeric(ns)
  x <- matrix(0, ns, P)
  mu <- replicate(m, matrix(0, P, K), simplify = F)
  
  state[1] <- sample(mvect, 1, prob=delta)
  dur <- sample(1:ld[state[1]],1,prob=d[[state[1]]])
  
  e=matrix(NA, ns, P)
  for(t in 1:dur){
    state[t] <- state[1]
    e[t,]=rNPN(1, rho = rho[[state[t]]])
    
    x[t,] <- e[t,]
  }
  
  total <- dur
  
  while(total<ns){
    state[total+1] <- sample(mvect,1,prob=gamma[state[total],])
    dur <- sample(1:ld[state[total+1]],1,prob=d[[state[total+1]]])
    for(t in 1:dur){
      if(total+t>ns) break
      state[total+t] <- state[total+1]
      e[total+t,]=rNPN(1, rho = rho[[state[total+t]]])
      
      x[total+t,] <- e[total+t,]
    }
    total <- total + dur
  }
  
  return(list(series=x, state=state, error=e))
}

coverage.summary <- function(fit, boot, theta, d_true, gamma_sim, sojourn.distribution, alpha = 0.05) {
  
  P = dim(fit[[1]]$omegaT)[1]
  S = fit[[1]]$S
  nsim = length(sapply(fit, length))
  R = length(sapply(boot[[1]], length))
  M = fit[[1]]$M
  if(sojourn.distribution == "poisson") b = 1
  if(sojourn.distribution == "nbinom") b = 2
  
  covered_omega = array(NA, c(P, P, S, nsim))
  covered_d = array(NA, c(M, S, nsim))
  covered_sojourn = matrix(NA, S, nsim)
  covered_gamma = array(NA, c(S, S, nsim))
  qqnorm = qnorm(1-alpha/2) # qt(1-alpha/2, df = R)
  for (sim in 1:nsim) {
    Ra = which(sapply(1:R, function(r) !is.null(results_list.boot[[sim]][[r]])))
    
    for (r in 1:R) {
      if(!is.null(boot[[sim]][[r]])) {
        # order the bootstrap estimates
        state.order = order(apply(boot[[sim]][[r]]$omega[2:4,1,], 2, which.max))
        # state.order = order(apply(boot[[sim]][[r]]$d, 2, which.max), decreasing = T)
        if(sojourn.distribution == "poisson") state.order = order(boot[[sim]][[r]]$sojourn$lambda_poi, decreasing = T)
        if(sojourn.distribution == "nbinom") state.order = order(boot[[sim]][[r]]$sojourn$mu, decreasing = T)
        boot[[sim]][[r]]$d = boot[[sim]][[r]]$d[,state.order]
        boot[[sim]][[r]]$gamma = boot[[sim]][[r]]$gamma[state.order,state.order]
        boot[[sim]][[r]]$omegaT = boot[[sim]][[r]]$omegaT[,,state.order]
        boot[[sim]][[r]]$omegaT = sapply(1:S, function(s) cov2cor(boot[[sim]][[r]]$omegaT[,,s]), simplify = "array")
      }
    }
    omega.array = simplify2array(sapply(Ra, function(r) boot[[sim]][[r]]$omegaT, simplify = F))
    d.array = simplify2array(sapply(Ra, function(r) boot[[sim]][[r]]$d, simplify = F))
    # sojourn.array = simplify2array(sapply(Ra, function(r) boot[[sim]][[r]]$sojourn$lambda_poi, simplify = F))
    gamma.array = simplify2array(sapply(Ra, function(r) boot[[sim]][[r]]$gamma, simplify = F))
    
    omega.sd = apply(omega.array, 1:3, sd)
    d.sd = apply(d.array, 1:2, sd)
    gamma.sd = apply(gamma.array, 1:2, sd)
    
    # order the fit estimates
    state.order = order(apply(fit[[sim]]$omega[2:4,1,], 2, which.max))
    # state.order = order(apply(fit[[sim]]$d, 2, which.max), decreasing = T)
    if(sojourn.distribution == "poisson") state.order = order(fit[[sim]]$sojourn$lambda_poi, decreasing = T)
    if(sojourn.distribution == "nbinom") state.order = order(fit[[sim]]$sojourn$mu, decreasing = T)
    gamma_hat = fit[[sim]]$gamma[state.order,state.order]
    rho_hat = fit[[sim]]$omegaT[,,state.order]
    rho_hat = simplify2array(sapply(1:S, function(s) as.matrix(forceSymmetric(cov2cor(rho_hat[,,s]))), simplify = F))
    sigma_hat = fit[[sim]]$S_tilde[,,state.order]
    d_hat = fit[[sim]]$d[,state.order]
    d_hat = as.list(as.data.frame(d_hat))

    for (s in 1:S) {
      tmp = (theta[,,s] >= rho_hat[,,s] - qqnorm * omega.sd[,,s]) & (theta[,,s] <= rho_hat[,,s] + qqnorm * omega.sd[,,s])
      covered_omega[,,s,sim] = tmp
      
      tmp = (d_true[[b]][[s]] >= d_hat[[s]] - qqnorm * d.sd[,s]) & (d_true[[b]][[s]] <= d_hat[[s]] + qqnorm * d.sd[,s])
      covered_d[,s,sim] = tmp * ifelse(d.sd[,s] == 0, NA, 1)
    }
    
    tmp = (gamma_sim >= gamma_hat - qqnorm * gamma.sd) & (gamma_sim <= gamma_hat + qqnorm * gamma.sd)
    covered_gamma[,,sim] = tmp
  }
  
  coverage_gamma = apply(covered_gamma, 1:2, mean)
  coverage_d = sapply(1:S, function(s) mean(covered_d[,s,], na.rm = T))
  coverage_omega = sapply(1:S, function(s) mean(apply(covered_omega[,,s,], 1:2, mean), na.rm = T))
  
  omegaA = omegaAc = matrix(NA, S, 1)
  for (s in 1:S) {
    tmp = ifelse(theta[,,s] != 0, apply(covered_omega[,,s,], 1:2, mean), NA)
    omegaA[s] = mean(tmp[lower.tri(tmp)], na.rm = T)
    
    tmp = ifelse(theta[,,s] == 0, apply(covered_omega[,,s,], 1:2, mean), NA)
    omegaAc[s] = mean(tmp[lower.tri(tmp)], na.rm = T)
  }
  
  out = list()
  out$coverage_gamma = coverage_gamma
  out$coverage_d = coverage_d
  out$coverage_omega = coverage_omega
  out$coverage_omegaA = omegaA
  out$coverage_omegaAc = omegaAc
  
  return(out)
}

# ksmoothed-nonparametric of the sojourn distribution
dksmoothed = function(d, ksmooth.thresh = 1e-20) {
  S = ncol(d)
  M = nrow(d)
  d.new = apply(matrix(d+1e-100,ncol=S),2,function(x) x/sum(x))
  for(s in 1:S) {
    d.new[,s] = approx(density(which(d.new[,s]>ksmooth.thresh),weights=d.new[which(d.new[,s]>ksmooth.thresh),s],from=1,n=M),xout=1:M)$y
    d.new[is.na(d.new[,s]),s] = 0
    d.new[,s] = (d.new[,s]+1e-300)/sum(d.new[,s])
  }
  
  return(d.new)
}