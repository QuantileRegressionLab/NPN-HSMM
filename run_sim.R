rm(list = ls())
library(markovchain)
library(cluster)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(mclust)
source("MainFunctions.R")


# set the true model parameters
N <- c(300,500,1000)
S <- 2 # number of states
P <- 30 # number of response
K <- 1 # the intercept
M <- 30 # maximum time spent in each state

mu = list(matrix(0, P, K), matrix(0, P, K))

theta = array(0,dim=c(P,P,S))
covarianza = array(0,dim=c(P,P,S))
for(j in 1:P){
  for(k in 1:P){
    if(j == k){
      theta[j,k,1] = 1
      theta[j,k,2] = 1
    }
    if(abs(j-k) == 1){
      theta[j,k,1] = 0.4
    }
    if(abs(j-k) == 2){
      theta[j,k,2] = 0.4
    }
  }
}
covarianza[,,1] = solve(theta[,,1])
covarianza[,,2] = solve(theta[,,2])

sigma_sim <- list(covarianza[,,1], covarianza[,,2])
init_sim <- c(1,0) # initial distribution
gamma_sim <- matrix(c(0,1,1,0),S,S,byrow=TRUE) # transition probability matrix

dl <- 1
d <- lapply(1:dl, function(x) {list()}) # sojourn distributions
# d[[1]] <- list(c(0.2,0.3,0.5), c(0.1,0.2,0.7))
d[[1]] <- list(shift.poi(1:M,10,1), shift.poi(1:M,5,1))
d[[2]] <- list(shift.nb(1:M,8,0.5), shift.nb(1:M,4,0.6))
d[[3]] <- list(dgeom(1:M-1,0.2), dgeom(1:M-1,0.1))

sojourn = c("poisson", "nbinom", "geometric")

error <- c("n", "outliers")


grids = list(seq(7,17,length.out=100), seq(10,19,length.out=100), seq(16,25,length.out=100))

nsim = 300

#################### N300 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N300 + POISSON + MULTIVARIATE GAUSSIAN")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N300_dPoisson_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N500 + POISSON + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N500_dPoisson_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + POISSON + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N1000_dPoisson_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + POISSON + OUTLIERS #########################################

print("N300 + POISSON + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N300_dPoisson_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + OUTLIERS #########################################

print("N500 + POISSON + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N500_dPoisson_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + OUTLIERS #########################################

print("N1000 + POISSON + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_poi_1 = rep(0,nsim)
lambda_poi_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[min_FP_HSMM]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_poi_1,lambda_poi_2)
file_name = paste("N1000_dPoisson_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N300 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N300_dNegBin_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N500 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N500_dNegBin_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N1000 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N1000_dNegBin_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + OUTLIERS #########################################

print("N300 + NEG BIN + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N300_dNegBin_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + OUTLIERS #########################################

print("N500 + NEG BIN + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N500_dNegBin_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + OUTLIERS #########################################

print("N1000 + NEG BIN + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
lambda_NB_1 = rep(0,nsim)
lambda_NB_2 = rep(0,nsim)
p_NB_1 = rep(0,nsim)
p_NB_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[min_FP_HSMM]
  lambda_NB_2[sim] = aaa$sojourn$mu[max_FP_HSMM]
  p_NB_1[sim] = aaa$sojourn$prob[min_FP_HSMM]
  p_NB_2[sim] = aaa$sojourn$prob[max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,lambda_NB_1,lambda_NB_2,p_NB_1,p_NB_2)
file_name = paste("N1000_dNegBin_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N300 + GEOM + MULTIVARIATE GAUSSIAN")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N300_dGeom_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N500 + GEOM + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N500_dGeom_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + GEOM + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N1000_dGeom_eMVNorm", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + OUTLIERS ####################################################

print("N300 + GEOM + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N300_dGeom_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + OUTLIERS ####################################################

print("N500 + GEOM + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N500_dGeom_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + OUTLIERS ####################################################

print("N1000 + GEOM + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0,nsim)
t.time_HMM = rep(0,nsim)
t.iter_HSMM = rep(0,nsim)
t.iter_HMM = rep(0,nsim)
lambda_HSMM = rep(0,nsim)
lambda_HMM = rep(0,nsim)
ARI_HSMM = rep(0,nsim)
ARI_HMM = rep(0,nsim)
errorRate_HSMM = rep(0,nsim)
errorRate_HMM = rep(0,nsim)
M_1_HSMM = rep(0,nsim)
M_2_HSMM = rep(0,nsim)
M_1_HMM = rep(0,nsim)
M_2_HMM = rep(0,nsim)
FP_1_HSMM = rep(0,nsim)
FP_2_HSMM = rep(0,nsim)
FP_1_HMM = rep(0,nsim)
FP_2_HMM = rep(0,nsim)
TFP_1_HSMM = rep(0,nsim)
TFP_2_HSMM = rep(0,nsim)
TFP_1_HMM = rep(0,nsim)
TFP_2_HMM = rep(0,nsim)
FN_1_HSMM = rep(0,nsim)
FN_2_HSMM = rep(0,nsim)
FN_1_HMM = rep(0,nsim)
FN_2_HMM = rep(0,nsim)
TFN_1_HSMM = rep(0,nsim)
TFN_2_HSMM = rep(0,nsim)
TFN_1_HMM = rep(0,nsim)
TFN_2_HMM = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
for(sim in 1:nsim){
  print(sim)
  set.seed(sim)
  
  data_gen <- hsmm.multi.gen(ns = N[a], P = P, K = K, m = S, delta = init_sim, gamma = gamma_sim, 
                             mu = mu, rho = sigma_sim, d = d[[b]], error = error[c])
  Y = data_gen$series
  state = data_gen$state
  states_init = pam(x = Y, k = S)
  states_init$clustering = sample(S, N[a], replace = T)
  # fit the underlying Markov chain
  hmm_init = markovchainFit(states_init$cluster)
  init = rep(0, S)
  init[states_init$cluster[1]] = 1
  # init[state[1]] = 1
  A = hmm_init$estimate@transitionMatrix
  A = A[-which(A %in% diag(A))]
  gamma_HMM = hmm_init$estimate@transitionMatrix
  gamma_HSMM = matrix(c(0,1,1,0),S,S,byrow=TRUE)
  
  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
    sigma = replicate(S, diag(P))
  )
  
  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_HSMM = rep(0,100)
  ICL_HMM = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa = try(EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = i, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
    if (inherits(bbb, "try-error")) {
      ICL_HMM[j] = Inf
    } else {
      ICL_HMM[j] = bbb$ICL
    }
  }
  minim_HSMM = which.min(ICL_HSMM)
  minim_HMM = which.min(ICL_HMM)
  lambda_HSMM[sim] = grid[minim_HSMM]
  lambda_HMM[sim] = grid[minim_HMM]
  aaa = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  bbb = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while(ARI_HSMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M-1,.1), dgeom(1:M-1,.1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(Y = Y, S = 2, Par = Par_HSMM, M = M, sojourn.distribution = sojourn[b], lambda = lambda_HSMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while(ARI_HMM[sim] < 0.3 & iter < 10){
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    gamma_HMM = hmm_init$estimate@transitionMatrix
    Par_HMM = list(
      gamma = gamma_HMM,
      init = init,
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HMM(Y = Y, S = 2, Par = Par_HMM, lambda = lambda_HMM[sim], pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2)
    if(ccc$loglik > llk){
      bbb = ccc
      llk = ccc$loglik
      ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
    }
  }
  t.time_HSMM[sim] = aaa$time
  t.time_HMM[sim] = bbb$time
  t.iter_HSMM[sim] = aaa$iter
  t.iter_HMM[sim] = bbb$iter
  errorRate_HSMM[sim] = classError(apply(aaa$u, 1, which.max), state)$errorRate
  errorRate_HMM[sim] = classError(apply(bbb$u, 1, which.max), state)$errorRate
  termine1 = sum(aaa$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(aaa$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(aaa$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(aaa$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(aaa$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HSMM = which.min(c(termine1, termine2))
  max_FP_HSMM = which.max(c(termine1, termine2))
  termine1 = sum(bbb$omega[,,1] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,2] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,2] != 0) / 2
  termine2 = sum(bbb$omega[,,2] != 0 & theta[,,1] == 0) / 2 + sum(bbb$omega[,,1] != 0 & theta[,,2] == 0) / 2 + sum(bbb$omega[,,2] == 0 & theta[,,1] != 0) / 2 + sum(bbb$omega[,,1] == 0 & theta[,,2] != 0) / 2
  min_FP_HMM = which.min(c(termine1, termine2))
  max_FP_HMM = which.max(c(termine1, termine2))
  FP_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] != 0 & theta[,,2] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] != 0 & theta[,,1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] != 0 & theta[,,2] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P*(P-1)/2 - P + 2)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P*(P-1)/2 - P + 2)
  FN_1_HSMM[sim] = sum(aaa$omega[,,min_FP_HSMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,,max_FP_HSMM] == 0 & theta[,,2] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,,min_FP_HMM] == 0 & theta[,,1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,,max_FP_HMM] == 0 & theta[,,2] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P-1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P-2)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P-1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P-2)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,,min_FP_HSMM]) - theta[,,1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,,max_FP_HSMM]) - theta[,,2])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,,min_FP_HMM]) - theta[,,1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,,max_FP_HMM]) - theta[,,2])^2)
  p_geom_1[sim] = aaa$d[1,min_FP_HSMM]
  p_geom_2[sim] = aaa$d[1,max_FP_HSMM]
}

matrix = cbind(t.time_HSMM,t.time_HMM,t.iter_HSMM,t.iter_HMM,lambda_HSMM,lambda_HMM,ARI_HSMM,ARI_HMM,errorRate_HSMM,errorRate_HMM,M_1_HSMM,M_1_HMM,M_2_HSMM,M_2_HMM,FP_1_HSMM,FP_1_HMM,FP_2_HSMM,FP_2_HMM,TFP_1_HSMM,TFP_1_HMM,TFP_2_HSMM,TFP_2_HMM,FN_1_HSMM,FN_1_HMM,FN_2_HSMM,FN_2_HMM,TFN_1_HSMM,TFN_1_HMM,TFN_2_HSMM,TFN_2_HMM,p_geom_1,p_geom_2)
file_name = paste("N1000_dGeom_eOutliers", ".csv", sep = "")
write.csv(matrix, file = file_name)