rm(list = ls())
library(markovchain)
library(cluster)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(mclust)
library(PearsonDS)
source("MainFunctions.R")
source("em_glasso.R")


# set the true model parameters
N <- c(300,500,1000)
S <- 3 # number of states
P <- 30 # number of response
K <- 1 # the intercept
M <- 30 # maximum time spent in each state

mu = list(matrix(0, P, K), matrix(0, P, K), matrix(0, P, K))

theta = array(0,dim=c(P,P,S))
covarianza = array(0,dim=c(P,P,S))
for(j in 1:P){
  for(k in 1:P){
    if(j == k){
      theta[j,k,1] = 1
      theta[j,k,2] = 1
      theta[j,k,3] = 1
    }
    if(abs(j-k) == 1){
      theta[j,k,1] = 0.4
    }
    if(abs(j-k) == 2){
      theta[j,k,2] = 0.4
    }
    if(abs(j-k) == 3){
      theta[j,k,3] = 0.4
    }
  }
}
covarianza[,,1] = solve(theta[,,1])
covarianza[,,2] = solve(theta[,,2])
covarianza[,,3] = solve(theta[,,3])

sigma_sim <- list(covarianza[,,1], covarianza[,,2], covarianza[,,3])
init_sim <- c(1,0,0) # initial distribution
gamma_sim <- matrix(c(0,0.8,0.2,0.4,0,0.6,0.7,0.3,0),S,S,byrow=TRUE) # transition probability matrix
gamma_sim <- matrix(0.5,S,S)
diag(gamma_sim) <- 0

dl <- 1
d <- lapply(1:dl, function(x) {list()}) # sojourn distributions
# d[[1]] <- list(c(0.2,0.3,0.5), c(0.1,0.2,0.7))
d[[1]] <- list(shift.poi(1:M,15,1), shift.poi(1:M,10,1), shift.poi(1:M,5,1))
d[[2]] <- list(shift.nb(1:M,10,0.5), shift.nb(1:M,7,0.6), shift.nb(1:M,4,0.7))
d[[3]] <- list(dgeom(1:M-1,0.3), dgeom(1:M-1,0.15), dgeom(1:M-1,0.1))

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
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dPoisson_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N500 + POISSON + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dPoisson_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + POISSON + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dPoisson_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + POISSON + OUTLIERS #########################################

print("N300 + POISSON + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dPoisson_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + OUTLIERS #########################################

print("N500 + POISSON + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dPoisson_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + OUTLIERS #########################################

print("N1000 + POISSON + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dPoisson_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N300 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dNegBin_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N500 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dNegBin_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N1000 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dNegBin_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + OUTLIERS #########################################

print("N300 + NEG BIN + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dNegBin_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + OUTLIERS #########################################

print("N500 + NEG BIN + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dNegBin_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + OUTLIERS #########################################

print("N1000 + NEG BIN + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dNegBin_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N300 + GEOM + MULTIVARIATE GAUSSIAN")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dGeom_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N500 + GEOM + MULTIVARIATE GAUSSIAN")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dGeom_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + GEOM + MULTIVARIATE GAUSSIAN")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dGeom_eMVNorm_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + OUTLIERS ####################################################

print("N300 + GEOM + OUTLIERS")
a = 1  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N300_dGeom_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + OUTLIERS ####################################################

print("N500 + GEOM + OUTLIERS")
a = 2  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N500_dGeom_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + OUTLIERS ####################################################

print("N1000 + GEOM + OUTLIERS")
a = 3  # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3  # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2  # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_glasso = rep(0,nsim)
t.iter_glasso = rep(0,nsim)
lambda_glasso = rep(0,nsim)
ARI_glasso = rep(0,nsim)
errorRate_glasso = rep(0,nsim)
M_1_glasso = rep(0,nsim)
M_2_glasso = rep(0,nsim)
M_3_glasso = rep(0,nsim)
FP_1_glasso = rep(0,nsim)
FP_2_glasso = rep(0,nsim)
FP_3_glasso = rep(0,nsim)
TFP_1_glasso = rep(0,nsim)
TFP_2_glasso = rep(0,nsim)
TFP_3_glasso = rep(0,nsim)
FN_1_glasso = rep(0,nsim)
FN_2_glasso = rep(0,nsim)
FN_3_glasso = rep(0,nsim)
TFN_1_glasso = rep(0,nsim)
TFN_2_glasso = rep(0,nsim)
TFN_3_glasso = rep(0,nsim)
p_geom_1 = rep(0,nsim)
p_geom_2 = rep(0,nsim)
p_geom_3 = rep(0,nsim)
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
  
  Sigma.s = replicate(S, diag(P), simplify = F)
  mu.s = matrix(0, S, P)
  
  # FITTING DELLA SERIE STORICA SIMULATA
  
  ICL_glasso = rep(0,100)
  j = 0
  for(i in grid){
    j = j+1
    aaa =   result <- try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = i, err = 1e-4, iterMax = 5e2, traceEM = F))
    if (inherits(aaa, "try-error")) {
      ICL_glasso[j] = Inf
    } else {
      ICL_glasso[j] = aaa$pen.criteria$ICL
    }
  }
  
  minim_glasso = which.min(ICL_glasso)
  lambda_glasso[sim] = grid[minim_glasso]
  aaa = em.glasso(Y = Y, K = S, delta = init, gamma = gamma_HMM, mu=mu.s, Sigma=Sigma.s, rho = lambda_glasso[sim], err = 1e-4, iterMax = 5e2, traceEM = F)
  ARI_glasso[sim] = adjustedRandIndex(apply(aaa$post, 2, which.max), state)
  llk = aaa$loglik
  
  aaa$omega = simplify2array(aaa$Theta)
  t.time_glasso[sim] = aaa$timetot
  t.iter_glasso[sim] = aaa$iterations
  errorRate_glasso[sim] = classError(apply(aaa$post, 2, which.max), state)$errorRate
  state.order = order(apply(aaa$omega[2:4,1,], 2, which.max))
  # state.order = order(1-diag(aaa$gamma), decreasing = T)
  FP_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] != 0 & theta[,,1] == 0) / 2
  FP_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] != 0 & theta[,,2] == 0) / 2
  FP_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] != 0 & theta[,,3] == 0) / 2
  TFP_1_glasso[sim] = 1 - FP_1_glasso[sim] / (P*(P-1)/2 - P + 1)
  TFP_2_glasso[sim] = 1 - FP_2_glasso[sim] / (P*(P-1)/2 - P + 2)
  TFP_3_glasso[sim] = 1 - FP_3_glasso[sim] / (P*(P-1)/2 - P + 3)
  FN_1_glasso[sim] = sum(aaa$omega[,,state.order[1]] == 0 & theta[,,1] != 0) / 2
  FN_2_glasso[sim] = sum(aaa$omega[,,state.order[2]] == 0 & theta[,,2] != 0) / 2
  FN_3_glasso[sim] = sum(aaa$omega[,,state.order[3]] == 0 & theta[,,3] != 0) / 2
  TFN_1_glasso[sim] = 1 - FN_1_glasso[sim] / (P-1)
  TFN_2_glasso[sim] = 1 - FN_2_glasso[sim] / (P-2)
  TFN_3_glasso[sim] = 1 - FN_3_glasso[sim] / (P-3)
  M_1_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[1]]) - theta[,,1])^2)
  M_2_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[2]]) - theta[,,2])^2)
  M_3_glasso[sim] = sum((cov2cor(aaa$omega[,,state.order[3]]) - theta[,,3])^2)
  p_geom_1[sim] = 1-diag(aaa$gamma)[state.order[1]]
  p_geom_2[sim] = 1-diag(aaa$gamma)[state.order[2]]
  p_geom_3[sim] = 1-diag(aaa$gamma)[state.order[3]]
}

matrix = cbind(t.time_glasso,t.iter_glasso,lambda_glasso,ARI_glasso,errorRate_glasso,M_1_glasso,M_2_glasso,M_3_glasso,FP_1_glasso,FP_2_glasso,FP_3_glasso,TFP_1_glasso,TFP_2_glasso,TFP_3_glasso,FN_1_glasso,FN_2_glasso,FN_3_glasso,TFN_1_glasso,TFN_2_glasso,TFN_3_glasso,p_geom_1,p_geom_2,p_geom_3)
file_name = paste("N1000_dGeom_eOutliers_glasso_K3", ".csv", sep = "")
write.csv(matrix, file = file_name)