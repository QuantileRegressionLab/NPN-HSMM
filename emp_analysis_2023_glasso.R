rm(list = ls())
graphics.off()
gc()
library(cluster)
library(readr)
library(readxl)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(markovchain)
library(mclust)
library(doBy)
library(qgraph)
library(viridis)
library(moments)
library(igraph)
library(xts)
library(PearsonDS)
source("MainFunctions.R")

###########################################################################################
# Fit the model
S_grid = c(2,3,4,5,6)
N = nrow(Y)
P = ncol(Y)
M = 30
grid = seq(1,13,length.out=200)
ICL = matrix(NA, length(grid), length(S_grid))
R = 20 # number of restart
tmp = lapply(1:R, function(x) {list()})
tmp.llk = matrix(NA, R, 1)
output.HMM = lapply(1:length(S_grid), function(x) {lapply(1:length(grid), function(x) {list()})})

for (s in 1:length(S_grid)) {
  print(s)
  S = S_grid[s]
  for (i in 1:length(grid)) {
    print(i)
    lambda = grid[i]
    for (r in 1:R) {
      set.seed(r)
      states_init = pam(x = Y, k = S)
      states_init$clustering = sample(S, N, replace = T)
      # fit the underlying Markov chain
      hmm_init = markovchainFit(states_init$cluster)
      init = rep(0, S)
      init[states_init$cluster[1]] = 1
      #init[1] = 1
      gamma = hmm_init$estimate@transitionMatrix
      
      Par = list(
        gamma = gamma,
        init = init,
        sigma = replicate(S, diag(P))
      )
      
      tmp[[r]] = try(EM_HMM(Y = Y, S = S, Par = Par, lambda = lambda, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
      tmp.llk[r] = as.numeric(tmp[[r]]$loglik)
    }
    output.HMM[[s]][[i]] = tmp[[which.max(tmp.llk)]]
  }
}


# determine the optimal number of states and Lasso regularization parameter
HMM.crit = matrix(NA, length(S_grid), length(grid))
for(s in 1:length(S_grid)) {
  for(i in 1:length(grid)) {
    if(!is.null(output.HMM[[s]][[i]]$ICL)) HMM.crit[s,i] = output.HMM[[s]][[i]]$ICL
  }
}
colnames(HMM.crit) = grid
rownames(HMM.crit) = S_grid
HMM.crit
which(HMM.crit == min(HMM.crit, na.rm = T), arr.ind = T)

##########################################################
##########################################################
##########################################################
source("em_glasso.R")

# Fit the model
S_grid = c(2,3,4,5,6)
N = nrow(Y)
P = ncol(Y)
M = 30
grid = seq(1,13,length.out=200)
ICL = matrix(NA, length(grid), length(S_grid))
R = 20 # number of restart
tmp = lapply(1:R, function(x) {list()})
tmp.llk = matrix(NA, R, 1)
output.glasso = lapply(1:length(S_grid), function(x) {lapply(1:length(grid), function(x) {list()})})

for (s in 1:length(S_grid)) {
  print(s)
  S = S_grid[s]
  for (i in 1:length(grid)) {
    print(i)
    lambda = grid[i]
    for (r in 1:R) {
      set.seed(r)
      states_init = pam(x = Y, k = S)
      states_init$clustering = sample(S, N, replace = T)
      # fit the underlying Markov chain
      hmm_init = markovchainFit(states_init$cluster)
      init = rep(0, S)
      init[states_init$cluster[1]] = 1
      #init[1] = 1
      gamma = hmm_init$estimate@transitionMatrix
      
      Sigma.s = replicate(S, diag(P), simplify = F)
      mu.s = matrix(0, S, P)
      
      tmp[[r]] = try(em.glasso(Y = Y, K = S, delta = init, gamma = gamma, mu=mu.s, Sigma=Sigma.s, rho = lambda, err = 1e-4, iterMax = 5e2, traceEM = F))
      tmp.llk[r] = as.numeric(tmp[[r]]$loglik)
    }
    output.glasso[[s]][[i]] = tmp[[which.max(tmp.llk)]]
  }
}


# determine the optimal number of states and Lasso regularization parameter
glasso.crit = matrix(NA, length(S_grid), length(grid))
for(s in 1:length(S_grid)) {
  for(i in 1:length(grid)) {
    if(!is.null(output.glasso[[s]][[i]]$pen.criteria$ICL)) glasso.crit[s,i] = output.glasso[[s]][[i]]$pen.criteria$ICL
  }
}
colnames(glasso.crit) = grid
rownames(glasso.crit) = S_grid
glasso.crit
which(glasso.crit == min(glasso.crit, na.rm = T), arr.ind = T)

save.image("glasso_results.RData")
