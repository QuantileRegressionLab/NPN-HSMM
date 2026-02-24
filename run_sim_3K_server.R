# rm(list = ls())
# gc()
# graphics.off()

lib = getwd()
repos = "http://cran.uk.r-project.org"
.libPaths(c(.libPaths(), lib))


library(markovchain)
library(cluster)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(mclust)
library(PearsonDS)
source("MainFunctions.R")


# set the true model parameters
N <- c(300, 500, 1000)
S <- 3 # number of states
P <- 30 # number of response
K <- 1 # the intercept
M <- 30 # maximum time spent in each state

mu = list(matrix(0, P, K), matrix(0, P, K), matrix(0, P, K))

theta = array(0, dim = c(P, P, S))
covarianza = array(0, dim = c(P, P, S))
for (j in 1:P) {
  for (k in 1:P) {
    if (j == k) {
      theta[j, k, 1] = 1
      theta[j, k, 2] = 1
      theta[j, k, 3] = 1
    }
    if (abs(j - k) == 1) {
      theta[j, k, 1] = 0.4
    }
    if (abs(j - k) == 2) {
      theta[j, k, 2] = 0.4
    }
    if (abs(j - k) == 3) {
      theta[j, k, 3] = 0.4
    }
  }
}
covarianza[,, 1] = solve(theta[,, 1])
covarianza[,, 2] = solve(theta[,, 2])
covarianza[,, 3] = solve(theta[,, 3])

sigma_sim <- list(covarianza[,, 1], covarianza[,, 2], covarianza[,, 3])
init_sim <- c(1, 0, 0) # initial distribution
gamma_sim <- matrix(
  c(0, 0.8, 0.2, 0.4, 0, 0.6, 0.7, 0.3, 0),
  S,
  S,
  byrow = TRUE
) # transition probability matrix
gamma_sim <- matrix(0.5, S, S)
diag(gamma_sim) <- 0

dl <- 1
d <- lapply(1:dl, function(x) {
  list()
}) # sojourn distributions
# d[[1]] <- list(c(0.2,0.3,0.5), c(0.1,0.2,0.7))
d[[1]] <- list(
  shift.poi(1:M, 15, 1),
  shift.poi(1:M, 10, 1),
  shift.poi(1:M, 5, 1)
)
d[[2]] <- list(
  shift.nb(1:M, 10, 0.5),
  shift.nb(1:M, 7, 0.6),
  shift.nb(1:M, 4, 0.7)
)
d[[3]] <- list(dgeom(1:M - 1, 0.3), dgeom(1:M - 1, 0.15), dgeom(1:M - 1, 0.1))

sojourn = c("poisson", "nbinom", "geometric")

error <- c("n", "outliers")


grids = list(
  seq(7, 17, length.out = 100),
  seq(10, 19, length.out = 100),
  seq(16, 25, length.out = 100)
)

nsim = 300

#################### N300 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N300 + POISSON + MULTIVARIATE GAUSSIAN")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = replicate(S, dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N300_dPoisson_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N500 + POISSON + MULTIVARIATE GAUSSIAN")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N500_dPoisson_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + POISSON + MULTIVARIATE GAUSSIAN")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N1000_dPoisson_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + POISSON + OUTLIERS #########################################

print("N300 + POISSON + OUTLIERS")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N300_dPoisson_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + POISSON + OUTLIERS #########################################

print("N500 + POISSON + OUTLIERS")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N500_dPoisson_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + POISSON + OUTLIERS #########################################

print("N1000 + POISSON + OUTLIERS")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_poi_1 = rep(0, nsim)
lambda_poi_2 = rep(0, nsim)
lambda_poi_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$lambda_poi, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_poi_1[sim] = aaa$sojourn$lambda_poi[state.order[1]]
  lambda_poi_2[sim] = aaa$sojourn$lambda_poi[state.order[2]]
  lambda_poi_3[sim] = aaa$sojourn$lambda_poi[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_poi_1,
  lambda_poi_2,
  lambda_poi_3
)
file_name = paste("N1000_dPoisson_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N300 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N300_dNegBin_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N500 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N500_dNegBin_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N1000 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N1000_dNegBin_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + NEG BIN + OUTLIERS #########################################

print("N300 + NEG BIN + OUTLIERS")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N300_dNegBin_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + NEG BIN + OUTLIERS #########################################

print("N500 + NEG BIN + OUTLIERS")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N500_dNegBin_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + NEG BIN + OUTLIERS #########################################

print("N1000 + NEG BIN + OUTLIERS")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
lambda_NB_1 = rep(0, nsim)
lambda_NB_2 = rep(0, nsim)
lambda_NB_3 = rep(0, nsim)
p_NB_1 = rep(0, nsim)
p_NB_2 = rep(0, nsim)
p_NB_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$sojourn$mu, decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  lambda_NB_1[sim] = aaa$sojourn$mu[state.order[1]]
  lambda_NB_2[sim] = aaa$sojourn$mu[state.order[2]]
  lambda_NB_3[sim] = aaa$sojourn$mu[state.order[3]]
  p_NB_1[sim] = aaa$sojourn$prob[state.order[1]]
  p_NB_2[sim] = aaa$sojourn$prob[state.order[2]]
  p_NB_3[sim] = aaa$sojourn$prob[state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  lambda_NB_1,
  lambda_NB_2,
  lambda_NB_3,
  p_NB_1,
  p_NB_2,
  p_NB_3
)
file_name = paste("N1000_dNegBin_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N300 + GEOM + MULTIVARIATE GAUSSIAN")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N300_dGeom_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N500 + GEOM + MULTIVARIATE GAUSSIAN")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N500_dGeom_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + GEOM + MULTIVARIATE GAUSSIAN")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N1000_dGeom_eMVNorm_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N300 + GEOM + OUTLIERS ####################################################

print("N300 + GEOM + OUTLIERS")
a = 1 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N300_dGeom_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N500 + GEOM + OUTLIERS ####################################################

print("N500 + GEOM + OUTLIERS")
a = 2 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N500_dGeom_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)

#################### N1000 + GEOM + OUTLIERS ####################################################

print("N1000 + GEOM + OUTLIERS")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 3 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 2 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
t.time_HSMM = rep(0, nsim)
t.time_HMM = rep(0, nsim)
t.iter_HSMM = rep(0, nsim)
t.iter_HMM = rep(0, nsim)
lambda_HSMM = rep(0, nsim)
lambda_HMM = rep(0, nsim)
ARI_HSMM = rep(0, nsim)
ARI_HMM = rep(0, nsim)
errorRate_HSMM = rep(0, nsim)
errorRate_HMM = rep(0, nsim)
M_1_HSMM = rep(0, nsim)
M_2_HSMM = rep(0, nsim)
M_3_HSMM = rep(0, nsim)
M_1_HMM = rep(0, nsim)
M_2_HMM = rep(0, nsim)
M_3_HMM = rep(0, nsim)
FP_1_HSMM = rep(0, nsim)
FP_2_HSMM = rep(0, nsim)
FP_3_HSMM = rep(0, nsim)
FP_1_HMM = rep(0, nsim)
FP_2_HMM = rep(0, nsim)
FP_3_HMM = rep(0, nsim)
TFP_1_HSMM = rep(0, nsim)
TFP_2_HSMM = rep(0, nsim)
TFP_3_HSMM = rep(0, nsim)
TFP_1_HMM = rep(0, nsim)
TFP_2_HMM = rep(0, nsim)
TFP_3_HMM = rep(0, nsim)
FN_1_HSMM = rep(0, nsim)
FN_2_HSMM = rep(0, nsim)
FN_3_HSMM = rep(0, nsim)
FN_1_HMM = rep(0, nsim)
FN_2_HMM = rep(0, nsim)
FN_3_HMM = rep(0, nsim)
TFN_1_HSMM = rep(0, nsim)
TFN_2_HSMM = rep(0, nsim)
TFN_3_HSMM = rep(0, nsim)
TFN_1_HMM = rep(0, nsim)
TFN_2_HMM = rep(0, nsim)
TFN_3_HMM = rep(0, nsim)
p_geom_1 = rep(0, nsim)
p_geom_2 = rep(0, nsim)
p_geom_3 = rep(0, nsim)
for (sim in 1:nsim) {
  print(sim)
  set.seed(sim)

  data_gen <- hsmm.multi.gen(
    ns = N[a],
    P = P,
    K = K,
    m = S,
    delta = init_sim,
    gamma = gamma_sim,
    mu = mu,
    rho = sigma_sim,
    d = d[[b]],
    error = error[c]
  )
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
  if (S == 2) {
    gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
  } else {
    gamma_HSMM = gamma_HMM
    diag(gamma_HSMM) = 0
    gamma_HSMM = gamma_HSMM / rowSums(gamma_HSMM)
  }

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  Par_HMM = list(
    gamma = gamma_HMM,
    init = init,
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = rep(0, 100)
  ICL_HMM = rep(0, 100)
  j = 0
  for (i in grid) {
    j = j + 1
    aaa = try(
      EM_HSMM(
        Y = Y,
        S = S,
        Par = Par_HSMM,
        M = M,
        sojourn.distribution = sojourn[b],
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
    if (inherits(aaa, "try-error")) {
      ICL_HSMM[j] = Inf
    } else {
      ICL_HSMM[j] = aaa$ICL
    }
    bbb = try(
      EM_HMM(
        Y = Y,
        S = S,
        Par = Par_HMM,
        lambda = i,
        pen_EBIC = 0.5,
        seed = 1,
        err = 1e-4,
        iterMax = 5e2
      ),
      silent = TRUE
    )
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
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  bbb = EM_HMM(
    Y = Y,
    S = S,
    Par = Par_HMM,
    lambda = lambda_HMM[sim],
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  ARI_HMM[sim] = adjustedRandIndex(apply(bbb$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM[sim] < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = cbind(dgeom(1:M - 1, .1), dgeom(1:M - 1, .1), dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM[sim] = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }
  llk = bbb$loglik
  iter = 0
  while (ARI_HMM[sim] < 0.3 & iter < 10) {
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
    ccc = EM_HMM(
      Y = Y,
      S = S,
      Par = Par_HMM,
      lambda = lambda_HMM[sim],
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
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
  state.order = order(aaa$d[1, ], decreasing = T)
  FP_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  FP_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] != 0 & theta[,, 1] == 0) / 2
  FP_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] != 0 & theta[,, 2] == 0) / 2
  FP_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] != 0 & theta[,, 3] == 0) / 2
  TFP_1_HSMM[sim] = 1 - FP_1_HSMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HSMM[sim] = 1 - FP_2_HSMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HSMM[sim] = 1 - FP_3_HSMM[sim] / (P * (P - 1) / 2 - P + 3)
  TFP_1_HMM[sim] = 1 - FP_1_HMM[sim] / (P * (P - 1) / 2 - P + 1)
  TFP_2_HMM[sim] = 1 - FP_2_HMM[sim] / (P * (P - 1) / 2 - P + 2)
  TFP_3_HMM[sim] = 1 - FP_3_HMM[sim] / (P * (P - 1) / 2 - P + 3)
  FN_1_HSMM[sim] = sum(aaa$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HSMM[sim] = sum(aaa$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HSMM[sim] = sum(aaa$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  FN_1_HMM[sim] = sum(bbb$omega[,, state.order[1]] == 0 & theta[,, 1] != 0) / 2
  FN_2_HMM[sim] = sum(bbb$omega[,, state.order[2]] == 0 & theta[,, 2] != 0) / 2
  FN_3_HMM[sim] = sum(bbb$omega[,, state.order[3]] == 0 & theta[,, 3] != 0) / 2
  TFN_1_HSMM[sim] = 1 - FN_1_HSMM[sim] / (P - 1)
  TFN_2_HSMM[sim] = 1 - FN_2_HSMM[sim] / (P - 2)
  TFN_3_HSMM[sim] = 1 - FN_3_HSMM[sim] / (P - 3)
  TFN_1_HMM[sim] = 1 - FN_1_HMM[sim] / (P - 1)
  TFN_2_HMM[sim] = 1 - FN_2_HMM[sim] / (P - 2)
  TFN_3_HMM[sim] = 1 - FN_3_HMM[sim] / (P - 3)
  M_1_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HSMM[sim] = sum((cov2cor(aaa$omega[,, state.order[3]]) - theta[,, 3])^2)
  M_1_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[1]]) - theta[,, 1])^2)
  M_2_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[2]]) - theta[,, 2])^2)
  M_3_HMM[sim] = sum((cov2cor(bbb$omega[,, state.order[3]]) - theta[,, 3])^2)
  p_geom_1[sim] = aaa$d[1, state.order[1]]
  p_geom_2[sim] = aaa$d[1, state.order[2]]
  p_geom_3[sim] = aaa$d[1, state.order[3]]
}

matrix = cbind(
  t.time_HSMM,
  t.time_HMM,
  t.iter_HSMM,
  t.iter_HMM,
  lambda_HSMM,
  lambda_HMM,
  ARI_HSMM,
  ARI_HMM,
  errorRate_HSMM,
  errorRate_HMM,
  M_1_HSMM,
  M_1_HMM,
  M_2_HSMM,
  M_2_HMM,
  M_3_HSMM,
  M_3_HMM,
  FP_1_HSMM,
  FP_1_HMM,
  FP_2_HSMM,
  FP_2_HMM,
  FP_3_HSMM,
  FP_3_HMM,
  TFP_1_HSMM,
  TFP_1_HMM,
  TFP_2_HSMM,
  TFP_2_HMM,
  TFP_3_HSMM,
  TFP_3_HMM,
  FN_1_HSMM,
  FN_1_HMM,
  FN_2_HSMM,
  FN_2_HMM,
  FN_3_HSMM,
  FN_3_HMM,
  TFN_1_HSMM,
  TFN_1_HMM,
  TFN_2_HSMM,
  TFN_2_HMM,
  TFN_3_HSMM,
  TFN_3_HMM,
  p_geom_1,
  p_geom_2,
  p_geom_3
)
file_name = paste("N1000_dGeom_eOutliers_K", S, ".csv", sep = "")
write.csv(matrix, file = file_name)
