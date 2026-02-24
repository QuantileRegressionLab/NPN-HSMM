rm(list = ls())
library(markovchain)
library(cluster)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(mclust)
library(PearsonDS)
library(foreach)
library(parallel)
library(doParallel)
source("MainFunctions.R")


# make clusters
cl = detectCores()
registerDoParallel(cl)


# set the true model parameters
N <- c(300, 500, 1500)
S <- 2 # number of states
P <- 30 # number of response
K <- 1 # the intercept
M <- 30 # maximum time spent in each state

mu = list(matrix(0, P, K), matrix(0, P, K))

theta = array(0, dim = c(P, P, S))
covarianza = array(0, dim = c(P, P, S))
for (j in 1:P) {
  for (k in 1:P) {
    if (j == k) {
      theta[j, k, 1] = 1
      theta[j, k, 2] = 1
    }
    if (abs(j - k) == 1) {
      theta[j, k, 1] = 0.4
    }
    if (abs(j - k) == 2) {
      theta[j, k, 2] = 0.4
    }
  }
}
covarianza[,, 1] = solve(theta[,, 1])
covarianza[,, 2] = solve(theta[,, 2])

sigma_sim <- list(covarianza[,, 1], covarianza[,, 2])
init_sim <- c(1, 0) # initial distribution
gamma_sim <- matrix(c(0, 1, 1, 0), S, S, byrow = TRUE) # transition probability matrix

dl <- 1
d <- lapply(1:dl, function(x) {
  list()
}) # sojourn distributions
# d[[1]] <- list(c(0.2,0.3,0.5), c(0.1,0.2,0.7))
d[[1]] <- list(shift.poi(1:M, 10, 1), shift.poi(1:M, 5, 1))
d[[2]] <- list(shift.nb(1:M, 8, 0.5), shift.nb(1:M, 4, 0.6))
d[[3]] <- list(dgeom(1:M - 1, 0.2), dgeom(1:M - 1, 0.1))
d_true <- d

sojourn = c("poisson", "nbinom", "geometric")

error <- c("n", "outliers")


# grids = list(seq(7,17,length.out=100), seq(10,19,length.out=100), seq(16,25,length.out=100))
grids = list(
  seq(7, 17, length.out = 100),
  seq(10, 19, length.out = 100),
  seq(16, 46, length.out = 100)
)

nsim = 300

#################### N1000 + POISSON + MULTIVARIATE GAUSSIAN ##################################

print("N1000 + POISSON + MULTIVARIATE GAUSSIAN")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 1 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
results_list = list()
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
    d = d_true[[b]],
    error = error[c]
  )
  Y = data_gen$series
  state = data_gen$state

  ################################################################################
  # Initialization of model parameters (HMM and HSMM)
  #
  # - Initial state sequence (states_init$clustering) is obtained by
  #   randomly assigning observations to the $K$ latent states according to a
  #   Multinomial distribution with equal probabilities 1/K.
  # - The off-diagonal elements of the transition matrix (gamma) are computed as
  #   proportions of transition from the generated partition.
  # - Emission covariance matrices (sigma) are initialized as diagonal matrices.
  # - Sojourn-time distributions (d) are are estimated from the initial
  #   partition assuming a Geometric distribution as in HMMs.

  # - This initialization is the same in every setting
  ################################################################################

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
  gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = replicate(S, dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = foreach(
    i = 1:length(grid),
    .combine = c,
    .packages = c("mvtnorm", "mclust", "mhsmm", "car", "glasso")
  ) %dopar%
    {
      aaa = try(
        EM_HSMM(
          Y = Y,
          S = S,
          Par = Par_HSMM,
          M = M,
          sojourn.distribution = sojourn[b],
          lambda = grid[i],
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-4,
          iterMax = 5e2
        ),
        silent = TRUE
      )
      if (inherits(aaa, "try-error")) {
        aaa$ICL = Inf
      } else {
        aaa$ICL
      }
      aaa$ICL
    }
  minim_HSMM = which.min(ICL_HSMM)
  lambda_HSMM = grid[minim_HSMM]
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM,
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = replicate(S, dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM,
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }

  aaa$lambda_HSMM = lambda_HSMM
  results_list[[sim]] = aaa
}

save.image("N1500_dPoisson_eMVNorm.RData")

########################################
# bootstrap resamples
R = 100 # number of resamples
results_list.boot <- list()
for (sim in 1:nsim) {
  print(sim)

  # order the estimates
  # state.order = order(apply(results_list[[sim]]$omega[2:4,1,], 2, which.max))
  # state.order = order(results_list[[sim]]$sojourn$lambda_poi, decreasing = T)
  state.order = 1:S
  delta = results_list[[sim]]$init[state.order]
  gamma = results_list[[sim]]$gamma[state.order, state.order]
  sigma = results_list[[sim]]$S_tilde[,, state.order]
  rho = results_list[[sim]]$omegaT[,, state.order]
  rho = sapply(
    1:S,
    function(s) as.matrix(forceSymmetric(cov2cor(rho[,, s]))),
    simplify = F
  )
  sigma = sapply(1:S, function(s) cov2cor(solve(rho[[s]])), simplify = F)
  # sigma = sapply(1:S, function(s) cov2cor(sigma[,,s]), simplify = F)
  d = results_list[[sim]]$d[, state.order]
  d = as.list(as.data.frame(d))
  lambda_HSMM = results_list[[sim]]$lambda_HSMM

  results_list.tmp = foreach(
    r = 1:R,
    .packages = c("mvtnorm", "mclust", "mhsmm", "car", "glasso")
  ) %dopar%
    {
      print(r)
      set.seed(r)

      data_boot = boot.hsmm.multi.gen(
        ns = N[a],
        P = P,
        K = K,
        m = S,
        delta = delta,
        gamma = gamma,
        rho = sigma,
        d = d
      )
      Y_boot = data_boot$series

      Par_HSMM = list(
        gamma = gamma,
        init = delta,
        d = results_list[[sim]]$d,
        sigma = results_list[[sim]]$sigma
      )
      fit.boot = try(
        EM_HSMM(
          Y = Y_boot,
          S = S,
          Par = Par_HSMM,
          M = M,
          sojourn.distribution = sojourn[b],
          lambda = lambda_HSMM,
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-4,
          iterMax = 5e2
        ),
        silent = TRUE
      )
      if (inherits(fit.boot, "try-error")) {
        fit.boot = NULL
      }
      fit.boot
    }
  results_list.boot[[sim]] = results_list.tmp
}

####################
save.image("N1500_dPoisson_eMVNorm_boot.RData")
####################
coverage.summary(
  fit = results_list,
  boot = results_list.boot,
  theta = theta,
  d_true = d_true,
  gamma_sim = gamma_sim,
  sojourn.distribution = "poisson",
  alpha = 0.05
)


####################
####################
####################

#################### N1000 + NEG BIN + MULTIVARIATE GAUSSIAN #########################################

print("N1000 + NEG BIN + MULTIVARIATE GAUSSIAN")
a = 3 # a=1: N=300; a=2: N=500; a=3: N=1000
b = 2 # b=1: Shifted Poisson; b=2: Neg Bin; b=3: Geometrica
c = 1 # c=1: Multivariata Gaussiana; c=2: Multivariata Gaussiana con outliers
grid = grids[[a]]
results_list = list()
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
    d = d_true[[b]],
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
  gamma_HSMM = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)

  Par_HSMM = list(
    gamma = gamma_HSMM,
    init = init,
    d = replicate(S, dgeom(1:M - 1, .1)),
    sigma = replicate(S, diag(P))
  )

  # FITTING DELLA SERIE STORICA SIMULATA

  ICL_HSMM = foreach(
    i = 1:length(grid),
    .combine = c,
    .packages = c("mvtnorm", "mclust", "mhsmm", "car", "glasso")
  ) %dopar%
    {
      aaa = try(
        EM_HSMM(
          Y = Y,
          S = S,
          Par = Par_HSMM,
          M = M,
          sojourn.distribution = sojourn[b],
          lambda = grid[i],
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-4,
          iterMax = 5e2
        ),
        silent = TRUE
      )
      if (inherits(aaa, "try-error")) {
        aaa$ICL = Inf
      } else {
        aaa$ICL
      }
      aaa$ICL
    }
  minim_HSMM = which.min(ICL_HSMM)
  lambda_HSMM = grid[minim_HSMM]
  aaa = EM_HSMM(
    Y = Y,
    S = S,
    Par = Par_HSMM,
    M = M,
    sojourn.distribution = sojourn[b],
    lambda = lambda_HSMM,
    pen_EBIC = 0.5,
    seed = 1,
    err = 1e-4,
    iterMax = 5e2
  )
  ARI_HSMM = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
  llk = aaa$loglik
  iter = 0
  while (ARI_HSMM < 0.3 & iter < 10) {
    iter = iter + 1
    states_init$clustering = sample(S, N[a], replace = T)
    hmm_init = markovchainFit(states_init$cluster)
    init = rep(0, S)
    init[states_init$cluster[1]] = 1
    Par_HSMM = list(
      gamma = gamma_HSMM,
      init = init,
      d = replicate(S, dgeom(1:M - 1, .1)),
      sigma = replicate(S, diag(P))
    )
    ccc = EM_HSMM(
      Y = Y,
      S = S,
      Par = Par_HSMM,
      M = M,
      sojourn.distribution = sojourn[b],
      lambda = lambda_HSMM,
      pen_EBIC = 0.5,
      seed = 1,
      err = 1e-4,
      iterMax = 5e2
    )
    if (ccc$loglik > llk) {
      aaa = ccc
      llk = ccc$loglik
      ARI_HSMM = adjustedRandIndex(apply(aaa$u, 1, which.max), state)
    }
  }

  aaa$lambda_HSMM = lambda_HSMM
  results_list[[sim]] = aaa
}

save.image("N1500_dNegBin_eMVNorm.RData")

########################################
# bootstrap resamples
R = 100 # number of resamples
results_list.boot <- list()
for (sim in 1:nsim) {
  print(sim)

  # order the estimates
  # state.order = order(apply(results_list[[sim]]$omega[2:4,1,], 2, which.max))
  # state.order = order(results_list[[sim]]$sojourn$lambda_poi, decreasing = T)
  state.order = 1:S
  delta = results_list[[sim]]$init[state.order]
  gamma = results_list[[sim]]$gamma[state.order, state.order]
  sigma = results_list[[sim]]$S_tilde[,, state.order]
  rho = results_list[[sim]]$omegaT[,, state.order]
  rho = sapply(
    1:S,
    function(s) as.matrix(forceSymmetric(cov2cor(rho[,, s]))),
    simplify = F
  )
  sigma = sapply(1:S, function(s) cov2cor(solve(rho[[s]])), simplify = F)
  # sigma = sapply(1:S, function(s) cov2cor(sigma[,,s]), simplify = F)
  d = results_list[[sim]]$d[, state.order]
  d = as.list(as.data.frame(d))
  lambda_HSMM = results_list[[sim]]$lambda_HSMM

  results_list.tmp = foreach(
    r = 1:R,
    .packages = c("mvtnorm", "mclust", "mhsmm", "car", "glasso")
  ) %dopar%
    {
      print(r)
      set.seed(r)

      data_boot = boot.hsmm.multi.gen(
        ns = N[a],
        P = P,
        K = K,
        m = S,
        delta = delta,
        gamma = gamma,
        rho = sigma,
        d = d
      )
      Y_boot = data_boot$series

      Par_HSMM = list(
        gamma = gamma,
        init = delta,
        d = results_list[[sim]]$d,
        sigma = results_list[[sim]]$sigma
      )
      fit.boot = try(
        EM_HSMM(
          Y = Y_boot,
          S = S,
          Par = Par_HSMM,
          M = M,
          sojourn.distribution = sojourn[b],
          lambda = lambda_HSMM,
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-4,
          iterMax = 5e2
        ),
        silent = TRUE
      )
      if (inherits(fit.boot, "try-error")) {
        fit.boot = NULL
      }
      fit.boot
    }
  results_list.boot[[sim]] = results_list.tmp
}


####################
save.image("N1500_dNegBin_eMVNorm_boot.RData")
####################
coverage.summary(
  fit = results_list,
  boot = results_list.boot,
  theta = theta,
  d_true = d_true,
  gamma_sim = gamma_sim,
  sojourn.distribution = "nbinom",
  alpha = 0.05
)
