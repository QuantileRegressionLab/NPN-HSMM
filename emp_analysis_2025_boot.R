rm(list = ls())
gc()
library(markovchain)
library(cluster)
library(MASS)
library(car)
library(glasso)
library(mhsmm)
library(mclust)
library(foreach)
library(parallel)
library(doParallel)

# load("analisi_empirica_2025_M30_R20_rev.RData")
# load("analisi_empirica_2025_M30_R30_rev.RData")
load("analisi_empirica_2025_M30_R50_rev.RData")
# load("analisi_empirica_2025_M30_R50_v2.RData")
# ! For GitHub
# load("analisi_empirica_2025_M30_R30_rev_woMSCI.RData")
source("MainFunctions.R")


HSMM.crit = matrix(NA, length(S_grid), length(grid))
for (s in 1:length(S_grid)) {
  for (i in 1:length(grid)) {
    if (!is.null(output.HSMM[[s]][[i]]$ICL)) {
      HSMM.crit[s, i] = output.HSMM[[s]][[i]]$ICL
    }
  }
}
colnames(HSMM.crit) = grid
rownames(HSMM.crit) = S_grid
HSMM.crit
idx.opt = which(HSMM.crit == min(HSMM.crit, na.rm = T), arr.ind = T)
idx.opt
##########################################################
##########################################################
##########################################################
S = S_grid[idx.opt[1]]
lambda_HSMM = grid[idx.opt[2]]

aaa = output.HSMM[[idx.opt[1]]][[idx.opt[2]]]


# make clusters
ncores = detectCores()
cl = makeCluster(ncores)
registerDoParallel(cl)

########################################
# bootstrap resamples
R = 500 # number of resamples

delta = aaa$init
gamma = aaa$gamma
rho = aaa$omegaT
rho = sapply(
  1:S,
  function(s) as.matrix(forceSymmetric(cov2cor(rho[,, s]))),
  simplify = F
)
sigma = sapply(
  1:S,
  function(s) cov2cor(as.matrix(nearPD(solve(rho[[s]]))$mat)),
  simplify = F
)
# sigma = sapply(1:S, function(s) cov2cor(sigma[,,s]), simplify = F)
d = aaa$d
# d = dksmoothed(d = d)
d = as.list(as.data.frame(d))

results_list.boot = foreach(
  r = 1:R,
  .packages = c("mvtnorm", "mclust", "mhsmm", "car", "glasso")
) %dopar%
  {
    tries = 0

    while (tries <= 10) {
      tries = tries + 1
      # set.seed(tries * r)

      data_boot = boot.hsmm.multi.gen(
        ns = N,
        P = P,
        K = 1,
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
        d = do.call(cbind, d),
        sigma = simplify2array(sigma)
      )
      oo.se = tryCatch(
        EM_HSMM(
          Y = Y_boot,
          S = S,
          Par = Par_HSMM,
          M = M,
          sojourn.distribution = "np",
          lambda = lambda_HSMM,
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-8,
          iterMax = 2e3
        ),
        error = function(e) {
          e
        }
      )
      if (
        !is(oo.se, "error") &
          !any(oo.se$d == 1) &
          !any(table(oo.se$predicted_state) <= 10)
      ) {
        break
      }
    }
    if (inherits(oo.se, "error")) {
      oo.se = NULL
    }

    fit.boot = oo.se
    # if(inherits(fit.boot, "try-error")) fit.boot = NULL
    # return(fit.boot)
  }


# estimate the standard errors based on the bootstrap resamples
gamma.array = array(NA, c(S, S, R))
d.array = array(NA, c(M, S, R))
omega.array = array(NA, c(P, P, S, R))
for (r in 1:R) {
  if (!is.null(results_list.boot[[r]])) {
    if (!any(results_list.boot[[r]]$d == 1)) {
      # non convergence
      # state.order = order(
      #   table(results_list.boot[[r]]$predicted_state),
      #   decreasing = T
      # )
      state.order = 1:S

      # state.order = order(
      #   table(aaa$predicted_state)#,
      #   # decreasing = T
      # )
      gamma.array[,, r] = results_list.boot[[r]]$gamma[state.order, state.order]
      omega.array[,,, r] = sapply(
        1:S,
        function(s) cov2cor(results_list.boot[[r]]$omegaT[,, state.order[s]]),
        simplify = "array"
      )
      d.array[,, r] = (results_list.boot[[r]]$d[, state.order])
      d.array[,, r] = dksmoothed(results_list.boot[[r]]$d[, state.order])
    }
  }
}

gamma.sd = apply(gamma.array, 1:2, sd.trim, trim = 0.05, na.rm = T)
gamma.sd

d.sd = apply(d.array, 1:2, sd.trim, trim = 0.05, na.rm = T)

omega.sd = apply(omega.array, 1:3, sd.trim, trim = 0.05, na.rm = T)


(sum(aaa$omega[,, 1] != 0) - P) / 2
(sum(aaa$omega[,, 2] != 0) - P) / 2
(sum(aaa$omega[,, 3] != 0) - P) / 2

alpha = 0.05
qqnorm = qnorm(1 - alpha / 2)
aaa$omega.sig = array(NA, c(P, P, S))
for (s in 1:S) {
  est = cov2cor(aaa$omegaT[,, s])

  # Confidence interval bounds
  lower = est - qqnorm * omega.sd[,, s]
  upper = est + qqnorm * omega.sd[,, s]

  # Indicator if 0 is inside the interval
  contains_zero = (lower <= 0 & upper >= 0)

  aaa$omega.sig[,, s] = (contains_zero == F) * cov2cor(aaa$omegaT[,, s]) # est
}

(sum(aaa$omega.sig[,, 1] != 0) - P) / 2
(sum(aaa$omega.sig[,, 2] != 0) - P) / 2
(sum(aaa$omega.sig[,, 3] != 0) - P) / 2

# boot_times = as.numeric(sapply(results_list.boot, function(x) x$time))
boot_times <- vapply(
  results_list.boot,
  function(x) if (length(x$time) == 0L) NA_real_ else as.numeric(x$time),
  FUN.VALUE = numeric(1)
)
mean(boot_times, na.rm = TRUE)
sd(boot_times, na.rm = TRUE)
sum(boot_times, na.rm = TRUE)

# iters = sapply(results_list.boot, function(x) x$iter)
iters <- vapply(
  results_list.boot,
  function(x) if (length(x$iter) == 0L) NA_real_ else as.numeric(x$iter),
  FUN.VALUE = numeric(1)
)
mean(iters, na.rm = TRUE)
sd(iters, na.rm = TRUE)
sum(iters, na.rm = TRUE)
sum(iters <= 500, na.rm = T)

####################
# save.image("analisi_empirica_2025_M30_R20_boot_rev.RData")
# save.image("analisi_empirica_2025_M30_R30_boot_rev.RData")
ts_tag <- format(Sys.time(), "%Y%m%d")
fname <- sprintf("analisi_empirica_2025_M30_R50_boot_%s.RData", ts_tag)
save.image(fname)
# ! For GitHub
# save.image("analisi_empirica_2025_M30_R50_boot_rev_woMSCI.RData")
####################
