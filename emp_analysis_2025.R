# rm(list = ls()) #! ok per github
# graphics.off()
# gc()

# ! levare per github
lib = getwd()
repos = "http://cran.uk.r-project.org"
.libPaths(c(.libPaths(), lib))

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

#####################################################################
# Load the data here
load("df_all_1725.RData")
# ! For GitHub:
# load("df_woMSCI_1725.RData")

#####################################################################
# Summary statistics
means = rep(0, ncol(Y))
std.dev = rep(0, ncol(Y))
skews = rep(0, ncol(Y))
kurts = rep(0, ncol(Y))
mins = rep(0, ncol(Y))
medians = rep(0, ncol(Y))
maxs = rep(0, ncol(Y))

for (i in 1:ncol(Y)) {
  means[i] = mean(Y[, i])
  std.dev[i] = sd(Y[, i])
  skews[i] = skewness(Y[, i])
  kurts[i] = kurtosis(Y[, i])
  mins[i] = min(Y[, i])
  medians[i] = median(Y[, i])
  maxs[i] = max(Y[, i])
}

summary(means)
summary(std.dev)
summary(skews)
summary(kurts)
summary(mins)
summary(medians)
summary(maxs)

###########################################################################################
# Fit the model
# S_grid = c(2, 3, 4)
S_grid = c(2, 3, 4, 5, 6)
N = nrow(Y)
P = ncol(Y)
M = 30
grid = seq(3, 13, length.out = 200)
ICL = matrix(NA, length(grid), length(S_grid))
R = 50 # number of restart
tmp = lapply(1:R, function(x) {
  list()
})
tmp.llk = array(NA, c(length(S_grid), length(grid), R))
output.HSMM = lapply(1:length(S_grid), function(x) {
  lapply(1:length(grid), function(x) {
    list()
  })
})

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
      A = diag(hmm_init$estimate@transitionMatrix)
      A = 1 - A
      if (S == 2) {
        gamma = matrix(c(0, 1, 1, 0), S, S, byrow = TRUE)
      } else if (S > 2) {
        gamma = hmm_init$estimate@transitionMatrix
        diag(gamma) = 0
        gamma = gamma / rowSums(gamma)
      } else {
        gamma = 1
        A = 1
      }

      Par = list(
        #mu = matrix(0, P, S),
        gamma = gamma,
        init = init,
        d = sapply(1:S, function(i) dgeom(1:M - 1, A[i])),
        sigma = replicate(S, diag(P))
      )

      tmp[[r]] = try(
        EM_HSMM(
          Y = Y,
          S = S,
          Par = Par,
          M = M,
          sojourn.distribution = "poisson",
          lambda = lambda,
          pen_EBIC = 0.5,
          seed = 1,
          err = 1e-8,
          iterMax = 2e3
        ),
        silent = TRUE
      )

      if (inherits(tmp[[r]], "try-error")) {
        tmp.llk[s, i, r] = -Inf # tmp[[r]] = NULL
      } else {
        tmp.llk[s, i, r] = as.numeric(tmp[[r]]$loglik)
      }
    }
    output.HSMM[[s]][[i]] = tmp[[which.max(tmp.llk[s, i, ])]]
  }
}

llk.mat = array(NA, c(length(S_grid), length(grid), R))
for (s in 1:length(S_grid)) {
  for (i in 1:length(grid)) {
    llk.mat[s, i, ] = abs(
      (max(tmp.llk[s, i, ]) - tmp.llk[s, i, ]) / max(tmp.llk[s, i, ])
    ) <
      10^-3
  }
}

rowMeans(llk.mat[1, , ])
rowMeans(llk.mat[2, , ])
apply(llk.mat, 1, mean)
mean(llk.mat)
# (max(tmp.llk) - tmp.llk)/max(tmp.llk) < 10^-3

save.image(paste0("analisi_empirica_2025_M", M, "_R", R, ".RData"))
# ! For GitHub
# save.image(paste0("analisi_empirica_2025_M", M, "_R", R, "_woMSCI.RData"))
