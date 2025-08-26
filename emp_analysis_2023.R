S_grid = c(2,3,4,5,6)
N = nrow(Y)
P = ncol(Y)
M = 30
grid = seq(3,13,length.out=200)
ICL = matrix(NA, length(grid), length(S_grid))
R = 20 # number of restart
tmp = lapply(1:R, function(x) {list()})
tmp.llk = matrix(NA, R, 1)
output.HSMM = lapply(1:length(S_grid), function(x) {lapply(1:length(grid), function(x) {list()})})

for (s in 1:length(S_grid)) {
  S = S_grid[s]
  for (i in 1:length(grid)) {
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
      if(S == 2) {
        gamma = matrix(c(0,1,1,0),S,S,byrow=TRUE)
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
        d = sapply(1:S, function(i) dgeom(1:M-1, A[i])),
        sigma = replicate(S, diag(P))
      )
      
      tmp[[r]] = try(EM_HSMM(Y = Y, S = S, Par = Par, M = M, sojourn.distribution = "poisson", lambda = lambda, pen_EBIC = 0.5, seed = 1, err = 1e-4, iterMax = 5e2))
      tmp.llk[r] = as.numeric(tmp[[r]]$loglik)
    }
    output.HSMM[[s]][[i]] = tmp[[which.max(tmp.llk)]]
  }
}

# save the environment
save.image("analisi_empirica_2023_M30_R20_rev.RData")

# determine the optimal number of states and Lasso regularization parameter
HSMM.crit = matrix(NA, length(S_grid), length(grid))
for(s in 1:length(S_grid)) {
  for(i in 1:length(grid)) {
    if(!is.null(output.HSMM[[s]][[i]]$ICL)) HSMM.crit[s,i] = output.HSMM[[s]][[i]]$ICL
  }
}
colnames(HSMM.crit) = grid
rownames(HSMM.crit) = S_grid
HSMM.crit
which(HSMM.crit == min(HSMM.crit, na.rm = T), arr.ind = T)

##########################################################
##########################################################
##########################################################

aaa = output.HSMM[[2]][[92]]