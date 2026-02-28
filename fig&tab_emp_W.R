rm(list = ls())
graphics.off()
gc()
library(cluster)
library(readr)
library(readxl)
library(MASS)
library(dplyr)
library(ggplot2)
library(car)
library(glasso)
library(mhsmm)
library(markovchain)
library(mclust)
library(doBy)
library(qgraph)
library(viridis)
library(igraph)
library(xts)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(ggcorrplot)
library(corpcor)
source("MainFunctions.R")

# Load the data here
# load("df_all_1725.RData")
# load("analisi_empirica_2025_M30_R30_boot_rev.RData")
# load("analisi_empirica_2025_M30_R20_boot_rev.RData")
# load("analisi_empirica_2025_M30_R50_boot_rev3.RData")
# load("analisi_empirica_2025_M30_R50_boot_20260224.RData")
# load("analisi_empirica_2025_M30_R50_boot_rev_Luca.RData")
# ! For GitHub:
load("df_woMSCI_1725.RData")
load("analisi_empirica_2025_M30_R50_boot_woMSCI.RData")
#########################################################################################################################

######## State relabeling based on the volatility ordering used in the first version of the paper #############
K <- dim(aaa$sigma)[3]

vol <- sapply(1:K, function(k) sum(diag(aaa$sigma[,, k])))
low_k <- which.min(vol)
high_k <- which.max(vol)
med_k <- setdiff(1:K, c(low_k, high_k))

perm <- c(low_k, high_k, med_k)
perm

relabel_hsmm <- function(aaa, perm) {
  out <- aaa

  if (!is.null(out$init)) {
    out$init <- out$init[perm]
  }
  if (!is.null(out$gamma)) {
    out$gamma <- out$gamma[perm, perm, drop = FALSE]
  }
  if (!is.null(out$d)) {
    out$d <- out$d[, perm, drop = FALSE]
  }
  if (!is.null(out$u)) {
    out$u <- out$u[, perm, drop = FALSE]
  }
  if (!is.null(out$emission_prob)) {
    out$emission_prob <- out$emission_prob[, perm, drop = FALSE]
  }
  if (!is.null(out$fden)) {
    out$fden <- out$fden[, perm, drop = FALSE]
  }
  if (!is.null(out$sigma)) {
    out$sigma <- out$sigma[,, perm, drop = FALSE]
  }
  if (!is.null(out$S_tilde)) {
    out$S_tilde <- out$S_tilde[,, perm, drop = FALSE]
  }
  if (!is.null(out$omega)) {
    out$omega <- out$omega[,, perm, drop = FALSE]
  }
  if (!is.null(out$omegaT)) {
    out$omegaT <- out$omegaT[,, perm, drop = FALSE]
  }
  if (!is.null(out$omega.sig)) {
    out$omega.sig <- out$omega.sig[,, perm, drop = FALSE]
  }

  out
}

aaa <- relabel_hsmm(aaa, perm)

d.sd <- d.sd[, perm, drop = FALSE]
gamma.sd <- gamma.sd[perm, perm, drop = FALSE]
omega.sd <- omega.sd[,, perm, drop = FALSE]

#numero archi totali
(sum(aaa$omega[,, 1] != 0) - P) / 2
(sum(aaa$omega[,, 2] != 0) - P) / 2
(sum(aaa$omega[,, 3] != 0) - P) / 2
#numero archi significativi
(sum(aaa$omega.sig[,, 1] != 0) - P) / 2
(sum(aaa$omega.sig[,, 2] != 0) - P) / 2
(sum(aaa$omega.sig[,, 3] != 0) - P) / 2

###### GRAFICI ######

labels = c(
  "BNB",
  "BTC",
  "DOGE",
  "ETH",
  "XRP",
  "SPX",
  "STOXX",
  "JPY",
  "GBP",
  "EUR",
  "CNY",
  "CHF",
  "CL",
  "HO",
  "NG",
  "RB",
  "BZ",
  "MAT",
  "REIT",
  "IND",
  "IT"
)
colrs = c("tomato", "dodgerblue2", "gold1", "darkorchid2", "green4")

g <- g_adj <- deg <- eig <- eig_adj <- e <- df.e <- df.e.b <- M.parcorr <- betw <- list()
Theta_mmdl <- list()
Sigma_mmdl <- list()

for (j in 1:S) {
  Theta_mmdl[[j]] <- aaa$omega.sig[,, j]
}

d = P
Adj_mmdl <- array(data = NA, dim = c(d, d, S))
for (j in 1:S) {
  Adj_mmdl[,, j] <- ifelse(
    as.matrix(Theta_mmdl[[j]]) != 0 &
      row(as.matrix(Theta_mmdl[[j]])) != col(as.matrix(Theta_mmdl[[j]])),
    1,
    0
  )
}

w_dist <- function(w) 1 / pmax(abs(w), .Machine$double.eps)
for (j in 1:S) {
  M.parcorr[[j]] <- wi2net(Theta_mmdl[[j]])
  g[[j]] <- graph_from_adjacency_matrix(
    as.matrix(M.parcorr[[j]]),
    mode = "upper",
    weighted = TRUE,
    diag = FALSE
  ) #grafo pesato
  g_adj[[j]] <- graph_from_adjacency_matrix(
    as.matrix(Adj_mmdl[,, j]),
    mode = "upper",
    weighted = F,
    diag = FALSE
  ) #grafo non pesato
  e[[j]] <- as_edgelist(g[[j]])
  # Count the number of degree for each node:
  deg[[j]] <- degree(g[[j]], mode = "all")
  w <- abs(E(g[[j]])$weight)
  eig[[j]] <- eigen_centrality(g[[j]], weights = w)$vector
  eig_adj[[j]] <- eigen_centrality(g_adj[[j]])$vector
  betw[[j]] <- edge_betweenness(g[[j]], weights = w_dist(w), directed = FALSE)
  df.e[[j]] <- as.data.frame(cbind(e[[j]], E(g[[j]])$weight))
  df.e.b[[j]] <- as.data.frame(cbind(
    e[[j]],
    10 * (1e-05 + betw[[j]]) / max(betw[[j]] + 1e-05)
  )) # add 1e-05 to avoid zeros
  V(g[[j]])[1:5]$color <- colrs[1] #crypto
  V(g[[j]])[6:7]$color <- colrs[2] #stock
  V(g[[j]])[8:12]$color <- colrs[3] #currency
  V(g[[j]])[13:17]$color <- colrs[4] #energy
  V(g[[j]])[18:21]$color <- colrs[5] # macroeconomy
  V(g[[j]])$label.font = 2
}

####### Centrality measures #####

centr_clo(g[[1]], mode = "all")
centr_clo(g[[2]], mode = "all")
centr_clo(g[[3]], mode = "all")

centr_betw(g[[1]], directed = FALSE)
centr_betw(g[[2]], directed = FALSE)
centr_betw(g[[3]], directed = FALSE)

#### Number of edges ####
sum(deg[[1]]) / 2
sum(deg[[2]]) / 2
sum(deg[[3]]) / 2

edge_density(g[[1]])
edge_density(g[[2]])
edge_density(g[[3]])
#
group_list <- list(
  "Cryptocurrencies" = 1:5,
  "Stock" = 6:7,
  "Currency" = 8:12,
  "Energy" = 13:17,
  "Macroeconomy" = 18:21
)

group_cols <- colrs
vertex_names <- labels

####### State 1: grafi con nodi proporzionali alla degree centrality ##
postscript("graph_K1.revB.eps", width = 980, height = 1010)
qgraph(
  df.e[[1]], ######si aggiunge e toglie la legenda ora solo inserendo legend=T
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Undirected graphs
  labels = labels,
  color = V(g[[1]])$color,
  vsize = deg[[1]] / max(deg[[1]]) * 12, # Regola la grandezza dei nodi
  label.cex = .9, # Grandezza dei label dentro i nodi.
  label.color = 'black', # string on label colors
  label.prop = 1.2,
  repulsion = .9,
  legend = F,
  legend.mode = "style2",
  legend.cex = .75,
  groups = group_list,
  nodeNames = vertex_names
)


### State 1: grafi con nodi proporzionali alla eigenvector centrality e archi prop alla edge-betweeness ##
# postscript("graph_K1.revB_eig.eps", width = 980, height = 1010)
# postscript("graph_K1.rev25.v2.eps", width = 980, height = 1010)
cairo_pdf(
  file = "graph_K1.rev25.v2.pdf",
  width = 980,
  height = 1010
)
qgraph(
  df.e.b[[1]], ######si aggiunge e toglie la legenda ora solo inserendo legend=T
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Directed graphs
  edge.color = "grey30", # Colore degli archi
  labels = labels, ######versione precedente !!!
  color = V(g[[1]])$color, ######versione precedente !!!
  vsize = 13 * (eig[[1]] / max(eig[[1]]))^(1 / 3), #farla non lineare, magari quadratica
  # label.cex = 1.3, # Grandezza dei label dentro i nodi. ###prima era 0.9
  label.cex = 1.1, #dati 25
  label.color = 'black', # string on label colors
  label.prop = 1.2, # proportion of the width of the node that the label scales
  repulsion = .9, #vecchio è 0.9
  # repulsion = 0.8, #dati25,
  edge.width = 3 * (betw[[1]] / max(betw[[1]]))^(1 / 3),
  legend = F,
  legend.mode = "style2",
  legend.cex = .75,
  groups = group_list,
  # color=group_cols, ##versione con legenda, da eliminare nella versione precedente
  nodeNames = vertex_names
)
dev.off()

####### State 2: grafi con nodi proporzionali alla degree centrality ##
postscript("graph_K2.revB.eps", width = 980, height = 1010)
qgraph(
  df.e[[2]],
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Directed graphs
  labels = labels,
  color = V(g[[2]])$color,
  vsize = deg[[2]] / max(deg[[2]]) * 10, # Regola la grandezza dei nodi
  label.cex = 1.5, # Grandezza dei label dentro i nodi. ###prima era 0.9
  label.color = 'black', # string on label colors
  label.prop = 0.9, # proportion of the width of the node that the label scales
  repulsion = .9,
  edge.width = abs(df.e[[2]][, 3]) * 2,
  legend = F,
  legend.mode = "style2",
  legend.cex = .45,
  groups = group_list
)
dev.off()

### State 2: grafi con nodi proporzionali alla eigenvector centrality e archi prop alla edge-betweeness ##
# postscript("graph_K2.revB_eig.eps", width = 980, height = 1010)
postscript("graph_K2.rev25.v2.eps", width = 980, height = 1010)
qgraph(
  df.e.b[[2]],
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Directed graphs
  edge.color = "grey30", # Colore degli archi
  labels = labels,
  color = V(g[[2]])$color,
  vsize = 13 * (eig[[2]] / max(eig[[2]]))^(1 / 3),
  # label.cex = 1.3, # Grandezza dei label dentro i nodi.
  label.cex = 1.1, #dati 25
  label.color = 'black', # string on label colors
  label.prop = 1.2, # proportion of the width of the node that the label scales
  repulsion = .9, #vecchio è 0.9
  # repulsion = .8, #vecchio è 0.9, dati25
  edge.width = 3 * (betw[[2]] / max(betw[[2]]))^(1 / 3),
  legend = F,
  legend.mode = "style2",
  legend.cex = .45,
  groups = group_list
)
dev.off()

#1010x980 se salvo da export

####### State 3: grafi con nodi proporzionali alla degree centrality ##
postscript("graph_K3.revB.eps", width = 980, height = 1010)
qgraph(
  df.e[[3]],
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Directed graphs
  labels = labels,
  color = V(g[[3]])$color,
  vsize = deg[[3]] / max(deg[[3]]) * 10, # Regola la grandezza dei nodi
  label.cex = 1.5, # Grandezza dei label dentro i nodi. ###prima era 0.9
  label.color = 'black', # string on label colors
  label.prop = 0.9, # proportion of the width of the node that the label scales
  repulsion = .9,
  edge.width = abs(df.e[[3]][, 3]) * 2,
  legend = F,
  legend.mode = "style2",
  legend.cex = .45,
  groups = group_list
)
dev.off()


### State 3: grafi con nodi proporzionali alla eigenvector centrality e archi prop alla edge-betweeness ##
# postscript("graph_K3.revB_eig.eps", width = 980, height = 1010)
postscript("graph_K3.rev25.v2.eps", width = 980, height = 1010)
qgraph(
  df.e.b[[3]], ######si aggiunge e toglie la legenda ora solo inserendo legend=T
  layout = "spring", # Senza questo le variabili sono disposte a cerchio
  fade = F, # Fa sì che i nodi non siano trasparenti
  directed = F, # Directed graphs
  edge.color = "gray30", # Colore degli archi
  labels = labels,
  color = V(g[[3]])$color,
  vsize = 13 * (eig[[3]] / max(eig[[3]]))^(1 / 3),
  # label.cex = 1.3, # Grandezza dei label dentro i nodi. ###prima era 0.9
  label.cex = 1.1, #dati 25
  label.color = 'black', # string on label colors
  label.prop = 1.2, # proportion of the width of the node that the label scales
  repulsion = 0.9, #vecchio è 0.9
  # repulsion = 0.8, #dati25, vecchio è 0.9
  # edge.width = 3*(betw[[3]]/max(betw[[3]]))^(1/3),
  edge.width = 3 * (betw[[3]] / max(betw[[3]]))^(1 / 3),
  legend = F,
  legend.mode = "style2",
  legend.cex = .75,
  groups = group_list
)
dev.off()

############################################################################################

#### Correlazione media per ogni stato:
mean(abs(as.vector(Theta_mmdl[[1]][upper.tri(Theta_mmdl[[1]], diag = F)])))
mean(abs(as.vector(Theta_mmdl[[2]][upper.tri(Theta_mmdl[[2]], diag = F)])))
mean(abs(as.vector(Theta_mmdl[[3]][upper.tri(Theta_mmdl[[3]], diag = F)])))
# mean(abs(as.vector(Theta_mmdl[[4]][upper.tri(Theta_mmdl[[4]], diag = F)])))

###### Predicted sojourn distributions ######
personal_colors <- c(
  "State 1" = "#F08080", # Rosso
  "State 2" = "#66CDAA", # Verde
  "State 3" = "#1874CD"
) # Blu


d_df <- data.frame(
  day = 1:M,
  value = c(aaa$d[, 1], aaa$d[, 2], aaa$d[, 3]),
  State = rep(c("State 1", "State 2", "State 3"), each = M)
)


# ksmoothed-nonparametric
d.new = dksmoothed(d = aaa$d)

d_df <- data.frame(
  day = 1:M,
  value = c(d.new),
  value_min = pmax(c(d.new - qnorm(0.975) * d.sd), 0),
  value_max = c(d.new + qnorm(0.975) * d.sd),
  State = rep(c("State 1", "State 2", "State 3"), each = M)
)


# ggsave("Pred_sojtimes25.eps", height = 3.18, width = 11.68, units = "in")
ggsave("Pred_sojtimes25.v2.eps", height = 3.18, width = 7.5, units = "in")
dd <- ggplot(
  d_df,
  aes(
    x = day,
    y = value,
    fill = factor(State, levels = c("State 1", "State 2", "State 3"))
  )
) +
  geom_bar(stat = "identity", width = 1) +
  geom_errorbar(
    aes(ymin = value_min, ymax = value_max),
    width = 0.5,
    position = position_dodge(.9)
  ) +
  facet_wrap(
    ~ factor(State, levels = c("State 1", "State 2", "State 3")),
    ncol = 3,
    nrow = 1,
    scales = "free_y"
  ) +
  xlab("Sojourn Time") +
  ylab("Probability") +
  ylim(0, .35) +
  # scale_fill_manual(values = c("State 1" = "red", "State 2" = "blue", "State 3" = "green")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    strip.text.x = element_text(size = 14),
    plot.margin = unit(c(0, 0.25, 0, 0), "inches") # top, right, bottom, left
  ) +
  scale_fill_manual(values = personal_colors)
print(dd)
dev.off()
#789x608 versione Luca dati 25
# 1121x305

###### Predicted posterior probabilities ######
posterior_df <- data.frame(
  day = as.Date(rownames(Y)),
  value = c(aaa$u[, 1], aaa$u[, 2], aaa$u[, 3]),
  State = rep(c("State 1", "State 2", "State 3"), each = N)
)

posterior_max <- data.frame(
  day = as.Date(rownames(Y)),
  value = apply(aaa$u, 1, which.max)
)

p <- ggplot(posterior_df, aes(x = day, y = value)) +
  geom_line() +
  facet_wrap(
    ~ factor(State, levels = c("State 1", "State 2", "State 3")),
    nrow = 3
  ) +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  aes(colour = State) +
  xlab("") +
  geom_smooth(method = "loess", color = "black") +
  xlab("Date") +
  ylab("Posterior Probabilities") +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 8),
    strip.text.x = element_text(
      size = 14
    ),
    plot.margin = unit(
      c(0, 0.25, 0, 0), #top,right,bottom,left
      "inches"
    )
  )
p
# 1121x305

####### Predicted sequence of hidden states over time #####
d_list <- list(aaa$d[, 1], aaa$d[, 2], aaa$d[, 3])
Gamma_star <- hsmm2hmm(aaa$gamma, d_list)
# emissions_star <- do.call(
#   cbind,
#   lapply(1:ncol(aaa$fden), function(i) {
#     matrix(rep(aaa$fden[, i], M), nrow = nrow(aaa$fden), ncol = M)
#   })
# )
emissions_star <- do.call(
  cbind,
  lapply(1:ncol(aaa$emission_prob), function(i) {
    matrix(
      rep(aaa$emission_prob[, i], M),
      nrow = nrow(aaa$emission_prob),
      ncol = M
    )
  })
)
initials_star <- c(
  aaa$init[1],
  rep(0, M - 1),
  rep(aaa$init[2], M),
  rep(aaa$init[3], M)
)


## with max a posteriori
max_post <- apply(aaa$u, 1, which.max)
## with Viterbi algorithm
seq_Viterbi <- Viterbi(
  Y[, 1],
  transProbs = Gamma_star,
  emissionProbs = emissions_star,
  initial_distribution = initials_star
)
Viterbi_transformed <- findInterval(seq_Viterbi, c(M + 1, 2 * M + 1)) + 1


sum(max_post != Viterbi_transformed)


pred.post_df <- data.frame(
  day = as.Date(rownames(Y)),
  value = c(max_post, Viterbi_transformed),
  Algorithm = rep(c("Max a Posteriori", "Viterbi"), each = N)
)


p1 <- ggplot(pred.post_df, aes(x = day, y = value)) +
  geom_line() +
  facet_wrap(
    ~ factor(Algorithm, levels = c("Max a Posteriori", "Viterbi")),
    nrow = 2
  ) +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  scale_y_continuous(breaks = 1:3) +
  aes(colour = Algorithm) +
  xlab("") +
  xlab("Date") +
  ylab("States") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16),
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.y = element_text(size = 12),
    strip.text.x = element_text(
      size = 14
    ),
    plot.margin = unit(
      c(0, 0.25, 0, 0), #bottom,top,left,right
      "inches"
    )
  )
p1

library(gridExtra)
grid.arrange(p, p1, nrow = 2)
#1121*610

table(Viterbi_transformed) / N
table(max_post) / N


###### All together returns plot #######
date = as.Date(rownames(Y))
ret <- Y

ret_crypto_df <- data.frame(
  date,
  ret[, 1],
  ret[, 2],
  ret[, 3],
  ret[, 4],
  ret[, 5],
  State = factor(Viterbi_transformed)
)

ret_stock_df <- data.frame(
  date,
  ret[, 6],
  ret[, 7],
  State = factor(Viterbi_transformed)
)

ret_currency_df <- data.frame(
  date,
  ret[, 8],
  ret[, 9],
  ret[, 10],
  ret[, 11],
  ret[, 12],
  State = factor(Viterbi_transformed)
)

ret_energy_df <- data.frame(
  date,
  ret[, 13],
  ret[, 14],
  ret[, 15],
  ret[, 16],
  ret[, 17],
  State = factor(Viterbi_transformed)
)

ret_macroeconomy_df <- data.frame(
  date,
  ret[, 18],
  ret[, 19],
  ret[, 20],
  ret[, 21],
  State = factor(Viterbi_transformed)
)


labels = c(
  "BNB",
  "BTC",
  "DOGE",
  "ETH",
  "XRP",
  "SP",
  "STOXX",
  "JPY",
  "GBP",
  "EUR",
  "CNY",
  "CHF",
  "CL",
  "HO",
  "NG",
  "RB",
  "BZ",
  "MAT",
  "REIT",
  "IND",
  "IT"
)


colnames(ret_crypto_df) <- c(
  "Date",
  c("BNB", "BTC", "DOGE", "ETH", "XRP"),
  "State"
)


colnames(ret_stock_df) <- c("Date", c("SP", "STOXX"), "State")


colnames(ret_currency_df) <- c(
  "Date",
  c("JPY", "GBP", "EUR", "CNY", "CHF"),
  "State"
)


colnames(ret_energy_df) <- c("Date", c("CL", "HO", "NG", "RB", "BZ"), "State")


colnames(ret_macroeconomy_df) <- c(
  "Date",
  c("MAT", "REIT", "IND", "IT"),
  "State"
)

############ Creation of time intervals for the horizontal bars ############
state_int <- ret_crypto_df %>%
  group_by(State) %>%
  mutate(start_date = Date, end_date = lead(Date, default = last(Date))) %>%
  ungroup()
############

gg1 <- melt(
  ret_crypto_df,
  id.vars = c("Date", "State"),
  variable.name = "Cryptocurrency",
  value.name = "Returns"
)
gg2 <- melt(
  ret_stock_df,
  id.vars = c("Date", "State"),
  variable.name = "Stock",
  value.name = "Returns"
)
gg3 <- melt(
  ret_currency_df,
  id.vars = c("Date", "State"),
  variable.name = "Currency",
  value.name = "Returns"
)
gg4 <- melt(
  ret_energy_df,
  id.vars = c("Date", "State"),
  variable.name = "Energy",
  value.name = "Returns"
)
gg5 <- melt(
  ret_macroeconomy_df,
  id.vars = c("Date", "State"),
  variable.name = "Macroeconomy",
  value.name = "Returns"
)

personal_colors <- c(
  "1" = "#F08080", # Rosso
  "2" = "#66CDAA", # Verde
  "3" = "#1874CD"
) # Blu


p1 <- ggplot(gg1, aes(x = Date, y = Returns, color = Cryptocurrency)) +
  geom_rect(
    data = state_int,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(gg1$Returns) - 0.2 * (max(gg1$Returns) - min(gg1$Returns)),
      ymax = min(gg1$Returns) - 0.05 * (max(gg1$Returns) - min(gg1$Returns)),
      fill = State
    ),
    inherit.aes = F
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'),
    legend.key.height = unit(0.8, 'cm'),
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = "%Y", breaks = "1 year") +
  xlab("") +
  scale_fill_manual(values = personal_colors)

p1

p2 <- ggplot(gg2, aes(x = Date, y = Returns, color = Stock)) +
  geom_rect(
    data = state_int,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(gg2$Returns) - 0.2 * (max(gg2$Returns) - min(gg2$Returns)),
      ymax = min(gg2$Returns) - 0.05 * (max(gg2$Returns) - min(gg2$Returns)),
      fill = State
    ),
    inherit.aes = F
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("") +
  scale_fill_manual(values = personal_colors)

p2

p3 <- ggplot(gg3, aes(x = Date, y = Returns, color = Currency)) +
  geom_rect(
    data = state_int,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(gg3$Returns) - 0.2 * (max(gg3$Returns) - min(gg3$Returns)),
      ymax = min(gg3$Returns) - 0.05 * (max(gg3$Returns) - min(gg3$Returns)),
      fill = State
    ),
    inherit.aes = F
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("") +
  scale_fill_manual(values = personal_colors)

p3

p4 <- ggplot(gg4, aes(x = Date, y = Returns, color = Energy)) +
  geom_rect(
    data = state_int,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(gg4$Returns) - 0.2 * (max(gg4$Returns) - min(gg4$Returns)),
      ymax = min(gg4$Returns) - 0.05 * (max(gg4$Returns) - min(gg4$Returns)),
      fill = State
    ),
    inherit.aes = F
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("") +
  scale_fill_manual(values = personal_colors)

p4

p5 <- ggplot(gg5, aes(x = Date, y = Returns, color = Macroeconomy)) +
  geom_rect(
    data = state_int,
    aes(
      xmin = start_date,
      xmax = end_date,
      ymin = min(gg5$Returns) - 0.2 * (max(gg5$Returns) - min(gg5$Returns)),
      ymax = min(gg5$Returns) - 0.05 * (max(gg5$Returns) - min(gg5$Returns)),
      fill = State
    ),
    inherit.aes = F
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("") +
  scale_fill_manual(values = personal_colors)
p5

ggarrange(p1, p2, p3, p4, p5, nrow = 5, ncol = 1, legend = "top")
# 2000x1500

###### All together price plot #######

normalTS <- function(x) (x - min(x)) / (max(x) - min(x))

library(ggplot2)
library(reshape2)
library(gridExtra)
# load("df_prices1726.RData")
load("df_prices1725.RData")
date <- as.Date(rownames(prices.df))
N <- dim(prices.df)[1]
prices <- as.matrix(prices.df)
prices_norm <- apply(prices, 2, normalTS)

price_crypto_df <- data.frame(
  date,
  prices_norm[, 1],
  prices_norm[, 2],
  prices_norm[, 3],
  prices_norm[, 4],
  prices_norm[, 5]
)

price_stock_df <- data.frame(date, prices_norm[, 6], prices_norm[, 7])

price_currency_df <- data.frame(
  date,
  prices_norm[, 8],
  prices_norm[, 9],
  prices_norm[, 10],
  prices_norm[, 11],
  prices_norm[, 12]
)

price_energy_df <- data.frame(
  date,
  prices_norm[, 13],
  prices_norm[, 14],
  prices_norm[, 15],
  prices_norm[, 16],
  prices_norm[, 17]
)

price_macroeconomy_df <- data.frame(
  date,
  prices_norm[, 18],
  prices_norm[, 19],
  prices_norm[, 20],
  prices_norm[, 21]
)

labels = c(
  "BNB",
  "BTC",
  "DOGE",
  "ETH",
  "XRP",
  "SP",
  "STOXX",
  "JPY",
  "GBP",
  "EUR",
  "CNY",
  "CHF",
  "CL",
  "HO",
  "NG",
  "RB",
  "BZ",
  "MAT",
  "REIT",
  "IND",
  "IT"
)
colnames(price_crypto_df) <- c("Date", c("BNB", "BTC", "DOGE", "ETH", "XRP"))
colnames(price_stock_df) <- c("Date", c("SP", "STOXX"))
colnames(price_currency_df) <- c("Date", c("JPY", "GBP", "EUR", "CNY", "CHF"))
colnames(price_energy_df) <- c("Date", c("CL", "HO", "NG", "RB", "BZ"))
colnames(price_macroeconomy_df) <- c("Date", c("MAT", "REIT", "IND", "IT"))

gg1 <- melt(
  price_crypto_df,
  id = c("Date"),
  variable.name = "Cryptocurrency",
  value.name = "Price"
)
gg2 <- melt(
  price_stock_df,
  id = c("Date"),
  variable.name = "Stock",
  value.name = "Price"
)
gg3 <- melt(
  price_currency_df,
  id = c("Date"),
  variable.name = "Currency",
  value.name = "Price"
)
gg4 <- melt(
  price_energy_df,
  id = c("Date"),
  variable.name = "Energy",
  value.name = "Price"
)
gg5 <- melt(
  price_macroeconomy_df,
  id = c("Date"),
  variable.name = "Macroeconomy",
  value.name = "Price"
)

p1 <- ggplot(gg1, aes(x = Date, y = Price, color = Cryptocurrency)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("")

p2 <- ggplot(gg2, aes(x = Date, y = Price, color = Stock)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("")


p3 <- ggplot(gg3, aes(x = Date, y = Price, color = Currency)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("")

p4 <- ggplot(gg4, aes(x = Date, y = Price, color = Energy)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("")

p5 <- ggplot(gg5, aes(x = Date, y = Price, color = Macroeconomy)) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22),
    axis.title.x = element_text(size = 28),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0) +
  geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") +
  xlab("Date")

library(Cairo)
library(ggpubr)
pp <- ggarrange(p1, p2, p3, p4, p5, nrow = 5, ncol = 1, legend = "top")
CairoPDF(
  file = "pricePlot.pdf",
  width = 28,
  height = 21
)

print(pp)
dev.off()


###### Heatmap ICL ######
mat <- t(HSMM.crit)
df <- melt(mat)
colnames(df) <- c("Lambda", "State", "Value")

# Find the minimum value
min_value <- min(df$Value)
df_min <- df %>% filter(Value == min_value)

# Do the heatmap
ggplot(df, aes(x = State, y = Lambda, fill = Value)) +
  labs(y = expression(alpha), x = "K") +
  geom_tile() +
  scale_fill_gradientn(colors = viridis::viridis(20)) +
  geom_segment(
    data = df_min,
    aes(x = State - .5, xend = State + .5, y = Lambda, yend = Lambda),
    color = "red",
    size = 1
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(
    fill = guide_colorbar(
      barheight = unit(4.5, "inches"),
      barwidth = unit(0.3, "inches")
    )
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.title.x = element_text(size = 28, face = "italic"),
    axis.title.y = element_text(size = 28),
    axis.text.x = element_text(size = 20, angle = 0),
    axis.text.y = element_text(size = 20),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  )
#760x540

######## Partial correlation plot #######
X <- as.matrix(Y)
P <- corpcor::cor2pcor(corpcor::cor.shrink(X))
diag(P) <- 1
ord <- hclust(dist(1 - abs(P)))$order
P <- P[ord, ord]
colnames(P) <- labels[ord]
rownames(P) <- labels[ord]

# Plot
p <- ggcorrplot(
  P,
  type = "full",
  lab = FALSE,
  show.diag = TRUE,
  colors = c("red", "white", "blue"),
  outline.col = "white"
) +
  scale_y_discrete(limits = rev(colnames(P))) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(
    fill = guide_colorbar(
      barheight = unit(4.5, "inches"),
      barwidth = unit(0.3, "inches")
    )
  ) +
  theme(
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18, angle = 90, hjust = 1, vjust = 1),
    legend.key.size = unit(1.8, 'cm'), #change legend key size
    legend.key.height = unit(0.8, 'cm'), #change legend key height
    legend.key.width = unit(2.5, 'cm')
  )

p
#557x660

#643x760

###################################################################################################################
