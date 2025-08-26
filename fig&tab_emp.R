rm(list = ls())
graphics.off()
gc()
library(cluster)
library(readr)
library(readxl)
library(MASS)
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
source("MainFunctions.R")

#####################################################################
load("Y_all.RData")
###################################################################################################################
# GRAFICI

labels = c("BNB", "BTC", "DOGE", "ETH", "XRP", "SPX", "STOXX", "JPY", "GBP", "EUR", "CNY", "CHF", "CL", "HO", "NG", "RB", "BZ", "MAT", "REIT", "IND", "IT")
colrs = c("tomato", "dodgerblue2", "gold1", "darkorchid2", "green4")


g <- deg <- eig <- e <- df.e <- M.parcorr <- list()
Theta_mmdl <- list()
Sigma_mmdl <- list()

for (j in 1:S) {
  Theta_mmdl[[j]] <- aaa$omega[,,j]
}

d = P
Adj_mmdl <- array(data = NA, dim = c(d, d, S))
for (j in 1:S) {
  Adj_mmdl[,,j] <- ifelse(as.matrix(Theta_mmdl[[j]])!=0 & row(as.matrix(Theta_mmdl[[j]]))!=col(as.matrix(Theta_mmdl[[j]])),1,0)
}

for (j in 1:S) {
  M.parcorr[[j]] <- wi2net(round(Theta_mmdl[[j]],2))
  # g[[j]] <- graph.adjacency(as.matrix(M.parcorr[[j]]), mode="upper", weighted=TRUE, diag=FALSE)
  g[[j]] <- graph_from_adjacency_matrix(as.matrix(M.parcorr[[j]]), mode="upper", weighted=TRUE, diag=FALSE)
  e[[j]] <- get.edgelist(g[[j]])
  df.e[[j]] <- as.data.frame(cbind(e[[j]],E(g[[j]])$weight))
  # Count the number of degree for each node:
  deg[[j]] <- degree(g[[j]], mode="all")
  # eig[[j]] <- eigen_centrality(g[[j]])
  V(g[[j]])[1:5]$color<-colrs[1] #crypto
  V(g[[j]])[6:7]$color<-colrs[2] #stock
  V(g[[j]])[8:12]$color<-colrs[3] #currency
  V(g[[j]])[13:17]$color<-colrs[4] #energy
  V(g[[j]])[18:21]$color<-colrs[5] # macroeconomy
  V(g[[j]])$label.font=2
}


group_list <- list("Cryptocurrencies"=1:5,
                   "Stock"=6:7,
                   "Currency"=8:12,
                   "Energy"=13:17,
                   "Macroeconomy" = 18:21)

group_cols <- colrs
vertex_names <- labels

setwd("/Users/bia/PostDoc/Ricerca/HM_NonPar_Emilio/figures")
postscript("graph_K1.eps", width = 980, height = 1010)
qgraph(df.e[[1]], 
       layout = "spring", # Senza questo le variabili sono disposte a cerchio
       fade=F,  # Fa sì che i nodi non siano trasparenti
       directed = F, # Directed graphs
       # edge.color = "red", # Colore degli archi
       # cut=2, # tutti gli archi hanno lo stesso spessore
       labels = labels, 
       color = V(g[[1]])$color,
       # vsize = 9.5,  # Regola la grandezza dei nodi
       vsize = deg[[1]]/max(deg[[1]])*10,  # Regola la grandezza dei nodi
       # pie = deg[[1]]/max(deg[[1]]), # Mette ai bordi dei nodi una grafico a torta
       # pieBorder = 0.25,
       # pieColor = "lightblue",
       label.cex = .9,  # Grandezza dei label dentro i nodi
       label.color = 'black', # string on label colors
       label.prop = 0.9, # proportion of the width of the node that the label scales
       repulsion = .9,
       edge.width = abs(df.e[[1]][,3])*2,
       legend = F,
       legend.mode="style2",
       legend.cex=.20,
       groups=group_list
)
dev.off()

postscript("graph_K2.eps", width = 980, height = 1010)
qgraph(df.e[[2]], 
       layout = "spring", # Senza questo le variabili sono disposte a cerchio
       fade=F,  # Fa sì che i nodi non siano trasparenti
       directed = F, # Directed graphs
       # edge.color = "red", # Colore degli archi
       # cut=2, # tutti gli archi hanno lo stesso spessore
       labels = labels, 
       color = V(g[[2]])$color,
       # vsize = 9.5,  # Regola la grandezza dei nodi
       vsize = deg[[2]]/max(deg[[2]])*10,  # Regola la grandezza dei nodi
       # pie = deg[[1]]/max(deg[[1]]), # Mette ai bordi dei nodi una grafico a torta
       # pieBorder = 0.25,
       # pieColor = "lightblue",
       label.cex = .9,  # Grandezza dei label dentro i nodi
       label.color = 'black', # string on label colors
       label.prop = 0.9, # proportion of the width of the node that the label scales
       repulsion = .9,
       edge.width = abs(df.e[[2]][,3])*2,
       legend = F,
       legend.mode="style2",
       legend.cex=.20,
       groups=group_list
)
dev.off()

postscript("graph_K3.eps", width = 980, height = 1010)
qgraph(df.e[[3]], 
       layout = "spring", # Senza questo le variabili sono disposte a cerchio
       fade=F,  # Fa sì che i nodi non siano trasparenti
       directed = F, # Directed graphs
       # edge.color = "red", # Colore degli archi
       # cut=2, # tutti gli archi hanno lo stesso spessore
       labels = labels, 
       color = V(g[[3]])$color,
       # vsize = 9.5,  # Regola la grandezza dei nodi
       vsize = deg[[3]]/max(deg[[3]])*10,  # Regola la grandezza dei nodi
       # pie = deg[[1]]/max(deg[[1]]), # Mette ai bordi dei nodi una grafico a torta
       # pieBorder = 0.25,
       # pieColor = "lightblue",
       label.cex = .9,  # Grandezza dei label dentro i nodi
       label.color = 'black', # string on label colors
       label.prop = 0.9, # proportion of the width of the node that the label scales
       repulsion = .9,
       edge.width = abs(df.e[[3]][,3])*2,
       legend = F,
       legend.mode="style2",
       legend.cex=.20,
       groups=group_list
)
dev.off()

############################################################################################

###### Predicted sojourn distributions ######
personal_colors<- c("State 1" = "#F08080",  # Rosso
                    "State 2" = "#66CDAA",# Verde
                    "State 3" = "#1874CD") # Blu


d_df <- data.frame(
  day = 1:M,
  value = c(aaa$d[,1], aaa$d[,2], aaa$d[,3]),
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

ggsave("Pred_sojtimes.eps", height = 3.18, width = 11.68, units = "in")
dd <- ggplot(d_df, aes(x = day, y = value, fill = factor(State, levels = c("State 1", "State 2", "State 3")))) +
  geom_bar(stat = "identity", width = 1) +
  geom_errorbar(aes(ymin= value_min,
                    ymax= value_max),
                width=0.5,
                position=position_dodge(.9)) +
  facet_wrap(~factor(State, levels = c("State 1", "State 2", "State 3")), ncol = 3, nrow = 1, scales = "free_y") +
  xlab("Sojourn Time") +
  ylab("Probability") +
  ylim(0, 0.3) +
  # scale_fill_manual(values = c("State 1" = "red", "State 2" = "blue", "State 3" = "green")) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 21),
    axis.title.y = element_text(size = 21),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    strip.text.x = element_text(size = 14),
    plot.margin = unit(c(0, 0.25, 0, 0), "inches") # top, right, bottom, left
  ) + scale_fill_manual(values = personal_colors)
print(dd)
dev.off()

# 1121x305

###### Predicted posterior probabilities ######
posterior_df <- data.frame(
  day = as.Date(rownames(Y)),
  value = c(aaa$u[,1], aaa$u[,2], aaa$u[,3]),
  State = rep(c("State 1", "State 2", "State 3"), each = N)
)

posterior_max <- data.frame(
  day = as.Date(rownames(Y)),
  value = apply(aaa$u, 1, which.max)
)

p <- ggplot(posterior_df, aes(x=day, y=value)) +
  geom_line() + facet_wrap(~factor(State, levels = c("State 1", "State 2", "State 3")), nrow = 3) +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + 
  aes(colour = State)  +
  xlab("") + geom_smooth(method="loess", color="black") + xlab("Date") + 
  ylab("Posterior Probabilities") +
  theme(legend.position = "none", axis.title.x = element_text(size=21),
        axis.title.y = element_text(size=21), axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=8),
        strip.text.x = element_text(
          size = 14),
        plot.margin = unit(c(0, 0.25, 0, 0), #top,right,bottom,left
                           "inches")) 
p 
# 1121x305

####### Predicted sequence of hidden states over time #####
d_list <- list(aaa$d[,1], aaa$d[,2], aaa$d[,3])
Gamma_star <- hsmm2hmm(aaa$gamma, d_list)
emissions_star <- do.call(cbind, lapply(1:ncol(aaa$u), function(i) matrix(rep(aaa$u[, i], M), nrow = nrow(aaa$u), ncol = M)))
initials_star <- c(aaa$init[1], rep(0, M-1), rep(aaa$init[2], M), rep(aaa$init[3], M))


## with max a posteriori
max_post <- apply(aaa$u, 1, which.max)
## with Viterbi algorithm
seq_Viterbi <- Viterbi(Y[,1], transProbs = Gamma_star, emissionProbs = emissions_star,
                       initial_distribution = initials_star)
Viterbi_transformed <- findInterval(seq_Viterbi, c(M+1, 2*M+1)) + 1

sum(max_post != Viterbi_transformed)

pred.post_df <- data.frame(
  day = as.Date(rownames(Y)),
  value = c(max_post, Viterbi_transformed),
  Algorithm = rep(c("Max a Posteriori", "Viterbi"), each = N)
  # value = c(max_post),
  # Algorithm = rep(c("Max a Posteriori"), each = N)
)

p1 <- ggplot(pred.post_df, aes(x=day, y=value)) +
  geom_line() + facet_wrap(~factor(Algorithm, levels = c("Max a Posteriori", "Viterbi")), nrow = 2) +
  # geom_line() + facet_wrap(~factor(Algorithm, levels = c("Max a Posteriori")), nrow = 1) +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + 
  scale_y_continuous(breaks = 1:3) +
  aes(colour = Algorithm)  +
  xlab("") + xlab("Date") + 
  ylab("States") +
  theme(legend.position = "none",
        axis.text.x  = element_text(size = 16),
        # plot.margin = margin(, .3, , .2, "cm"),
        axis.title.x = element_text(size=21),
        axis.title.y = element_text(size=21),
        axis.text.y = element_text(size=12),
        strip.text.x = element_text(
          size = 14),
        plot.margin = unit(c(0, 0.25, 0, 0), #bottom,top,left,right
                           "inches")) 
p1


grid.arrange(p,p1, nrow = 2)

#1121*610


table(Viterbi_transformed)/N
table(max_post)/N

###### Alltogheter returns plot #######

# normalTS <- function(x) (x - min(x))/(max(x) - min(x))
date = as.Date(rownames(Y))
# date = c[2:length(c)]

ret <- Y
# prices_norm <- apply(prices, 2, normalTS)

# ret_crypto_df <- data.frame(date, 
#                                ret[,1], ret[,2], ret[,3], ret[,4], ret[,5], State = factor(max_post))
ret_crypto_df <- data.frame(date, 
                            ret[,1], ret[,2], ret[,3], ret[,4], ret[,5], State = factor(Viterbi_transformed))

# ret_stock_df <- data.frame(date,
#                            ret[,6], ret[,7], State = factor(max_post))
ret_stock_df <- data.frame(date,
                           ret[,6], ret[,7], State = factor(Viterbi_transformed))

# ret_currency_df <- data.frame(date,
#                               ret[,8], ret[,9], ret[,10], ret[,11], ret[,12], State = factor(max_post))
ret_currency_df <- data.frame(date,
                              ret[,8], ret[,9], ret[,10], ret[,11], ret[,12], State = factor(Viterbi_transformed))

# ret_energy_df <- data.frame(date,
#                             ret[,13], ret[,14], ret[,15], ret[,16], ret[,17], State = factor(max_post))
ret_energy_df <- data.frame(date,
                            ret[,13], ret[,14], ret[,15], ret[,16], ret[,17], State = factor(Viterbi_transformed))

# ret_macroeconomy_df <- data.frame(date,
#                                   ret[,18], ret[,19], ret[,20], ret[,21], State = factor(max_post))
ret_macroeconomy_df <- data.frame(date,
                                  ret[,18], ret[,19], ret[,20], ret[,21], State = factor(Viterbi_transformed))

labels = c("BNB", "BTC", "DOGE", "ETH", "XRP", "SP", "STOXX", "JPY", "GBP", "EUR", "CNY", "CHF", "CL", "HO", "NG", "RB", "BZ", "MAT", "REIT", "IND", "IT")

colnames(ret_crypto_df) <- c("Date", c("BNB", "BTC", "DOGE", "ETH", "XRP"), "State")

colnames(ret_stock_df) <- c("Date", c("SP", "STOXX"), "State")

colnames(ret_currency_df) <- c("Date", c("JPY","GBP","EUR","CNY", "CHF"), "State")

colnames(ret_energy_df) <- c("Date", c("CL", "HO", "NG", "RB", "BZ"), "State")

colnames(ret_macroeconomy_df) <- c("Date", c("MAT", "REIT", "IND", "IT"), "State")


############ Creazione degli intervalli di tempo per le barre orizzontali ############
state_int <- ret_crypto_df %>%
  group_by(State) %>%
  mutate(start_date = Date, end_date = lead(Date, default = last(Date))) %>%
  ungroup()
############

gg1 <- melt(ret_crypto_df, id.vars = c("Date", "State"), variable.name="Cryptocurrency", value.name="Returns")
gg2 <- melt(ret_stock_df, id.vars = c("Date", "State"), variable.name="Stock", value.name="Returns")
gg3 <- melt(ret_currency_df, id.vars = c("Date", "State"), variable.name="Currency", value.name="Returns")
gg4 <- melt(ret_energy_df, id.vars = c("Date", "State"), variable.name = "Energy", value.name = "Returns")
gg5 <- melt(ret_macroeconomy_df, id.vars = c("Date", "State"), variable.name="Macroeconomy", value.name="Returns")

# Definisci una palette di colori personalizzati
# personal_colors<- c("1" = "#E41A1C",  # Rosso
#                 "2" = "#4DAF4A",# Verde
#                 "3" = "#4876FF") # Blu  

personal_colors<- c("1" = "#F08080",  # Rosso
                    "2" = "#66CDAA",# Verde
                    "3" = "#1874CD") # Blu  

p1 <- ggplot(gg1, aes(x = Date, y = Returns, color = Cryptocurrency)) + 
  geom_rect(data = state_int, 
    aes(xmin = start_date, xmax = end_date,
        ymin = min(gg1$Returns) - 0.2 * (max(gg1$Returns) - min(gg1$Returns)),
        ymax = min(gg1$Returns) - 0.05 * (max(gg1$Returns) - min(gg1$Returns)),
      # ymin = min(gg1$Returns) - 20, ymax = min(gg1$Returns) - 1,
      fill = State), inherit.aes = F) + 
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

p2 <- ggplot(gg2, aes(x=Date, y=Returns, color=Stock))+
  geom_rect(data = state_int, 
    aes(xmin = start_date, xmax = end_date,
        ymin = min(gg2$Returns) - 0.2 * (max(gg2$Returns) - min(gg2$Returns)),
        ymax = min(gg2$Returns) - 0.05 * (max(gg2$Returns) - min(gg2$Returns)),
        # ymin = min(gg2$Returns) - 5, ymax = min(gg2$Returns) - 1,
        fill = State), inherit.aes = F) + 
  theme(legend.text = element_text(size=18), legend.title = element_text(size=22),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28), axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.key.size = unit(1.8, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0)+ geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + xlab("") +
  scale_fill_manual(values = personal_colors)

p2

p3 <- ggplot(gg3, aes(x=Date, y=Returns, color=Currency))+
  geom_rect(data = state_int, 
            aes(xmin = start_date, xmax = end_date,
                ymin = min(gg3$Returns) - 0.2 * (max(gg3$Returns) - min(gg3$Returns)),
                ymax = min(gg3$Returns) - 0.05 * (max(gg3$Returns) - min(gg3$Returns)),
                # ymin = min(gg3$Returns) - 2.5, ymax = min(gg3$Returns) - 1,
                fill = State), inherit.aes = F) + 
  theme(legend.text = element_text(size=18), legend.title = element_text(size=22),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28), axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.key.size = unit(1.8, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0)+ geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + xlab("") +
  scale_fill_manual(values = personal_colors)

p3

p4 <- ggplot(gg4, aes(x=Date, y=Returns, color=Energy))+
  geom_rect(data = state_int, 
            aes(xmin = start_date, xmax = end_date,
                ymin = min(gg4$Returns) - 0.2 * (max(gg4$Returns) - min(gg4$Returns)),
                ymax = min(gg4$Returns) - 0.05 * (max(gg4$Returns) - min(gg4$Returns)),
                # ymin = min(gg4$Returns) - 10, ymax = min(gg4$Returns) - 1,
                fill = State), inherit.aes = F) + 
  theme(legend.text = element_text(size=18), legend.title = element_text(size=22),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28), axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.key.size = unit(1.8, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0)+ geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + xlab("") +
  scale_fill_manual(values = personal_colors)

p4

p5 <- ggplot(gg5, aes(x=Date, y=Returns, color=Macroeconomy))+
  geom_rect(data = state_int, #barra degli stati
            aes(xmin = start_date, xmax = end_date,
                ymin = min(gg5$Returns) - 0.2 * (max(gg5$Returns) - min(gg5$Returns)),
                ymax = min(gg5$Returns) - 0.05 * (max(gg5$Returns) - min(gg5$Returns)),
                # ymin = min(gg5$Returns) - 5, ymax = min(gg5$Returns) - 1,
                fill = State), inherit.aes = F) + 
  theme(legend.text = element_text(size=18), legend.title = element_text(size=22),
        axis.title.x = element_text(size=28),
        axis.title.y = element_text(size=28), axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20),
        legend.key.size = unit(1.8, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(2.5, 'cm')
  ) +
  geom_point(alpha = 0)+ geom_line() +
  scale_x_date(date_labels = c("%Y"), breaks = "1 year") + xlab("") +
  scale_fill_manual(values = personal_colors)
p5

ggarrange(p1, p2, p3, p4, p5, nrow = 5, ncol = 1, legend = "top")
# 2000x1500

###### Heatmap ICL ######
mat <- t(HSMM.icl)
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
  geom_segment(data = df_min, aes(x = State - .5, xend = State + .5, y = Lambda, yend = Lambda), color = "red", size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_colorbar(barheight = unit(4.5, "inches"), barwidth = unit(0.3, "inches"))) +
  theme(legend.text = element_text(size=18), legend.title = element_blank(),
        axis.title.x = element_text(size=28, face = "italic"),
        axis.title.y = element_text(size=28), axis.text.x = element_text(size=20, angle = 0),
        axis.text.y = element_text(size=20),
        legend.key.size = unit(1.8, 'cm'), #change legend key size
        legend.key.height = unit(0.8, 'cm'), #change legend key height
        legend.key.width = unit(2.5, 'cm')
  ) 
#760x540


###################################################################################################################