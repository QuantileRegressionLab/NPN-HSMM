rm(list = ls())
gc()
library(dplyr)
library(stringr)

# ---- 1) EXPECTED INPUT SHAPE -----------------------------------------------
# Provide a data frame 'df' with one row per (sojourn, dist, T).
# Columns required (means + ses):
# sojourn: one of c("Shifted Poisson","Shifted Negative Binomial","Geometric")
# dist   : one of c("Multivariate Gaussian","Multivariate Gaussian with outliers")
# T      : e.g., 300, 500, 1000
# ARI_mean, ARI_se,
# M1_mean, M1_se,   # ||M_1||_F
# M2_mean, M2_se,   # ||M_2||_F
# M3_mean, M3_se,   # ||M_3||_F
# TNR1_mean, TNR1_se, TNR2_mean, TNR2_se, TNR3_mean, TNR3_se,
# TPR1_mean, TPR1_se, TPR2_mean, TPR2_se, TPR3_mean, TPR3_se

example_df <- function(df) {
  df <- as.data.frame(df)
  head1 <- expand.grid(
    sojourn = c("Shifted Poisson","Shifted Negative Binomial","Geometric"),
    dist    = c("Multivariate Gaussian","Multivariate Gaussian with outliers"),
    T       = c(300, 500, 1000),
    KEEP.OUT.ATTRS = FALSE
  ) %>%
    arrange(sojourn, dist, T)
  
  data.frame(head1, df)
  
}

# ---- 2) FORMATTER -----------------------------------------------------------
fmt <- function(m, s, digits = 3) {
  paste0(formatC(m, format = "f", digits = digits),
         " (",
         formatC(s, format = "f", digits = digits),
         ")")
}

# ---- 3) TABLE BUILDER -------------------------------------------------------
make_sim_table <- function(df,
                           sojourn_order = c("Shifted Poisson",
                                             "Shifted Negative Binomial",
                                             "Geometric"),
                           dist_order    = c("Multivariate Gaussian",
                                             "Multivariate Gaussian with outliers"),
                           file_tex = NULL) {
  
  stopifnot(all(c("sojourn","dist","T",
                  "ARI_mean","ARI_se",
                  "M1_mean","M1_se",
                  "M2_mean","M2_se",
                  "M3_mean","M3_se",
                  "TNR1_mean","TNR1_se","TNR2_mean","TNR2_se","TNR3_mean","TNR3_se",
                  "TPR1_mean","TPR1_se","TPR2_mean","TPR2_se","TPR3_mean","TPR3_se") %in% names(df)))
  
  df <- df %>%
    mutate(
      sojourn = factor(sojourn, levels = sojourn_order),
      dist    = factor(dist,    levels = dist_order)
    ) %>%
    arrange(sojourn, dist, T) %>%
    mutate(
      ARI  = fmt(ARI_mean,  ARI_se),
      M1   = fmt(M1_mean,   M1_se),
      M2   = fmt(M2_mean,   M2_se),
      M3   = fmt(M3_mean,   M3_se),
      TNR1 = fmt(TNR1_mean, TNR1_se),
      TNR2 = fmt(TNR2_mean, TNR2_se),
      TNR3 = fmt(TNR2_mean, TNR3_se),
      TPR1 = fmt(TPR1_mean, TPR1_se),
      TPR2 = fmt(TPR2_mean, TPR2_se),
      TPR3 = fmt(TPR3_mean, TPR3_se)
    ) %>%
    select(sojourn, dist, T, ARI, M1, M2, M3, TNR1, TNR2, TNR3, TPR1, TPR2, TPR3)
  
  # Build the body by inserting \multicolumn header rows exactly like your LaTeX
  rows_list <- list()
  hline_after <- integer(0)  # collect row indices to add \hline after subheaders
  cur_row <- 0
  
  for (s in levels(df$sojourn)) {
    # bold left-aligned block header
    rows_list[[length(rows_list)+1]] <-
      data.frame(
        T    = sprintf("\\multicolumn{11}{l}{\\textbf{%s}}", s),
        ARI  = "", M1="", M2="", M3="", TNR1="", TNR2="", TNR3="", TPR1="", TPR2="", TPR3="",
        stringsAsFactors = FALSE
      )
    cur_row <- cur_row + 1
    
    df_s <- df %>% filter(sojourn == s)
    
    for (d in levels(df$dist)) {
      df_sd <- df_s %>% filter(dist == d)
      if (nrow(df_sd) == 0) next
      
      # centered subheader
      rows_list[[length(rows_list)+1]] <-
        data.frame(
          T    = sprintf("\\multicolumn{11}{c}{%s}", d),
          ARI  = "", M1="", M2="", M3="", TNR1="", TNR2="", TNR3="", TPR1="", TPR2="", TPR3="",
          stringsAsFactors = FALSE
        )
      cur_row <- cur_row + 1
      hline_after <- c(hline_after, cur_row)  # add \hline below this row
      
      # actual data rows
      body <- df_sd %>%
        transmute(
          T = as.character(T), ARI, M1, M2, M3, TNR1, TNR2, TNR3, TPR1, TPR2, TPR3
        )
      rows_list[[length(rows_list)+1]] <- body
      cur_row <- cur_row + nrow(body)
    }
  }
  
  out <- bind_rows(rows_list)
  
  # Build kable
  col_names <- c("$T$","ARI","$||M_1||_F$","$||M_2||_F$","$||M_3||_F$",
                 "$\\text{TNR}_1$","$\\text{TNR}_2$","$\\text{TNR}_3$","$\\text{TPR}_1$","$\\text{TPR}_2$","$\\text{TPR}_3$")
  
  kb <- knitr::kable(
    out,
    format = "latex",
    booktabs = TRUE,
    align = "c",
    col.names = col_names,
    escape = FALSE,
    linesep = ""
  )
  
  if (!is.null(file_tex)) {
    writeLines(as.character(kb), file_tex)
  }
  kb
}

# ---- 4) RUN -----------------------------------------------
N = c(300, 500, 1000)
sojourn = c("dPoisson", "dNegBin", "dGeom")
error <- c("eMVNorm", "eOutliers")
idx = expand.grid(N, error, sojourn)
idxl = nrow(idx)
idx = c(idx)
myres <- matrix(NA, idxl, 20)

for (i in 1:idxl) {
  file_path = paste0("Sim/N", paste(idx$Var1[i], idx$Var3[i], idx$Var2[i], "K3.csv", sep = "_"))
  tmp = read.csv(file_path)
  df.mean = as.data.frame(t(apply(tmp, 2, mean)))
  df.sd = as.data.frame(t(apply(tmp, 2, sd)))
  myres[i,] = c(df.mean$ARI_HSMM, df.sd$ARI_HSMM, df.mean$M_1_HSMM, df.sd$M_1_HSMM, df.mean$M_2_HSMM, df.sd$M_2_HSMM,
                df.mean$M_3_HSMM, df.sd$M_3_HSMM, df.mean$TFP_1_HSMM, df.sd$TFP_1_HSMM, df.mean$TFP_2_HSMM, df.sd$TFP_2_HSMM,
                df.mean$TFP_3_HSMM, df.sd$TFP_3_HSMM, df.mean$TFN_1_HSMM, df.sd$TFN_1_HSMM, df.mean$TFN_2_HSMM, df.sd$TFN_2_HSMM,
                df.mean$TFN_3_HSMM, df.sd$TFN_3_HSMM)
}

names.myres = c("ARI", "M1", "M2", "M3", "TNR1", "TNR2", "TNR3", "TPR1", "TPR2", "TPR3")
colnames(myres) = sapply(1:length(names.myres), function(i) paste(names.myres[i], c("mean", "se"), sep = "_"))


df <- example_df(myres)        # demo
tbl <- make_sim_table(df)
tbl   # prints LaTeX to console and writes a .tex file

###########################################
###########################################
###########################################
rm(list = ls())
gc()

N = c(300, 500, 1000)
sojourn = c("dPoisson", "dNegBin", "dGeom")
error <- c("eMVNorm", "eOutliers")
idx = expand.grid(N, error, sojourn)
idxl = nrow(idx)
idx = c(idx)
myres <- c()

for (i in 1:idxl) {
  file_path = paste0("Sim/N", paste(idx$Var1[i], idx$Var3[i], idx$Var2[i], "K3.csv", sep = "_"))
  tmp = read.csv(file_path)
  tmp.vec = c(mean(tmp$t.time_HSMM), mean(tmp$t.iter_HSMM))
  myres = cbind(myres, tmp.vec)
}

myres = matrix(t(myres), nrow = 3, byrow = F)
myres = cbind(N, myres)
xtable(myres, digits = 3)


names.myres = c("ARI", "M1", "M2", "M3", "TNR1", "TNR2", "TNR3", "TPR1", "TPR2", "TPR3")
colnames(myres) = sapply(1:length(names.myres), function(i) paste(names.myres[i], c("mean", "se"), sep = "_"))


df <- example_df(myres)        # demo
tbl <- make_sim_table(df)
tbl   # prints LaTeX to console and writes a .tex file

###########################################
###########################################
###########################################
rm(list = ls())
gc()
library(dplyr)
library(stringr)
library(Matrix)
library(xtable)

N = c(1500, 3000)
sojourn = c("dPoisson", "dNegBin")
states <- c("K2", "K3")
idx = expand.grid(N, states, sojourn)
idxl = nrow(idx)
idx = c(idx)
cov.list <- list()

for (i in 1:idxl) {
  file_path = paste0("Sim/N", paste(idx$Var1[i], idx$Var3[i], "eMVNorm", idx$Var2[i], "boot.RData", sep = "_"))
  load(file_path)
  sojourn.distribution = ifelse(idx$Var3[i] == "dPoisson", "poisson", "dbinom")
  tmp = coverage.summary(fit = results_list, boot = results_list.boot, theta = theta, d_true = d_true, gamma_sim = gamma_sim, sojourn.distribution = sojourn.distribution, alpha = 0.05)
  cov.list[[i]] = tmp

  print(paste(idx$Var1[i], idx$Var3[i], "eMVNorm", idx$Var2[i]))
  # xtable(t(matrix(cov.list[[i]])), digits = 3)

  gc()
  rm("results_list", "results_list.boot")
}
save(cov.list, file = "coverage.RData")

for (i in 1:idxl) {
  print(paste(idx$Var1[i], idx$Var3[i], "eMVNorm", idx$Var2[i]))
  cov.gamma.off = cov.list[[i]]$coverage_gamma
  cov.tmp = c(idx$Var1[i], cov.gamma.off[upper.tri(cov.gamma.off)], cov.gamma.off[lower.tri(cov.gamma.off)], cov.list[[i]]$coverage_d, cov.list[[i]]$coverage_omega)
  print(xtable(t(matrix(cov.tmp)), digits = 2))
}
