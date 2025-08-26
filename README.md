# README file

The `R` scripts in this repository are designed to evaluate the performance of the research detailed in *Nonparanormal hidden semi-Markov graphical models for analyzing financial markets interconnectivity* by E. Ferrante, B. Foroni, L. Merlo, and L. Petrella (2025). 


## Prerequisites
### Software Requirements

-   [R](https://cran.r-project.org/) version 4.5.1 or higher
-   [RStudio](https://rstudio.com/) version 2024.12.1+563 or higher

Some portions of the code have been parallelized using the following packages:
-   `foreach`
-   `parallel`
-   `doParallel`

### R Packages used (version in parentheses)

- MASS (7.3.65)
- mvtnorm (1.3.3)
- foreach (1.5.2)
- parallel (4.5.1)
- doParallel (1.0.17)
- mclust (6.1.1)
- cluster (2.1.8)
- markovchain (0.10.0)
- mhsmm (0.4.21)
- glasso (1.11)
- car (3.1.3)
- Matrix (1.7.3)
- dplyr (1.1.4)
- stringr (1.5.1)
- doBy (4.6.26)
- qgraph (1.9.8)
- igraph (2.1.4)
- viridis (0.6.5)
- ggplot2 (3.5.2)
- xts (0.14.1)
- reshape2 (1.4.4)
- gridExtra (2.3)
- ggpubr (0.6.0)

## Script description
### MainFunctions.R
This code contains the functions to run the `R` scripts in this repository. The main functions are: 
-    `hsmm.multi.gen` allows to simulate data from the proposed nonparanormal hidden Semi-Markov graphical model.
-    `EM_HSMM` contains the EM algorithm to fit the model for a given number of hidden states and value of the Lasso regularization parameter.
-    `boot.hsmm.multi.gen` allows to sample from a fitted nonparanormal model (used for parametric bootstrap).

### Simulations
The following scripts allows to reproduce the results of our simulation study:
-    `run_sim` and `run_sim_3K` contain the code to run all the considered scenarios in our simulation study with 2 and 3 hidden states.
-    `run_sim_glasso` and `run_sim_glasso_3K` contain the code to fit the hidden Markov Normal graphical model of Städler and Mukherjee (2013) for all the considered scenarios in our simulation study with 2 and 3 hidden states. The EM algorithm to implement this model can be found in the script `em_glasso`.
-    `run_sim_bootstrap` and `run_sim_bootstrap_3K` contain the code to obtain the confidence interval coverage based on parametric bootstrap, with 2 and 3 hidden states.
-    `tab_sim` reproduces the tables in the Supplementary Materials summarizing the results.

### Empirical application
-    `emp_analysis_2023` allows to implement the procedure for selecting simultaneously the number of latent states and the optimal value of the Lasso regularization parameter using the ICL. This script yields the final model with 3 states, and can be used to reproduce Figure 1 in the manuscript.
-    `emp_analysis_2023_boot` implements the parametric bootstrap approach described in the Supplementary Materials of the work. This script requires the output provided by `emp_analysis_2023.R` and yields the estimated standard deviation of each parameter estimates.
-    Once we obtained the estimates and standard errors, the script `fig&tab_emp.R` allows to reproduce Figures 1 (ICL Heatmap), 2 (state-colored time series of returns), 3 (sojourn distributions) and 4 (state-specific graphs) in the manuscript.
