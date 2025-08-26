# README file

The `R` scripts in this repository are designed to evaluate the performance of the research detailed in *Nonparanormal hidden semi-Markov graphical models for analyzing financial markets interconnectivity* by Emilio Ferrante, Foroni, Merlo, and Petrella (2025). 


## Prerequisites
### Software Requirements

-   [R](https://cran.r-project.org/) version 4.5.1 or higher
-   [RStudio](https://rstudio.com/) version 2024.12.1+563 or higher

### R Packages used (version in parentheses)

- MASS (7.3.65)
- mvtnorm (1.3.3)
- foreach (1.5.2)
- doParallel (1.0.17)
- parallel (4.5.1)
- mclust (6.1.1)
- cluster (2.1.8)
- markovchain (0.10.0)
- mhsmm (0.4.21)
- glasso (1.11)
- car (3.1.3)
- Matrix (1.7.3) 
- ggplot2 (3.5.2)

## Script description
### MainFunctions.R
This code contains the functions to run the `R` scripts in this repository. The main functions are: 
-    `hsmm.multi.gen` allows to simulate data from the proposed nonparanormal hidden Semi-Markov graphical model.
-    `EM_HSMM` contains the EM algorithm to fit the model for a given number of hidden states and value of the Lasso regularization parameter.
-    `boot.hsmm.multi.gen` allows to sample from a fitted nonparanormal model (used for parametric bootstrap).

### Simulations
The following scripts allows to reproduce the results of our simulation study:
-    `run_simulazioni` and `run_simulazioni_3K` contain the code to run all the considered scenarios in our simulation study with 2 and 3 hidden states.
-    `run_simulazioni_glasso` and `run_simulazioni_glasso_3K` contain the code to fit the hidden Markov Normal graphical model of Städler and Mukherjee (2013). for all the considered scenarios in our simulation study with 2 and 3 hidden states.
