# polnet:  A Statistical Analysis of Political Networks
[![Build Status](https://travis-ci.org/insongkim/polnet.svg?branch=master)](https://travis-ci.org/insongkim/polnet)

This R package provides a computationally efficient way of fitting
the Latent Space Network Model (LSNM) and bipartite Link Community Model (biLCM) developed by Kim and Kunisky (2018).

## Authors
[In Song Kim](http://web.mit.edu/insong/www/), Dmitriy Kunisky, [Sean (Shiyao) Liu](https://polisci.mit.edu/people/sean-shiyao-liu), [Sooahn Shin](http://sooahnshin.com/)

## Paper
[Mapping Political Communities: A Statistical Analysis of Lobbying Networks in Legislative Politics](http://web.mit.edu/insong/www/pdf/network.pdf)

## Installation

There is currently no package available for the LNSM or biLCM models on CRAN. It is possible to download the package through GitHub instead.

First, make sure the `devtools` package is installed.
``` r
if(!require(devtools)) install.packages("devtools")
```
This only has to be done if `devtools` is not already installed and only has to be done once.

Then, use the `install_github()` function from `devtools` to install the package.

``` r
library(devtools)
install_github("insongkim/polnet", dependencies = TRUE, ref = "master")
```

## Examples

### LSNM with variational inference

The example of LSNM (with variational inference) using synthetic data is as follows. The largest advantage for running a variational inference version of the LSNM model is its speed. As a result, this becomes useful in the case where the researcher is interested in identifying the 'clusters' of the units.

``` r
library(polnet)
set.seed(11)
sim.data <- random_LSNM_data_cluster(n.cluster=4, 
                                     group1.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5,    
                                     group2.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, 
                                     v=3, sigma_sq_L = 0.5, sigma_sq_P = 0.7, tau=c(0.5, 0.8))
res <- LSNM(sim.data$LSNM_data$A, D=2, method = "vi", iter=50000)
plot.compare.LSNM(res, sim.data$LSNM_data$Theta, sim.data$LSNM_data$Psi, sim.data$group1.popularity, 
                  sim.data$group2.popularity, sim.data$group1.cluster, sim.data$group2.cluster, 
                  legend_position = "center")
```
![](https://github.com/insongkim/repo-data/blob/master/polnet/lsnm_short_ex.png)

### LSNM with MCMC

If a researcher wants to do a full MCMC procedure with her data. We recommend the following procedure:

1. Run a variational inference version of the LSNM model. 

2. Fix the coordinates of 2 * d actors and run the full MCMC model.

The reason to fix 2 * d actors is to achieve identification by removing rotational/flipping/reflectional invariance, where d is the dimension of the model. If 2 * d cannot achieve identification for your data, you may increase the number of fixed actors until identification has been achieved.

We now illustrate this procedure with a synthetic data. We first run a variational inference version of LSNM on the synthetic dataset.

```r
library("polnet")
set.seed(11)
# Generate Synthetic Data
sim.data <- random_LSNM_data_cluster(n.cluster=4, 
                                     group1.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, 
                                     group2.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, 
                                     v=3, sigma_sq_L = 0.5, sigma_sq_P = 0.7, tau=c(0.5, 0.8))

# Quick Estimation with Variational Inference 
res.vb <- LSNM(sim.data$LSNM_data$A, D=2, method = "vi", iter=50000)
plot.compare.LSNM(res.vb, sim.data$LSNM_data$Theta, 
                  sim.data$LSNM_data$Psi, sim.data$group1.popularity, 
                  sim.data$group2.popularity, sim.data$group1.cluster, 
                  sim.data$group2.cluster, legend_position = "center")
```
![](https://github.com/insongkim/repo-data/blob/master/polnet/lsnm_vb_true.png)

We then choose the wildest row and column actors from the posterior estimates from the LSNM model (vi version). There are different ways to do this choice. We have embedded two algorithms in our choose.fix function. Here we use the o

```r
fix.actors <- choose.fix(res.vb, choose.method="axis")
```

We then run the full MCMC model with the coordinates of these chosen actors fixed. 

```r
res.mcmc <- LSNM(sim.data$LSNM_data$A, D=2, 
          method = "mcmc", iter=2000, 
          fixed.actor.object = fix.actors, 
          cores=4, control = list(max_treedepth = 20))
```

The following plot compares the MCMC posterior coordinates of actors with their true coordinates.
![](https://github.com/insongkim/repo-data/blob/master/polnet/lsnm_mcmc_true.png)


### biLCM

The example of biLCM using synthetic data is as follows.

``` r
set.seed(11)
sim.data <- random_biLCM_data(m = 100, n = 50, k = 4, kappa_weight = NULL, a = 10000, b = 1, alpha_membership = NULL, alpha_c = rep(0.5, 100), beta_membership = NULL, beta_c = rep(0.5, 50), non_zero = TRUE)
res <- biLCM(edges = sim.data$A, group1.id = NULL, group2.id = NULL, count.id = NULL, k = 4, tolerance = 1e-6, max.iter = 200)
plot.compare.biLCM(res, sim.data, group1 = TRUE, nth = 10)
```

The following plot compares true and estimated community distribution of the 10th actor in the first group (i.e., $(\alpha_{10,1},\alpha_{10,2},\alpha_{10,3},\alpha_{10,4})$). Note that the order of communities is permutation-invariant.
![](https://github.com/insongkim/repo-data/blob/master/polnet/bilcm_compare.png?raw=true)