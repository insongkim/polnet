# polnet:  A Statistical Analysis of Political Networks
[![Build Status](https://travis-ci.org/sooahn/polnet.svg?branch=development)](https://travis-ci.org/sooahn/polnet)

This R package provides a computationally efficient way of fitting
the Latent Space Network Model (LSNM) and bipartite Link Community Model (biLCM) developed by Kim and Kunisky (2018).

Authors
-------------------------
[In Song Kim](http://web.mit.edu/insong/www/), Dmitriy Kunisky, Sean (Shiyao) Liu, [Sooahn Shin](http://sooahnshin.com/)

Paper
-------------------------
[Mapping Political Communities: A Statistical Analysis of Lobbying Networks in Legislative Politics](http://web.mit.edu/insong/www/pdf/network.pdf)

Installation
 -------------------------

There is currently no package available for the LNSM or biLCM models on CRAN. It is possible to download the package through GitHub instead.

First, make sure the `devtools` package is installed.
``` r
if(!require(devtools)) install.packages("devtools")
```
This only has to be done if `devtools` is not already installed and only has to be done once.

Then, use the `install_github()` function from `devtools` to install the package.

``` r
library(devtools)
install_github("insongkim/polnet", dependencies = TRUE, ref = "development")
```

Examples
-------------------------

The example of LSNM (with variational inference) using synthetic data is as follows. This is useful in the case where the researcher is interested in identifying the 'clusters' of the units.

``` r
library(polnet)
set.seed(11)
sim.data <- random_LSNM_data_cluster(n.cluster=4, group1.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, group2.center=rbind(c(-0.5,-1), c(-1, 0.3), c(0.4, 1), c(0.2, -0.2))*5, v=3, sigma_sq_L = 0.5, sigma_sq_P = 0.7, tau=c(0.5, 0.8))
res <- LSNM(sim.data$LSNM_data$A, D=2, iter=50000)
plot.compare.LSNM(res, sim.data$LSNM_data$Theta, sim.data$LSNM_data$Psi, sim.data$group1.popularity, sim.data$group2.popularity, legend_position = "center")
```
