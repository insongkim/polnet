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
install_github("insongkim/polnet", dependencies = TRUE, ref = "development")
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

We then choose the wildest row and column actors from the posterior estimates from the LSNM model (vi version). The wildest row and column actors are those with the largest/smallest x-coordinates and largest/smallest y-coordinates. By doing so, we can fix 4 row actors and 4 column actors. There is no panacea algorithm that works for all datasets. We will elaborate another possible algorithm to choose row/column actors to fix in a later example.

```r
# Find the Fixed Actors by Selecting the Wildest Actor, Regardless of Octants
res.array <- as.array(res.vb$stan_fitted_model)[,1,]
res.array.mean <- colMeans(res.array)

row.max <- 1
row.min <- 1
col.max <- 1
col.min <- 1

row.embeddings <- res.array.mean[grep("^row_embedding\\[\\d+,\\d+\\]$", names(res.array.mean))]
max.row.x.pool <- row.embeddings[grep("^row_embedding\\[\\d+,1\\]$", names(row.embeddings))]
max.row.x.index <- as.numeric(gsub("^row_embedding\\[(\\d+),1\\]","\\1",names(max.row.x.pool[order(max.row.x.pool,decreasing=T)[1:row.max]])))
min.row.x.index <- as.numeric(gsub("^row_embedding\\[(\\d+),1\\]","\\1",names(max.row.x.pool[order(max.row.x.pool,decreasing=F)[1:row.max]])))

max.row.y.pool <- row.embeddings[grep("^row_embedding\\[\\d+,2\\]$", names(row.embeddings))]
max.row.y.index <- as.numeric(gsub("^row_embedding\\[(\\d+),2\\]","\\1",names(max.row.y.pool[order(max.row.y.pool,decreasing=T)[1:row.max]])))
min.row.y.index <- as.numeric(gsub("^row_embedding\\[(\\d+),2\\]","\\1",names(max.row.y.pool[order(max.row.y.pool,decreasing=F)[1:row.max]])))

fixed_row_index <- sort(unique(c(max.row.x.index, min.row.x.index,
                                 max.row.y.index, min.row.y.index)))
N_fixed_row <- length(fixed_row_index)

col.embeddings <- res.array.mean[grep("^col_embedding\\[\\d+,\\d+\\]$", names(res.array.mean))]
max.col.x.pool <- col.embeddings[grep("^col_embedding\\[1,\\d+\\]$", names(col.embeddings))]
max.col.x.index <- as.numeric(gsub("^col_embedding\\[1,(\\d+)\\]","\\1",names(max.col.x.pool[order(max.col.x.pool,decreasing=T)[1:col.max]])))
min.col.x.index <- as.numeric(gsub("^col_embedding\\[1,(\\d+)\\]","\\1",names(max.col.x.pool[order(max.col.x.pool,decreasing=F)[1:col.max]])))

max.col.y.pool <- col.embeddings[grep("^col_embedding\\[2,\\d+\\]$", names(col.embeddings))]
max.col.y.index <- as.numeric(gsub("^col_embedding\\[2,(\\d+)\\]","\\1",names(max.col.y.pool[order(max.col.y.pool,decreasing=T)[1:col.max]])))
min.col.y.index <- as.numeric(gsub("^col_embedding\\[2,(\\d+)\\]","\\1",names(max.col.y.pool[order(max.col.y.pool,decreasing=F)[1:col.max]])))

fixed_col_index <- sort(unique(c(max.col.x.index, min.col.x.index,
                                 max.col.y.index, min.col.y.index)))
N_fixed_col <- length(fixed_col_index)
```

We then run the full MCMC model with the coordinates of these chosen actors fixed. 

```r
## Generate Data Matrice for Fixed Parameters
regex_row_embedding <- paste0("^row_embedding\\[(",paste0(fixed_row_index, collapse="|"),"),\\d+\\]$")
fixed_row_embedding <- res.array.mean[grep(regex_row_embedding, names(res.array.mean), perl = TRUE)]
regex_row_embedding_reshape_1 <- paste0("^row_embedding\\[(",paste0(fixed_row_index, collapse="|"),"),1\\]$")
regex_row_embedding_reshape_2 <- paste0("^row_embedding\\[(",paste0(fixed_row_index, collapse="|"),"),2\\]$")

fixed_row_embedding <- cbind(fixed_row_embedding[grep(regex_row_embedding_reshape_1, names(fixed_row_embedding))],
                             fixed_row_embedding[grep(regex_row_embedding_reshape_2, names(fixed_row_embedding))])

rownames(fixed_row_embedding) <- gsub(",\\d+\\]$","\\]",rownames(fixed_row_embedding))

regex_col_embedding <- paste0("^col_embedding\\[\\d+,(",paste0(fixed_col_index, collapse="|"),")\\]$")
fixed_col_embedding <- res.array.mean[grep(regex_col_embedding, names(res.array.mean), perl = TRUE)]
regex_col_embedding_reshape_1 <- paste0("^col_embedding\\[1,(",paste0(fixed_col_index, collapse="|"),")\\]$")
regex_col_embedding_reshape_2 <- paste0("^col_embedding\\[2,(",paste0(fixed_col_index, collapse="|"),")\\]$")

fixed_col_embedding <- cbind(fixed_col_embedding[grep(regex_col_embedding_reshape_1, names(fixed_col_embedding))],
                             fixed_col_embedding[grep(regex_col_embedding_reshape_2, names(fixed_col_embedding))])

rownames(fixed_col_embedding) <- gsub("\\d+,(\\d+)\\]$","\\1]",rownames(fixed_col_embedding))
fixed_col_embedding <- t(fixed_col_embedding)

res <- LSNM(sim.data$LSNM_data$A,
            N_fixed_row=N_fixed_row, N_fixed_col=N_fixed_col, 
            fixed_row_index=fixed_row_index, fixed_row_embedding=fixed_row_embedding,
            fixed_col_index=fixed_col_index, fixed_col_embedding=fixed_col_embedding,
            D=2, cores=7, warmup=1000, iter=2000, chains=4, control = list(max_treedepth = 20), 
            method="mcmc")
```

The following plot compares the MCMC posterior coordinates of actors with their true coordinates.
![](https://github.com/insongkim/repo-data/blob/master/polnet/lsnm_mcmc_true.png)
