# polnet:  A Statistical Analysis of Political Networks
[![Build Status](https://travis-ci.org/insongkim/wfe.svg?branch=master)](https://travis-ci.org/insongkim/polnet)
<!-- [![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/wfe)](https://cran.r-project.org/package=wfe) -->

This R package provides a computationally efficient way of fitting
the Latent Space Network Model (LSNM) developed by Kim and Kunisky (2018).

Authors
-------------------------
[In Song Kim](http://web.mit.edu/insong/www/), [Dmitry Kunisky], Jacob Jaffe

Paper
-------------------------
[Mapping Political Communities: A Statistical Analysis of Lobbying Networks in Legislative Politics](http://web.mit.edu/insong/www/pdf/network.pdf)

Installation
 -------------------------

There is currently no package available for the LNSM model on CRAN. It is possible to download the package through GitHub instead.

First, make sure the `devtools` package is installed.
``` r
if(!require(devtools)) install.packages("devtools")
```
This only has to be done if `devtools` is not already installed and only has to be done once.

Then, use the `install_github()` function from `devtools` to install the package.

``` r
library(devtools)
install_github("insongkim/polnet",dependencies=TRUE)
```

Example
 -------------------------
This example uses randomly generated poisson data:
```r
#Load polnet package
datmat = matrix(rpois(500, 1), ncol = 10)
lnsmR = LSNM(edges = datmat, D = 2)
 ```
 - `edges` is the matrix of connection strength data, a matrix or an object coercible to a matrix by `as.matrix()`
 -   `D` is the dimensionality of the latent space model, defaults to 2
