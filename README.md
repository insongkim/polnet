# polnet:  A Statistical Analysis of Political Networks
[![Build Status](https://travis-ci.org/sooahn/polnet.svg?branch=test)](https://travis-ci.org/sooahn/polnet)

This R package provides a computationally efficient way of fitting
the Latent Space Network Model (LSNM) and bipartite Link Community Model (biLCM) developed by Kim and Kunisky (2018).

Authors
-------------------------
[In Song Kim](http://web.mit.edu/insong/www/), Dmitriy Kunisky, Jacob Jaffe

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
