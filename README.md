Functional Brain Age Score (FBAS)
================

## Introduction

This funtion implements the methodology described in the paper

## Installation

This package can be installed with the 'devtools' package:

``` r
library(devtools) 
devtools::install_github(repo='hwiyoungstat/FBAS') 
```

## Usage

**`FBAS`** is the primary function that implements the FBAS based mediation analysis

``` r
FBAS(x, y, m, lambda, rho)
```

-   Inputs
    -   `x` : Input
    -   `y` : Exposure
    -   `m` : Mediator
    -   `lambda` : Regularizer
    -   `rho` : Lagrangian penalty
-   Output
    -   Estimated coefficients for FBAS

We generate a toy example :

-   The first 5 mediators among 20 are active (nonzero effects)

``` r
set.seed(20200314)
library(MASS)
N <- 50
P <- 20
Q <- 5
R <- 6

alpha <- c(rep(-0.5,Q),rep(0,P-Q))
beta <- c(rep(-0.5,R),rep(0,P-R))

Sds <- rep(1,P)
Cormat <- 0*(diag(1,P)==0)+diag(1,P)
Sig <- diag(Sds)%*%Cormat%*%diag(Sds)
  
X <- rnorm(N)
M <- X%*%t(alpha)+ mvrnorm(n = N, rep(0,P), Sigma=Sig)
Y <- M%*%beta + rnorm(N)
```

Result

``` r
result <-FBAS(X,Y,M,0.03,1)
which(result[[1]]!=0)
```

    ## [1] 1 2 3 4 5
