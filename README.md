# QuantRegLaguerre
This repository contains the R-Package "QuantRegLaguerre" which provides implementations for the quantile regression methodology used in the paper
"Efficient Quantile Regression under Censoring Using Laguerre Polynomials" written by Ingrid Van Keilegom and myself. The package contains all functions from the
estimator, performing cross-validation, to finding the bootstrap estimate for the variance. Therefore it can be directly used for the analysis of a given data set.
For the theory we refer to our preprint which is to appear soon.

To install the package in your local R environment from this repository you can use the following code.
```
library(remotes)
install_github("akreiss/QuantRegLaguerre")
library(QuantRegLaguerre)
```

Next, we generate a data set. In this example we generate data according to DGP in the preprint.

```
n   <- 200 # Number of observations
tau <- 0.5 # Quantile of the noise

beta <- c(2,1) # True parameters

set.seed(2023)
epsilon <- rnorm(n,mean=0,sd=1)-qnorm(tau)
x <- rnorm(n)
T <- beta[1]+x*beta[2]+(0.2+2*(x-0.5)^2)*epsilon
C <- runif(n,min=0,max=7)
Y <- pmin(T,C)
Delta <- as.numeric(T==Y)
```

To perform the estimation, we can either manually specify the degrees of the approximations, e.g., as follows
```
lag_est <- laguerre_estimator(1,1,1,Cov=x,Y=Y,Delta=Delta,tau=tau)
```
Here, the first three ones indicate that we used polynomials of degree at most 1 for all three approximations. As an alternative, we can also use cross-validation
to find suitable degrees:
```
CV_est <- laguerre_cross_validation(Y,Delta,tau,0.2,Cov=x)
```
Finally, to be able to perform inference, we use the boostrap to find confidence intervals. In order to use bootstrap in QuantRegLaguerre it is necessary to first
run cross-validation and the provide the result to the bootstrap functions. This looks as follows:
```
lag_boot <- laguerre_bootstrap(Y,Delta,matrix(x,ncol=1),tau,0.2,CV_est,level=0.1)
```
