
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ExactMeanCov

<!-- badges: start -->

<!-- badges: end -->

`MomentMatching` ports the MATLAB package [Simulations with Exact Means
and
Covariances](https://la.mathworks.com/matlabcentral/fileexchange/24416-simulations-with-exact-means-and-covariances?s_tid=prof_contriblnk)
into R.

## Installation

You can install the development version of `MomentMatching` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Reckziegel/MomentMatching")
```

## Example

``` r
library(MomentMatching)
```

Assume that there are `20` risk-factors with a total of `200`
realizations from each risk-exposure.

``` r
N <- 20 # number of risk factors
J <- 200 # number of realizations/scenarios
```

Randomly generate a location parameter `M` and a scatter matrix `S` that
embeds the dispersion around the “true” parameter `M` and `S`.

``` r
set.seed(1234)

# ...vector of expected values M...
M <- matrix(runif(N), ncol = N) - 0.5

# ...covariance matrix S
A <- matrix(runif(N * N), ncol = N) - 0.5
S <- A %*% t(A)
```

Now, if we impose the additional assumption that the data generating
process (DGP) is normally distributed and we try to get some data from
the parameters at hand, we are left with:

``` r
# generate sample of size J from multivariate normal N~(M,S)
X <- MASS::mvrnorm(n = J, mu = M, Sigma = S) 
```

But this is just a blurred shadow of the “true” parameters generated
above, which can be seen by computing the absolute error of the location
and dispersion samples:

``` r
# Sample Moments
M_sample <- colMeans(X)
S_sample <- stats::cov(X)

max(abs(M - M_sample)) / max(abs(M)) # Sample Location Error
#> [1] 0.4293144
max(abs(S - S_sample)) / max(abs(S)) # Sample Dispersion Error
#> [1] 0.1227554
```

To minimize this error, Meucci (2009) suggests an affine transformation
that greatly helps to alleviate this problem:

``` r
# exact match between sample and population moments
X_ <- MvnRnd(M, S, J) 
```

Now, we can compute the moments of the simulated series `X_`, as we did
with `X`, to convince ourselves that Moment-Matching mechanism, indeed,
offers a simulation that better tracks the original data:

``` r
# Moments Matching
M_meucci <- colMeans(X_)
S_meucci <- (nrow(X_) - 1) / nrow(X_) * stats::cov(X_)

# Moment-Matching Errors
max(abs(M - M_meucci)) / max(abs(M)) # Location-Matching Error
#> [1] 0
max(abs(S - S_meucci)) / max(abs(S)) # Dispersion-Matching Error
#> [1] 2.539216e-15
```
