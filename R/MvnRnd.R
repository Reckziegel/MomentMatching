#' Simulations with Exact Means and Covariances
#'
#' This function generate scenarios from elliptical distributions in which the
#' sample moments match the populational moments.
#'
#' @param M A \code{1xN} matrix with the location parameter of each invariant (asset).
#' @param S A \code{NxN} matrix with the dispersion matrix of the invariants (assets).
#' @param J A numeric scalar with the desired number of scenarios.
#'
#' @return A J x N matrix with the simulated time-series.
#'
#' @export
#'
#' @references
#' Meucci, Attilio, Simulations with Exact Means and Covariances (June 7, 2009). Available at SSRN: \url{https://ssrn.com/abstract=1415699} or \url{http://dx.doi.org/10.2139/ssrn.1415699}.
#'
#' Attilio Meucci (2020). Simulations with Exact Means and Covariances, \href{https://www.mathworks.com/matlabcentral/fileexchange/24416-simulations-with-exact-means-and-covariances}{MATLAB Central File Exchange}. Retrieved October 11, 2020.
#'
#' @examples
#' rets  <- diff(log(EuStockMarkets))
#' mu    <- colMeans(rets)
#' sigma <- stats::cov(rets)
#' MvnRnd(M = mu, S = sigma, J = 50)
MvnRnd <- function(M, S, J) {

  if (!(is.numeric(J) && length(J) == 1)) {
    stop("J is not a number (a length one numeric vector).")
  }

  if (!is.matrix(M)) {
    M <- t(as.matrix(M))
  }

  if (nrow(M) == ncol(S)) {
    M <- t(M)
  }

  N <- length(M)


  # generate antithetic variables (mean = 0)
  Y <- MASS::mvrnorm(n = J / 2, mu = matlab::zeros(N, 1), Sigma = S)
  Y <- rbind(Y, -Y)

  # compute sample covariance: NOTE defined as "cov(Y,1)", not as "cov(Y)"
  S_ <- (nrow(Y) - 1) / nrow(Y) * stats::cov(Y)

  # solve Riccati equation using Schur method
  H <- rbind(
    cbind(matlab::zeros(N,N), -S_),
    cbind(-S, matlab::zeros(N,N))
  )

  U <- QZ::ordqz(H, keyword = "lhp")$Q

  U_lu <- U[1:N, 1:N]
  U_ld <- U[(N + 1):nrow(U), 1:N]

  B <- U_ld %*% solve(U_lu)

  # affine transformation to match mean and covariances
  X <- Y %*% B + matlab::repmat(t(M), J, 1)

  X

}
