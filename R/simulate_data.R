#' Generate correlated bivariate mixed outcome data
#'
#' \code{gen_mixed_data} returns randomly generated correlated bivariate mixed
#' outcomes, and covariate matrices to model them, based on design parameters
#' set in the function.
#'
#' A Gaussian copula is used to generate the correlated outcomes. Marginally,
#' the continuous outcome follows a normal distribution with identity link to
#' covariates, while the binary outcome follows a Bernoulli distribution with
#' logit link to covariates. Covariates are generated from a zero-mean unit
#' variance multivariate normal distribution, with an AR(1) correlation
#' structure.
#'
#' @param Beta.cont Vector of true regression coefficients for the continuous
#'   outcome.
#' @param Beta.bin Vector of true regression coefficients for the binary
#'   outcome.
#' @param N Number of pairs of correlated outcomes.
#' @param rho Gaussian copula parameter.
#' @param intercept Assume an intercept (for both outcomes)? (default TRUE). If
#'   TRUE, then the first coefficient in Beta.cont and Beta.bin are assumed to
#'   correspond to intercepts.
#' @param cov Specify if the covariate matrices for the continuous outcome and
#'   the binary outcome should share all covariates (set to "same"), share some
#'   covariates (set to "shared"), or share no covariates (set to "separate").
#' @param xcor Correlation parameter for AR(1) correlation structure of
#'   covariate design matrices (assumed same for both).
#' @param sigma_yc Marginal variance of continuous responses.
#' @return A list of generated data
#' \item{yc}{Vector of continuous outcomes.}
#' \item{yb}{Vector of binary outcomes.}
#' \item{X}{Covariate matrix for the continuous outcomes.}
#' \item{Z}{Covariate matrix for the binary outcomes.}
#' @examples
#' # default settings
#' gen_mixed_data(rnorm(5), rnorm(5), 10, 0.5)
#' # separate covariate matrices, non-unit continuous variance
#' gen_mixed_data(rnorm(5), rnorm(5), 10, 0.5, cov = "separate", sigma_yc = 2)
#' @export gen_mixed_data
gen_mixed_data <- function(Beta.cont, Beta.bin, N, rho, intercept = TRUE,
                           cov = "same", xcor = 0.25, sigma_yc = 1) {

  p <- length(Beta.cont)
  q <- length(Beta.bin)
  if (p < 5 | q < 5) stop("Please specify at least five coefficients per outcome.")
  # Generate covariates
  W1 <- mvtnorm::rmvnorm(N, mean = rep(0, p), sigma = corrmat(xcor, p, "AR1"))
  W2 <- mvtnorm::rmvnorm(N, mean = rep(0, q), sigma = corrmat(xcor, q, "AR1"))
  w3 <- stats::rnorm(N)
  w4 <- stats::rnorm(N)
  if (cov == "separate") {
    X <- W1
    Z <- W2
  } else if (cov == "same") {
    X <- Z <- W2
  } else if (cov == "shared") {
    X <- Z <- cbind(w4, w3, W1[, 1])
    X <- cbind(X, W1[, 2:(p - 2)])  # X = [w4 w3 w1(1) w1(2)...w1(p-2)]
    Z <- cbind(Z, W2[, 1:(q - 3)])  # Z = [w4 w3 w1(1) w2(1)...w2(q-3)]
  } else stop("Error in gen_mixed_data: Please set cov to either 'same',
              'shared' or 'separate'")
  if (intercept) X[, 1] <- Z[, 1] <- 1

  # Generate correlated uniform samples from Gaussian copula. Then apply
  # inverse normal cdf to one, inverse logistic cdf to other.
  # Cutoff logistic at 0 to get the binary outcome.
  norm.cop <- copula::normalCopula(rho, dim = 2, dispstr = "ex")
  # y <- numeric(N * 2)
  yc <- yb <- numeric(N)
  for (i in 1:N) {
    copdata <- copula::rCopula(1, norm.cop)
    # y[2 * i - 1] <- qnorm(copdata[1], mean = crossprod(X[i, ], Beta.cont),
    #                       sd = sigma_yc)
    # y[2 * i] <- (qlogis(copdata[2], location = crossprod(Z[i, ], Beta.bin)) >
    #                0) * 1  # 0 1 coding
    yc[i] <- stats::qnorm(copdata[1],
                          mean = crossprod(X[i, ], Beta.cont), sd = sigma_yc)
    yb[i] <- (stats::qlogis(copdata[2],
                            location = crossprod(Z[i, ], Beta.bin)) > 0) * 1
  }
  # Return generated data as list
  l <- list(yc = yc, yb = yb, X = X, Z = Z)
}
