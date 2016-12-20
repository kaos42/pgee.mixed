# sd without zeros (defined as abs < 1e-3)
nzsd <- function(x) {
  z <- x[abs(x) > 0.001]
  stats::sd(z)
}

# Sandwich formula for standard error, for mixed outcomes.
# X should be after trans_X_mixed
# Note that se's should be calculated using the stdized design matrices and
#   stdized coefficients. These se's need to be unstandardized.
sandwich.mixed <- function(y, X, Beta, alpha, lambda, N, wctype, eps,
                           intercept, phi, penalty, weights) {

  mu <- meanfn(X %*% Beta, family = "Mixed")
  v <- varfn(mu, family = "Mixed")
  H <- CppHess2(X, v, alpha, phi, N, weights)
  E <- Emat_wrap(Beta, lambda, "Mixed", eps, intercept, penalty)
  M <- CppM2(y, X, mu, v, alpha, phi, N, weights)
  p <- dim(H)[1]
  HE <- H + E + 1e-06 * diag(p)
  C <- solve(HE, M) %*% solve(HE)
}

# Sandwich formula for nonmmixed outcomes.
# As above, for standardized, hence results should be unstandardized later.
sandwich.nonmixed <- function(y, X, Beta, alpha, lambda, N, m, family,
                              wctype, eps, intercept, phi, penalty,
                              weights) {

  mu <- meanfn(X %*% Beta, family = family)
  v <- varfn(mu, family = family)
  R <- corrmat(alpha, m, wctype)
  H <- CppHess(X, v, R, phi, N, weights)
  E <- Emat_wrap(Beta, lambda, family, eps, intercept, penalty)
  M <- CppM(y, X, mu, v, R, phi, N, weights)
  p <- dim(H)[1]
  HE <- H + E + 1e-06 * diag(p)
  C <- solve(HE, M) %*% solve(HE)
}

# wrapper for all family types
# For family=="Mixed", X will be transformed. Else, X won't.
sandwich <- function(y, X, Beta, alpha, lambda, N, m, family,
                     wctype, eps, intercept, phi, penalty,
                     weights) {
  if (family != "Mixed") {
    sandwich.nonmixed(y, X, Beta, alpha, lambda, N, m, family,
                      wctype, eps, intercept, phi, penalty,
                      weights)
  } else {
    sandwich.mixed(y, X, Beta, alpha, lambda, N, wctype, eps,
                   intercept, phi, penalty, weights)
  }
}
