# library(Matrix)

# Mean (inverse link) function
meanfn <- function(eta, family) {
  if (family == "Gaussian") {
    mu <- eta
  } else if (family == "Binomial") {
    mu <- 1/(1 + exp(-eta))
  } else if (family == "Mixed") {
    # Separate out bin and cont, apply mean funcs, recombine
    id <- 1:length(eta)
    id_b <- (id%%2 == 0)
    id_c <- (id%%2 != 0)
    eta_c <- eta[id_c]
    eta_b <- eta[id_b]
    mu_c <- eta_c
    mu_b <- 1/(1 + exp(-eta_b))
    mu <- numeric(length(eta))
    mu[id_c] <- mu_c
    mu[id_b] <- mu_b
  } else stop("Error in meanfn: Wrong family specification")
  as.vector(mu)
}

# Variance function
varfn <- function(mu, family) {
  if (family == "Gaussian") {
    v <- rep(1, length(mu))
  } else if (family == "Binomial") {
    v <- mu * (1 - mu)
    # For numerical stability
    v <- ifelse(v < 1e-4, 1e-4, v)
  } else if (family == "Mixed") {
    # Separate out bin and cont, apply mean funcs, recombine
    id <- 1:length(mu)
    id_b <- (id%%2 == 0)
    id_c <- (id%%2 != 0)
    mu_c <- mu[id_c]
    mu_b <- mu[id_b]
    v_c <- rep(1, length(mu_c))
    v_b <- mu_b * (1 - mu_b)
    v_b <- ifelse(v_b < 1e-4, 1e-4, v_b)  # For numerical stability
    v <- numeric(length(mu))
    v[id_c] <- v_c
    v[id_b] <- v_b
  } else stop("Error in varfn: Wrong family specification")
  v
}

# penalty derivaritve functions, q
# Vectorized below to q_SCAD, q_MCP, etc. functions
foo_SCAD <- function(theta, lambda, a = 3.7) {
  if (a <= 2)
    stop("Error in foo_SCAD: \"a\" must be > 2")
  if (lambda < 0)
    stop("Error in foo_MCP:lambda must be >= 0")
  if (theta < 0) {
    err <- c("Error in foo_SCAD: ", "Illegal value for theta in SCAD penalty. ",
             "Make sure theta >= 0")
    stop(paste(err, sep = ""))
  }
  if (theta > lambda) {
    Q <- max((a * lambda - theta), 0)/(a - 1)
  } else {
    Q <- lambda
  }
}

q_SCAD <- Vectorize(foo_SCAD, "theta")

foo_MCP <- function(theta, lambda, a = 3) {
  if (a <= 1)
    stop("Error in foo_MCP: \"a\" must be > 1")
  if (lambda < 0)
    stop("Error in foo_MCP: lambda must be >= 0")
  if (theta <= a * lambda) {
    Q <- lambda - theta/a
  } else {
    Q <- 0
  }
}

q_MCP <- Vectorize(foo_MCP, "theta")

q_LASSO <- function(theta, lambda) {
  if (lambda < 0) stop("Error in foo_MCP:lambda must be >= 0")
  rep(lambda, length(theta))
}

# Return diagonal (vector) of the penalty matrix E.
# For family == "Mixed", Beta must be stacked as c(Beta_c, Beta_b)
Emat <- function(Beta, lambda, eps, penalty) {
  if (penalty == "SCAD") {
    E <- q_SCAD(abs(Beta), lambda = lambda)/(eps + abs(Beta))
  } else if (penalty == "MCP") {
    E <- q_MCP(abs(Beta), lambda = lambda)/(eps + abs(Beta))
  } else if (penalty == "LASSO") {
    E <- q_LASSO(abs(Beta), lambda = lambda)/(eps + abs(Beta))
  } else stop("error in Emat: Unknown penalty selected")
}

# Return the diagonal matrix of penalties E
# For family == "Mixed", lambda = c(lambda_c, lambda_b) must be a
#   vector of length 2
Emat_wrap <- function(Beta, lambda, family, eps, intercept, penalty) {
  if (family != "Mixed") {
    E <- diag(Emat(Beta, lambda, eps, penalty))
    if (intercept)
      E[1, 1] <- 0
    E
  } else {
    # Beta = (Beta_c, Beta_b)
    k <- length(Beta)
    Beta_c <- Beta[1:(k/2)]
    Beta_b <- Beta[(k/2 + 1):k]
    # lambda = (lambda_c, lambda_b)
    E_c <- Emat(Beta_c, lambda[1], eps, penalty)
    E_b <- Emat(Beta_b, lambda[2], eps, penalty)
    if (intercept) E_c[1] <- E_b[1] <- 0  # Dont penalize intercept
    E <- diag(c(E_c, E_b))
  }
  E
}

# Return a correlation matrix for the desired structure: Ind,
#  Compound Symmetry, or AR(1)
corrmat <- function(alpha = 0, m, type = "Ind") {
  if (type == "Ind") {
    R <- diag(m)
  } else if (type == "CS") {
    one <- rep(1, times = m)
    R <- (1 - alpha) * diag(m) + alpha * tcrossprod(one)
  } else if (type == "AR1") {
    times <- c(1:m)
    H <- abs(outer(times, times, "-"))
    R <- alpha^H
    R[cbind(1:m, 1:m)] <- R[cbind(1:m, 1:m)]
    return(R)
  } else {
    err <- c("Error in corrmat(): ",
             "Choose type either \"Ind\", \"CS\" or \"AR1\"")
    stop(paste(err, sep = ""))
  }
}

# Transform raw design matrices X and Z to form compatible for mixed analysis
# Returns a 2N x (p + q) matrix, where N = dim(X)[1] = dim(Z)[1],
#  p = dim(X)[2], and q = dim(Z)[2]
#              __           _
# x_i, z_i -> | x_i^T   0    |
#             |_  0   z_i^T _|
trans_X_mixed <- function(X, Z) {
  cbind(apply(X, 2, interweave_zeros, FALSE),
        apply(Z, 2, interweave_zeros, TRUE))
}

# Insert zeros between elements of a vector
# Used in trans_x_mixed()
interweave_zeros <- function(v, before) {
  n <- length(v)
  if (before) {
    c(rbind(0, v))
  } else {
    c(rbind(v, 0))
  }
}

# return standardized matrix with column means and sds.
stdMat <- function(X, intercept = TRUE) {
  if (intercept) {
    temp <- X[, 2:dim(X)[2]]
    Xstd <- cbind(1, scale(temp))
    cmean <- attr(scale(temp), "scaled:center")
    csd <- attr(scale(temp), "scaled:scale")
  } else {
    Xstd <- scale(X)
    cmean <- attr(Xstd, "scaled:center")
    csd <- attr(Xstd, "scaled:scale")
  }
  list(Xstd = Xstd, cmean = cmean, csd = csd)
}

# Unstandardize standardized regression coefficients. Assumes an intercept
unstandardize_coefs <- function(Beta.new, xmean, xsd, zmean, zsd) {
  if (is.null(zmean)) {
    # case for family != 'Mixed'
    a <- Beta.new[1]
    b <- Beta.new[2:length(Beta.new)]
    b <- b/xsd
    a <- a - sum(b * xmean)
    Beta.unstd <- c(a, b)
  } else {
    # case for family == 'Mixed'
    # Assumes intercepts are the first and (p+1)th elements
    a.cont <- Beta.new[1]
    b.cont <- Beta.new[2:(length(xmean) + 1)]
    a.bin <- Beta.new[length(xmean) + 2]
    b.bin <- Beta.new[(length(xmean) + 3):length(Beta.new)]
    b.cont <- b.cont/xsd
    b.bin <- b.bin/zsd
    a.cont <- a.cont - sum(b.cont * xmean)
    a.bin <- a.bin - sum(b.bin * zmean)
    Beta.unstd <- c(a.cont, b.cont, a.bin, b.bin)
  }
  Beta.unstd
}

# scale function for weighted design matrix
weighted.scale <- function(X, w = rep(1, nrow(X))) {
  if (identical(X[, 1], rep(1, nrow(X)))) {
    # intercept
    has_intercept <- TRUE
    X <- X[, -1]
  } else has_intercept <- FALSE
  m <- apply(X, 2, stats::weighted.mean, w)
  s <- apply(X, 2, weighted.sd, w)
  X <- scale(X, center = m, scale = s)
  if (has_intercept)
    X <- cbind(1, X)
  list(Xstd = X, cmean = m, csd = s)
}

# weighted standard deviation.
# In the trivial case, denominator is n, not n-1
weighted.sd <- function(x, w = rep(1, length(x))) {
  stopifnot(length(x) == length(w))
  s <- sqrt(sum((x * w - stats::weighted.mean(x, w))^2)/sum(w))
}

# check that intercept in design matrix is first column only, if exists
check_intercept <- function(X) {
  if(is.null(X)) {
    return(TRUE)
  } else {
    p <- dim(X)[2]
    for (j in 2:p) {
      if(identical(X[, j], rep(1, dim(X)[1]))) {
        return(FALSE)
      }
    }
    return(TRUE)
  }
}

# # Cleanup C++ code after package unloaded
# .onUnload <- function (libpath) {
#   library.dynam.unload("pgee", libpath)
# }
