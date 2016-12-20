#' pgee.mixed: Penalized Generalized Estimating Equations for Bivariate Mixed Outcomes
#'
#' Perform simultaneous estimation and variable selection for correlated bivariate
#' mixed outcomes (one continuous outcome and one binary outcome per cluster) using
#' penalized generalized estimating equations. In addition, clustered Gaussian and binary
#' outcomes can also be modeled. The SCAD, MCP, and LASSO penalties are supported.
#' Cross-validation can be performed to find the optimal regularization parameter(s).
#'
#' @references Deshpande, V., Dey, D. K., and Schifano, E. D. (2016). Variable
#' selection for correlated bivariate mixed outcomes using penalized generalized
#' estimating equations. Technical Report 16-23, Department of Statistics,
#' University of Connecticut, Storrs, CT.
#' @references Wang, L., Zhou, J., and Qu, A. (2012). Penalized generalized
#'   estimating equations for high-dimensional longitudinal data analysis.
#'   Biometrics, \strong{68}, 353–360.
#'
#' @docType package
#' @name pgee.mixed
NULL

#' @useDynLib pgee.mixed
#' @importFrom Rcpp sourceCpp
NULL

# return the sequence of coefficients to sum probabilities over for fdr
#   calculation
# b is the vector of stacked regression coefficients: b = c(b_c, b_b)
# last.cont should be the position of the last continuous regression coefficient
#   (Note that this assumes that all binary coefs follow the cont coefs)
# type: See FDR.mixed()
get_fdrseq <- function(b, last.cont, type, intercept = TRUE) {
  if (type == "all") {
    fdrseq <- 1:length(b)
    if (intercept)
      fdrseq <- setdiff(fdrseq, c(1, (last.cont + 1)))
  } else if (type == "continuous") {
    fdrseq <- 1:last.cont
    if (intercept)
      fdrseq <- setdiff(fdrseq, 1)
  } else if (type == "binary") {
    fdrseq <- (last.cont + 1):length(b)
    if (intercept)
      fdrseq <- setdiff(fdrseq, (last.cont + 1))
  } else stop(paste("Error in get_fdrseq(): argument \"type\" should be either",
                    "\"all\", \"continuous\", or \"binary\")", sep = ""))
  fdrseq
}

# Estimate false discovery rate for mixed outcomes
# b = c(b_c, b_b) are the estimated regression coefficients
# y is the vector of mixed outcomes, given as(y_c1, y_b1, ..., y_cN, y_bN)
# X is the transformed design matrix for mixed outcome analysis
#   (via trans_X_mixed())
# w are observation weights
# lambda = c(lambda_c, lambda_b)
#   are the continous and binary tuning parameters
# alpha is the estimated association parameter between the outcomes
# last.cont: See get_fdrseq(). The default value assumes that half of the b's
#   correspond to the continuous outcomes
# type: "all": all coefficients other than intercepts
#       "continuous": only continuous regression coefs
#       "binary": only binary regression coefs
FDR.mixed <- function(b, y, X, w = rep(1, length(y)/2), lambda, alpha,
                      last.cont = NULL, type, intercept = TRUE) {
  m <- 2
  N <- length(y)/m
  id.c <- seq(1, m * N, 2)
  id.b <- seq(2, m * N, 2)
  # Calculate residuals
  mu <- meanfn(X %*% b, family = "Mixed")
  r <- y - mu
  r.c <- r[id.c]
  r.b <- r[id.b]
  # Account for weights
  r.c <- r.c * sqrt(w)
  r.b <- r.b * sqrt(w)
  # Estimate variances
  sig2.c <- crossprod(r.c)/sum(w)
  sig2.b <- crossprod(r.b)/sum(w)
  # Calculate W matrix
  v <- varfn(mu, family = "Mixed")
  Rhat <- matrix(c(1, alpha, alpha, 1), nrow = m)
  W <- CppW2(X, v, Rhat, N, w)

  FDR <- 0
  if (is.null(last.cont)) last.cont <- length(b)/2
  if (last.cont > length(b))
    stop("Error in FDR.mixed(): last.cont cannot be greater than length(b)")
  # Compute block diagonal matrix of covariances
  V <- kronecker(diag(N), matrix(c(sig2.c,
                                   alpha * sqrt(sig2.c * sig2.b),
                                   alpha * sqrt(sig2.c * sig2.b),
                                   sig2.b),
                                 nrow = m))
  # Compute probabilities of false positives, and sum to get FDR
  fdrseq <- get_fdrseq(b, last.cont, type, intercept)
  for (k in fdrseq) {
    if (k <= last.cont) {
      l <- lambda[1]
    } else l <- lambda[2]
    pr.fd.k <- 2 * stats::pnorm(-l * sum(w)/sqrt(crossprod(W[, k], V) %*% W[, k]))
    FDR <- FDR + pr.fd.k
  }
  nz <- sum(abs(b[fdrseq]) > 0.001)
  if (nz == 0) return(0)  # Edge case, all coefficients are zero
  as.vector(FDR/nz)
}


# Estimate Pearson residuals
# y: Vector of outcomes. For family=="Mixed", should be of form
#    c(y_c1, y_b1,...,y_cN, y_bN).
# X: Design matrix. For family=="Mixed", should be transformed
#    via trans_X_mixed().
# family: "Gaussian", "Binomial", or "Mixed"
presid.est <- function(y, X, Beta, family="Gaussian") {
  mu <- meanfn(X %*% Beta, family)
  r <- (c(y) - c(mu)) / c(sqrt(varfn(mu, family)))
}

# Estimate dispersion parameter (continuous outcomes).
# Note that our phi is the reciprocal of Liang and Zeger (1986)'s phi,
#   as we multipy by phi whereas they divide by it in the estimating eqn
# pres: Pearson residuals
# p: number of regression coefficients
# m: cluster size
# w: sampling weights.
phi.est <- function(pres, p, w) {
  phi <- sum(w * (pres^2))/(sum(w))
}

# Weighted biserial correlation
# x: continuous variable
# y: binary (0/1) variable
# w: weights
weighted.biserial <- function(x, y, w) {
  x_1 <- x[y == 1]
  x_0 <- x[y == 0]
  w_1 <- w[y == 1]
  w_0 <- w[y == 0]
  # don't use n-1 because w's might add to 1 in denominator
  s_w_x <- sqrt(sum(w * (x - stats::weighted.mean(x, w))^2)/(sum(w)))
  r <- (stats::weighted.mean(x_1, w_1) - stats::weighted.mean(x_0, w_0)) *
        stats::weighted.mean(y, w) * (1 - stats::weighted.mean(y, w)) /
        stats::dnorm(stats::qnorm(stats::weighted.mean(y, w))) / s_w_x
}

# Estimate the working correlation matrix R - nonmixed case
# pres: Estimated Pearson residuals
# w: vector of cluster weights
# m: cluster size, assumed equal for all clusters
# p: Number of regression coefficients
# phi: dispersion
# type: Working correlation structure: Ind, CS or AR1
# TO DO: Incorporate weights
alpha.est <- function(pres, w, m, p, phi, type = "Ind") {
  if (type == "Ind") {
    alpha_hat <- 0
  } else if (type == "CS") {
    alpha_hat <- CppAlphaCS(pres, w, m, p, phi)
  } else if (type == "AR1") {
    alpha_hat <- CppAlphaAR1(pres, w, m, p, phi)
  } else {
    stop("Error in R.alpha.est():
         type should be either \"Ind\", \"CS\", or \"AR1\".")
  }
  # Bound alpha
  if (alpha_hat >= 1) {
    warning("Association parameter hit upper bound")
    alpha_hat <- 0.99
  }
  if (alpha_hat <= -1) {
    warning("Association parameter hit lower bound")
    alpha_hat <- -0.99
  }
  alpha_hat
  }

# Estimate working correlation matrix - mixed outcomes case
# biserial correlation b/w continuous residuals and binary outcomes
# y: mixed outcomes vector (y_c1, y_b1,...,y_cN, y_bN)
# eta: linear predictor vector corresponding to y
# w: weights
alpha.est.mixed <- function(y, eta, w) {
  # Separate out continuous and binary
  id <- 1:length(y)
  id_b <- (id%%2 == 0)
  id_c <- (id%%2 != 0)
  y_c <- y[id_c]
  y_b <- y[id_b]
  eta_c <- eta[id_c]
  eta_b <- eta[id_b]
  r_c <- y_c - meanfn(eta_c, "Gaussian")
  # Biserial correlation between cont residual and y_b
  alpha <- c(weighted.biserial(r_c, y_b, w))
  if (abs(alpha) > 0.8) {
    warning("Absolute value of alpha estimated larger than 0.8.
            This is almost impossible for mixed responses.")
  }
  alpha
}

#' Penalized Generalized Estimating Equations
#'
#' Estimate regression coefficients using Penalized Generalized Estimating
#' Equations (PGEEs). Linear and binary logistic models are currently supported.
#' In particular, can handle the case of bivariate correlated mixed outcomes, in
#' which each cluster consists of one continuous outcome and one binary outcome.
#'
#' \code{pgee.fit} estimates the regression coefficients for a single value of the
#' tuning paramter (or a single pair of tuning parameters in the mixed outcomes
#' case). To select optimal tuning parameter(s) via k-fold cross validation, see
#' \code{cv.pgee}.
#'
#' For bivariate mixed outcomes, the false discovery rate can be estimated.
#'
#' @param N Number of clusters.
#' @param m Cluster size. Assumed equal across all clusters. Should be set to 2
#'   for family=="Mixed".
#' @param X Design matrix. If family=="Mixed", then design matrix for continuous
#'   responses. For family!="Mixed", should have N*m rows. For family=="Mixed",
#'   should have N rows. For standardize=TRUE, the first column should be a
#'   column vector of ones, corresponding to the intercept.
#' @param Z Design matrix for binary responses for family=="Mixed". Should not
#'   be provided for other family types. If not provided for family=="Mixed", is
#'   set equal to X. For family!="Mixed", should have N*m rows. For
#'   family=="Mixed", should have N rows. For standardize=TRUE, the first
#'   column should be a column vector of ones, corresponding to the intercept.
#' @param y Response vector. Don't use this argument for family == "Mixed".
#'   Instead, use arguments yc and yb. Since the cluster size is assumed equal
#'   across clusters, the vector is assumed to have the form c(y_1,
#'   y_2,...,y_N), with y_i = c(y_i1,...,y_im).
#' @param yc Continuous response vector. Use only for family=="Mixed".
#' @param yb Binary (0/1) response vector. Use only for family=="Mixed".
#' @param wctype Working correlation type; one of "Ind", "CS", or "AR1". For
#'   family=="Mixed", "CS" and "AR1" are equivalent.
#' @param family "Gaussian", "Binomial", or "Mixed" (use the last for bivariate
#'   mixed outcomes). Note that for "Binomial", currently only binary outcomes
#'   are supported.
#' @param lambda Tuning parameter(s). A vector of two tuning parameters should
#'   be provided for family=="Mixed" (one for the continuous outcome
#'   coefficients, and one of the binary outcome coefficients). Otherwise, a
#'   single tuning parameter should be provided.
#' @param eps Disturbance in the Linear Quadratic Approximation algorithm.
#' @param maxiter The maximum number of iterations the Newton algorithm tries
#'   before declaring failure to converge.
#' @param tol.coef Converge of the Newton algorithm is declared if two
#'   conditions are met: The L1-norm of the difference of successive iterates
#'   should be less than tol.coef AND the L1-norm of the penalized score should
#'   be less than tol.score.
#' @param tol.score See \code{tol.coef}.
#' @param init Vector of initial values for regression coefficients. For
#'   family=="Mixed", should be c(init_c, init_b). Defaults to glm values.
#' @param standardize Standardize the design matrices prior to estimation?
#' @param penalty "SCAD", "MCP", or "LASSO".
#' @param weights Vector of cluster weights. All observations in a cluster are
#'   assumed to have the same weight.
#' @param FDR Should the false discovery rate be estimated for family=="Mixed"?
#'   Currently, FDR cannot be estimated for other family types.
#' @param fdr.corr Association parameter to use in FDR estimation. The default
#'   is to use the association parameter estimated from the PGEEs.
#' @param fdr.type Estimate the FDR for only the coefficients corresponding to
#'   the continuous outcomes ("continuous"), for only the coefficients
#'   corresponding to the binary outcomes ("binary"), or for all coefficients
#'   ("all", the default).
#' @return A list
#'   \item{coefficients}{Vector of estimated regression
#'   coefficients. For family=="Mixed", this takes the form c(coef_c, coef_b).}
#'   \item{vcov}{Sandwich formula based covariance matrix of estimated
#'   regression coefficients (other than the intercept(s)). The row/column
#'   names correspond to elements of \code{coefficients}.}
#'   \item{phi}{Estimated dispersion parameter.}
#'   \item{alpha}{Estimated association parameter.}
#'   \item{num_iterations}{Number of iterations the Newton algorithm ran.}
#'   \item{converge}{0=converged, 1=did not converge.}
#'   \item{PenScore}{Vector of penalized score functions at the estimated
#'   regression coefficients. If the algorithm converges, then these should be
#'   close to zero.}
#'   \item{FDR}{Estimated FDR for family=="Mixed", if requested.}
#' @examples
#' set.seed(100)
#' # Gaussian
#' N <- 100
#' m <- 10
#' p <- 10
#' y <- rnorm(N * m)
#' # If you want standardize = TRUE, you must provide an intercept.
#' X <- cbind(1, matrix(rnorm(N * m * (p - 1)), N * m, p - 1))
#' fit <- pgee.fit(X = X, y = y, N = N, m = m, lambda = 0.5, wctype = "CS",
#'             family = "Gaussian")
#' str(fit)
#' fit$coefficients
#' fit$vcov
#'
#' # Binary
#' y <- sample(0:1, N*m, replace = TRUE)
#' fit <- pgee.fit(X = X, y = y, N = N, m = m, lambda = 0.1, wctype = "CS",
#'             family = "Binomial")
#' str(fit)
#' fit$coefficients
#' fit$vcov
#'
#' # Bivariate mixed outcomes
#' # Generate some data
#' Bc <- c(2.0, 3.0, 1.5, 2.0, rep(0, times = p - 4))
#' Bb <- c(0.7, -0.7, -0.4, rep(0, times = p - 3))
#' dat <- gen_mixed_data(Bc, Bc, N, 0.5)
#' # Estimate regression coefficients and false discovery rate
#' fit <- pgee.fit(X = dat$X, yc = dat$yc, yb = dat$yb, N = N, m = 2,
#'             wctype = "CS", family = "Mixed", lambda = c(0.1, 0.05),
#'             FDR = TRUE)
#' str(fit)
#' fit$coefficients
#' fit$vcov
#' @export pgee.fit
pgee.fit <- function(N, m, X, Z = NULL, y = NULL, yc = NULL, yb = NULL,
                     wctype = "Ind", family = "Gaussian", lambda = 0,
                     eps = 1e-06, maxiter = 1000, tol.coef = 1e-03,
                     tol.score = 1e-03, init = NULL,
                     standardize = TRUE, penalty = "SCAD", weights = rep(1, N),
                     FDR = FALSE, fdr.corr = NULL, fdr.type = "all") {

  # Checks and basic setup ----------------------------------------------------
  valid_wc_types <- c("Ind", "CS", "AR1")
  valid_family_types <- c("Gaussian", "Binomial", "Mixed")
  if (!(wctype %in% valid_wc_types))
    stop("Error in pgee.fit():
         Invalid Working Correlation Matrix type specified")
  if (!(family %in% valid_family_types))
    stop("Error in pgee.fit(): Invalid family specified")
  if (family == "Mixed") {
    if (m != 2 | length(lambda) != 2) {
      stop("Error in pgee.fit(): Check m, lambda")
    }
  }

  # y checks
  if (family == "Mixed") {
    if (any(is.null(yc), is.null(yb))) {
      stop("Error in pgee.fit(): For family==\"Mixed\", use arguments yc and yb
            as response vectors.")
    }
    if (!is.null(y)) {
      warning("For family==\"Mixed\", argument y is ignored.
               Arguments yc and yb are used as response vectors.")
    }
    y <- c(rbind(yc, yb))
  } else {
    if (is.null(y)) {
      stop("Error in pgee.fit(): Please provide response vector using argument y.")
    }
    if (any(!is.null(yc), !is.null(yb))) {
      warning("For family \"Gaussian\" or \"Binomial\", arguments yc and yb
              are ignored. Argument y is used as the response vector.")
    }
  }

  # Weights, normalize
  if (methods::hasArg(weights)) {
    weights <- weights/sum(weights) * N
  }

  # Force intercept to be first column
  if (!check_intercept(X))
    stop("Error in pgee.fit(): In design matrix X, specify the column of ones
          for the intercept in the first column only.")
  if (!check_intercept(Z))
    stop("Error in pgee.fit(): In design matrix Z, specify the column of ones
          for the intercept in the first column only.")

  # For mixed, if Z isn't provided, use X
  if (is.null(Z) && family == "Mixed") Z <- X

  # Can't standardize if no intercept column
  if (standardize) {
    if (!identical(X[, 1], rep(1, dim(X)[1]))) {
      message("Can't standardize design matrix X without
              an intercept column. Adding it now.")
      X <- cbind(1, X)
    }
    if (!is.null(Z)) {
      if (!identical(Z[, 1], rep(1, dim(Z)[1]))) {
        message("Can't standardize design matrix Z without
            an intercept column. Adding it now.")
        Z <- cbind(1, Z)
      }
    }
  }

  # intercept flag, useful for Emat non-penalization of intercept
  intercept <- identical(X[, 1], rep(1, dim(X)[1]))
  # If X has an intercept, make sure Z has one
  if (intercept && !standardize && !is.null(Z) &&
      !identical(Z[, 1], rep(1, dim(Z)[1]))) {
    stop("In pgee.fit(): For family==\"Mixed\", current implementation does not
          support intercept for X and not for Z. Please add an intercept column
          to Z.")
  }
  # ---------------------------------------------------------------------------
  # Additional pre-algorithm setup --------------------------------------------
  # Standardize design matrices if required
  p.x <- dim(X)[2]
  if (standardize) {
    X.ustd <- X
    Z.ustd <- Z
    # l <- stdMat(X)
    # X <- l$Xstd
    # xmean <- l$cmean
    # xsd <- l$csd
    # try with weighted mean/sd
    if (family == "Mixed") {
      l <- weighted.scale(X, weights)
    } else {
      l <- weighted.scale(X, rep(weights, each = m))
    }

    X <- l$Xstd
    xmean <- l$cmean
    xsd <- l$csd
  }
  if (!is.null(Z)) {
    p.z <- dim(Z)[2]
    if (standardize) {
      # l <- stdMat(Z)
      # Z <- l$Xstd
      # zmean <- l$cmean
      # zsd <- l$csd
      # try with weighted mean/sd
      l <- weighted.scale(Z, weights)
      Z <- l$Xstd
      zmean <- l$cmean
      zsd <- l$csd
    }
  } else zmean <- zsd <- NULL

  # transform X for mixed
  if (family == "Mixed") {
    Xold <- X
    X <- trans_X_mixed(X, Z)
  }

  # Inital values
  # wctype 'Ind'
  Rhat <- diag(m)
  alpha_hat <- 0
  # Currently estimating alpha only once, with glm estimates
  if (family == "Mixed") {
    # For now, compute correlation once from initial values
    if (is.null(init)) {
      Bn0 <- stats::glm.fit(Xold, yc, weights = weights)$coefficients
      Bb0 <- stats::glm.fit(Z, yb, family = stats::binomial(link = "logit"),
                     weights = weights,
                     start = rep(0, dim(Z)[2]))$coefficients
      Beta_init <- c(Bn0, Bb0)
      rm(Bn0, Bb0)
    } else {
      Beta_init <- init
      rm(init)
    }
    if (wctype != "Ind") {
      # alpha estimation doesn't use dispersion, for now
      alpha_hat <- alpha.est.mixed(y, X %*% Beta_init, weights)
    }
    Rhat <- corrmat(alpha_hat, m, "CS")
  } else if (family == "Gaussian") {
    if (is.null(init)) {
      Beta_init <- stats::glm.fit(X, y,
                           weights = rep(weights, each = m))$coefficients
    } else {
      Beta_init <- init
      rm(init)
    }
    pres <- presid.est(y, X, Beta_init, "Gaussian")
    phi_hat <- phi.est(pres, p.x, rep(weights, each = m))
    if (wctype != "Ind") {
      alpha_hat <- alpha.est(pres, weights, m, p.x, phi_hat, wctype)
    }
    Rhat <- corrmat(alpha_hat, m, wctype)
  } else {  # family == "Binomial"
    if (is.null(init)) {
      Beta_init <- stats::glm.fit(X, y, family = stats::binomial(link = "logit"),
                           weights = rep(weights, each = m),
                           start = rep(0, dim(X)[2]))$coefficients
    } else {
      Beta_init <- init
      rm(init)
    }
    pres <- presid.est(y, X, Beta_init, "Binomial")
    phi_hat <- phi.est(pres, p.x, rep(weights, each = m))
    if (wctype != "Ind") {
      alpha_hat <- alpha.est(pres, weights, m, length(Beta_init), phi_hat, wctype)
    }
    Rhat <- corrmat(alpha_hat, m, wctype)
  }
  # ---------------------------------------------------------------------------
  # Newton iterative algorithm ------------------------------------------------
  # NR Hessian fraction, use successively smaller if diverges
  delta <- c(1, 0.5, 0.1, 0.05)
  for (id in seq_along(delta)) {
    # Reset for new run
    Beta.new <- Beta_init
    itercount <- 0
    converge <- 0  # 0=converged, 1=did not converge
    condition <- TRUE  # iterates until condition=FALSE
    while (condition) {
      Beta.old <- Beta.new

      # compute continuous dispersion
      if (family == "Mixed") {
        phi_hat <- phi.est(yc - Xold %*% Beta.old[1:p.x], p.x, weights)
        # bin disp assumed 1
      } else if (family == "Gaussian") {
        phi_hat <- phi.est(y - X %*% Beta.old, p.x, rep(weights, each = m))
      } else phi_hat <- 1  # Assuming dispersion of 1 for binary

      # Now update Beta
      mu <- meanfn(X %*% Beta.old, family = family)
      v <- varfn(mu, family = family)
      E <- Emat_wrap(Beta.old, lambda = lambda, family = family, eps = eps,
                     intercept = intercept, penalty)
      if (family == "Mixed") {
        H <- CppHess2(X, v, alpha_hat, phi_hat, N, weights)
        S <- CppScore2(X, v, alpha_hat, y, mu, N, phi_hat, weights)
      } else {
        H <- CppHess(X, v, Rhat, phi_hat, N, weights)
        S <- CppScore(X, v, Rhat, y, mu, N, phi_hat, weights)
      }
      Beta.new <- Beta.old + delta[id] * as.vector(CppIncBeta(Beta.old, S,
                                                              H, E, N))

      # update Rhat - this works badly if lambdas are not good, so just using
      #   initial alpha from glm estimates
      # alpha_hat <- R.alpha.est.mixed(y, X %*% Beta.new)
      # Rhat <- corrmat(alpha_hat, m, 'CS')

      itercount <- itercount + 1
      if (itercount > maxiter) {
        warning('Did not coverge before maxiter number of iterations')
        converge <- 1
        break
      }
      penscore <- as.vector(PenScore(Beta.old, S, E, N))
      condition <- (sum(abs(as.vector(Beta.old - Beta.new))) > tol.coef &
                    sum(abs(penscore)) > tol.score)

      # if diverges, go to smaller delta
      if (is.na(condition))
        break
    }
    if (identical(condition, FALSE))
      break
  }
  # if the smallest delta still didnt work, warning of non-convergence
  if (id == length(delta) & ((is.na(condition)) | condition)) {
    warning("In pgee.fit(): Beta diverges for smallest delta")
    converge <- 1
  }
  # ---------------------------------------------------------------------------
  # Post-algorithm ------------------------------------------------------------

  # Approaches to computing covariance matrix:
  #   1. Only non-zero: This is turning out hard to implement. Work on it later
  #   2. Everything, and get rid of intercept. Easier, doing this now.
  #   Approach 1 --------------------------------------------------------------
  # # Explicitly threshold
  # nzid <- abs(Beta.new) > 1e-03
  # Beta.new[!nzid] <- 0
  #
  # # Sandwich covariance matrix
  # # Don't care about intercept
  # if (intercept) {
  #   nzid[1] <- FALSE
  #   if (family == "Mixed") nzid[p.x + 1] <- FALSE
  # }
  # se <- sandwich(y, X[, nzid], Beta.new[nzid], alpha_hat, lambda, N, m,
  #                family, wctype, eps, FALSE, phi_hat, penalty, weights)
  # stopifnot(dim(se)[1] == dim(se)[2])
  # # unstandardize, if necessary
  # if (standardize) {
  #   if (family != "Mixed") {
  #     nzid2 <- nzid[-1]
  #     stopifnot(dim(se)[2] == length(xsd[nzid2]))
  #     se <- se / outer(xsd[nzid2], xsd[nzid2])
  #   } else {
  #     nzid2 <- nzid[c(-1, -(p.x + 1))]
  #     xzsd <- c(xsd, zsd)
  #     stopifnot(dim(se)[2] == length(xzsd[nzid2]))
  #     se <- se / outer(xzsd[nzid2], xzsd[nzid2])
  #   }
  # }
  # stopifnot(length(which(nzid)) == dim(se)[2])
  # rownames(se) <- colnames(se) <- which(nzid)
  # -----------------------------------------------------------------------------
  # Approach 2
  se <- sandwich(y, X, Beta.new, alpha_hat, lambda, N, m,
                 family, wctype, eps, intercept, phi_hat, penalty, weights)
  # Extract the non-intercept submatrix
  if(intercept) {
    if (family != "Mixed") {
      se <- se[-1, -1]
      se.ids <- 2:p.x
    } else {
      se <- se[c(-1, -(p.x+1)), c(-1, -(p.x+1))]
      se.ids <- c(2:p.x, (p.x+2):(p.x+p.z))
    }
  } else {
    if (family != "Mixed") {
      se.ids <- 1:p.x
    } else {
      se.ids <- 1:(p.x+p.z)
    }
  }
  # Unstandardize, if required
  if (standardize) {
    if (family != "Mixed") {
      se <- se / outer(xsd, xsd)
    } else {
      xzsd <- c(xsd, zsd)
      se <- se / outer(xzsd, xzsd)
    }
  }
  rownames(se) <- colnames(se) <- se.ids
  # FDR if required
  if (FDR) {
    # Only allowed for mixed response with intercept
    if (family != "Mixed") {
      print("FDR can be estimated only for family==\"Mixed\".")
    } else {
      if (any(is.nan(Beta.new))) {
        print("Warning: Can't compute FDR, NaNs in estimated coefficients")
      } else {
        if (is.null(fdr.corr)) fdr.corr <- alpha_hat
        FDRest <- FDR.mixed(Beta.new, y, X, weights, lambda, fdr.corr,
                            last.cont = , type = fdr.type)
      }
    }
  } else FDRest <- NA

  # unstandardize coefficients, if necessary.
  Beta.std <- Beta.new
  if (standardize) {
    Beta.new <- unstandardize_coefs(Beta.new, xmean, xsd, zmean, zsd)
  }

  # return results as list
  l <- list(coefficients = Beta.new, vcov = se, phi = phi_hat, alpha = alpha_hat,
            num_iterations = itercount, converge = converge,
            PenScore = penscore, FDR = FDRest)
}
