# Return a list, each element being the cluster ids for the test fold.
createFolds <- function(y, k = 10, list = TRUE, returnTrain = FALSE) {
  if (is.numeric(y)) {
    cuts <- floor(length(y)/k)
    if (cuts < 2)
      cuts <- 2
    if (cuts > 5)
      cuts <- 5
    y <- cut(y, unique(stats::quantile(y, probs = seq(0, 1, length = cuts))),
             include.lowest = TRUE)
  }
  if (k < length(y)) {
    y <- factor(as.character(y))
    numInClass <- table(y)
    foldVector <- vector(mode = "integer", length(y))
    for (i in 1:length(numInClass)) {
      seqVector <- rep(1:k, numInClass[i]%/%k)
      if (numInClass[i]%%k > 0)
        seqVector <- c(seqVector, sample(1:k, numInClass[i]%%k))
      foldVector[which(y == dimnames(numInClass)$y[i])] <- sample(seqVector)
    }
  } else foldVector <- seq(along = y)
  if (list) {
    out <- split(seq(along = y), foldVector)
    names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))),
                        sep = "")
    if (returnTrain)
      out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
  } else out <- foldVector
  out
}

# loss functions: http://web.mit.edu/lrosasco/www/publications/loss.pdf loss
# loss functions for binary classification
lossfn_clss <- function(y, eta, type, weights) {
  # 0/1 coding
  switch(type, deviance = stats::weighted.mean(log(1 + exp(-2 * eta * (2 * y - 1))),
                                        weights)/log(2),
               classification = stats::weighted.mean(1 * (eta * (2 * y - 1) < 0),
                                              weights))
}

# loss functions for linear regression.
lossfn_reg <- function(y, eta, type, weights) {
  t = y - eta
  switch(type, square = stats::weighted.mean(t^2, weights),
               absolute = stats::weighted.mean(abs(t), weights))
}

# Loss function wrapper for nonmixed outcomes.
# Currently not penalizing for phi.
lossfn <- function(y, X, Beta, phi, family, weights, type) {
  if (family == "Gaussian") {
    # -sum(dnorm(y, mean=X%*%Beta, sd=phi, log=TRUE))
    lossfn_reg(y, X %*% Beta, type, weights)
  } else if (family == "Binomial") {
    # -sum(dbinom(y, 1, meanfn(X%*%Beta, 'Binomial'), log = TRUE))
    lossfn_clss(y, X %*% Beta, type, weights)
  } else stop("Error in lossfn(): Wrong family specified")
}

# loss function for mixed outcomes.
# y, X, Z, Beta: see pgee.fit().
lossfn_mixed <- function(y, X, Beta, type_c, type_b, marginal, Z, weights) {
  # Separate out y_c and y_b
  y_c <- y[seq(1, length(y), 2)]
  y_b <- y[seq(2, length(y), 2)]
  # Separate out Beta_c and Beta_b
  p.x <- dim(X)[2]
  p.z <- dim(Z)[2]
  Beta_c <- Beta[1:p.x]
  Beta_b <- Beta[(p.x + 1):(p.x + p.z)]
  eta_c <- X %*% Beta_c
  eta_b <- Z %*% Beta_b
  if (marginal == 0) {
    # mixed
    loss <- lossfn_reg(y_c, eta_c, type_c, weights) +
            lossfn_clss(y_b, eta_b, type_b, weights)
  } else if (marginal == 1) {
    # continuous
    loss <- lossfn_reg(y_c, eta_c, type_c, weights)
  } else if (marginal == 2) {
    # binary
    loss <- lossfn_clss(y_b, eta_b, type_b, weights)
  } else stop("Error in lossfn_mixed():
              Argyment marginal should be either 0, 1, or 2")
  loss
}

# Select optimal tuning parameter using K-fold CV, nonmixed outcomes.
# Folds are at the cluster level, not the observation level.
get_lambda_cv <- function(y, X, N, m, K, grid_lambda, Beta0,
                          wctype, family, eps, maxiter, tol.coef, tol.score,
                          warm, weights, penalty, type_c, type_b) {
  flds <- createFolds(c(1:N), K)
  LossMat <- matrix(, nrow = length(grid_lambda), ncol = K)
  for (j in 1:K) {
    clustid <- rep(1:N, each = m)
    # Get training and test sets
    id.test <- clustid %in% unlist(flds[j])
    id.train <- clustid %in% unlist(flds[-j])
    y.test <- y[id.test]
    y.train <- y[id.train]
    X.test <- X[id.test, ]
    X.train <- X[id.train, ]
    N.test <- length(unlist(flds[j]))
    N.train <- length(unlist(flds[-j]))
    w.train <- weights[unlist(flds[-j])]
    w.test <- weights[unlist(flds[j])]

    for (i in seq_along(grid_lambda)) {
      # Estimate Beta's from training set
      res.cv <- pgee.fit(X = X.train, y = y.train, N = N.train, m = m,
                     wctype = wctype, family = family, lambda = grid_lambda[i],
                     init = Beta0, penalty = penalty, weights = w.train,
                     maxiter = maxiter, tol.coef = tol.coef,
                     tol.score = tol.score, eps = eps)
      Beta.cv <- res.cv$coefficients
      # For warm start, use previous value of Beta as next initial value
      if (warm) {
        Beta0 <- Beta.cv
      }
      phi.cv <- res.cv$phi
      # Get loss using training set and estimates parameters
      LossMat[i, j] <- ifelse(any(is.na(Beta.cv)), 1e+100,
                              lossfn(y.test, X.test, Beta.cv, phi.cv,
                                     family = family,
                                     weights = rep(w.test, each = m),
                                     type = ifelse(family == "Gaussian",
                                                   type_c, type_b)))
    }
  }
  avgloss <- rowMeans(LossMat)
  lambda <- grid_lambda[which(rowMeans(LossMat) == min(rowMeans(LossMat)))]
  lambda <- lambda[1]  #break ties, return a unique lambda
  l <- list(lambda.opt = lambda,
            lambda.loss = avgloss[which(grid_lambda == lambda)],
            lambda.grid = grid_lambda,
            LossMat = LossMat)
}

# Select optimal tuning parameter using K-fold CV, mixed outcomes.
# Folds are at the cluster level, not the observation level.
get_lambda_cv_mixed <- function(y, X, N, K, grid1, grid2, init, wctype,
                                eps, maxiter, tol.coef, tol.score,
                                warm, type_c, type_b, marginal, Z,
                                weights, penalty) {
  flds <- createFolds(c(1:N), K)
  # The Loss matrix for single has folds in columns and grids in row.
  # To modify for mixed outcomes, we shall make each row correspond to a
  # unique combination of the two grids.
  # We shall let the first grid vary slower than the second grid. So the
  # rows shall correspond to something like:
  # (lam_11,...,lam_1(g2),......,lam_(g1)1,...,lam_(g1)(g2))
  g1 <- length(grid1)
  g2 <- length(grid2)
  LossMat <- matrix(, nrow = g1 * g2, ncol = K)

  # store cold start value
  init_cold <- init

  for (j in 1:K) {
    print(paste("fold", j))
    if (K > 1) {
      clustid <- rep(1:N, each = 2)  # m = 2
      # Get training and test sets
      id.test <- clustid %in% unlist(flds[j])
      id.train <- clustid %in% unlist(flds[-j])
      y.test <- y[id.test]
      y.train <- y[id.train]
      # Hack to get it to work for X, which is N x p instead of old Nm x p
      id.test.X <- unlist(flds[j])
      id.test.X <- id.test.X[order(id.test.X)]
      id.train.X <- unlist(flds[-j])
      id.train.X <- id.train.X[order(id.train.X)]
      X.test <- X[id.test.X, ]
      X.train <- X[id.train.X, ]
      Z.test <- Z[id.test.X, ]
      Z.train <- Z[id.train.X, ]
      N.test <- length(unlist(flds[j]))
      N.train <- length(unlist(flds[-j]))
      # dont use id.train/id.test, they contain 2x repeated values
      w.train <- weights[unlist(flds[-j])]
      w.test <- weights[unlist(flds[j])]

    } else {
      # K=1 means use whole dataset, no train and test separate
      y.test <- y.train <- y
      X.train <- X.test <- X
      Z.train <- Z.test <- Z
      N.train <- N.test <- N
      w.train <- w.test <- weights
    }

    # Traverse in descending order and use warm starts (by default)
    for (i1 in order(grid1, decreasing = TRUE)) {
      for (i2 in order(grid2, decreasing = TRUE)) {

        print(paste("i1: ", i1, "   i2: ", i2, sep = ""))
        # Estimate Beta's from training set
        res.cv <- pgee.fit(X = X.train, Z = Z.train,
                       yc = y.train[seq(1, N.train*2, 2)],
                       yb = y.train[seq(2, N.train*2, 2)],
                       N = N.train, m = 2,
                       wctype = wctype, family = "Mixed",
                       lambda = c(grid1[i1], grid2[i2]), init = init,
                       penalty = penalty, weights = w.train, maxiter = maxiter,
                       tol.coef = tol.coef, tol.score = tol.score, eps = eps)
        Beta.cv <- res.cv$coefficients

        # Warm start: Use previous value until you "cross over",
        #   then revert to nearest. If can't converge, use cold start.
        if (warm) {
          if(res.cv$converge != 0) {
            init <- init_cold
          } else {
            # reverse traversal of grids
            if (i2 == g2)
              init_reset <- Beta.cv  # store this for the reset of i
            if (i2 == 1) {
              init <- init_reset  # update the value for reset of i
            } else {
              init <- Beta.cv  # usual update
            }
          }
        }

        # Get loss using training set and estimates parameters
        id <- g2 * (i1 - 1) + i2
        LossMat[id, j] <- ifelse(any(is.na(Beta.cv)), 1e+100,
                                 lossfn_mixed(y.test, X.test, Beta.cv,
                                              type_c = type_c, type_b = type_b,
                                              marginal = marginal, Z.test,
                                              weights = w.test))
      }
    }
  }
  avgloss <- rowMeans(LossMat)  # Averaging across folds
  minloss.id <- which(rowMeans(LossMat) == min(rowMeans(LossMat)))
  minloss.id <- minloss.id[1]  # to break ties
  # To solve minloss.id = g2(i1-1) + i2 for i1 and i2:
  i2 <- ifelse(minloss.id%%g2 == 0, g2, minloss.id%%g2)
  i1 <- (minloss.id - i2)/g2 + 1
  lambda1 <- grid1[i1]
  lambda2 <- grid2[i2]
  l <- list(lambda = c(lambda1, lambda2), lambda.loss = avgloss[minloss.id],
            LossMat = LossMat)
}

# Get optimal lambda via CV, then fit to entire data. Nonmixed outcomes.
cv.pgee.nonmixed <- function(X, y, N, m, K, wctype, family, grid_lambda,
                          eps, maxiter, tol.coef, tol.score, init,
                          penalty, warm, type_c, type_b, weights) {
  # checks
  stopifnot(family == "Gaussian" | family == "Binomial")
  stopifnot(length(y) == N*m)
  stopifnot(dim(X)[1] == N*m)
  stopifnot(all(K >= 1, K <= N))
  stopifnot(length(weights) == N)
  stopifnot(eps > 0)

  lambda.list <- get_lambda_cv(y, X, N, m, K, grid_lambda, init, wctype,
                               family, eps, maxiter, tol.coef, tol.score,
                               warm, weights, penalty, type_c, type_b)
  lambda.opt <- lambda.list$lambda.opt
  fit.results <- pgee.fit(X = X, y = y, N = N, m = m, wctype = wctype,
                      family = family, lambda = lambda.opt, eps = eps,
                      maxiter = maxiter, tol.coef = tol.coef,
                      tol.score = tol.score, init = init,
                      penalty = penalty, weights = weights)
  fit.results$lambda.loss <- lambda.list$lambda.loss
  fit.results$LossMat <- lambda.list$LossMat
  fit.results
}

# Get optimal lambda via CV, then fit to entire data. Mixed outcomes.
cv.pgee.mixed <- function(X, Z = X, yc, yb, N, K, wctype, grid1, grid2,
                          eps, maxiter, tol.coef, tol.score,
                          init, penalty, warm, type_c,
                          type_b, marginal, weights, FDR, fdr.corr, fdr.type) {
  # get both optimal lambdas
  lambda.list <- get_lambda_cv_mixed(c(rbind(yc, yb)), X, N, K, grid1, grid2,
                                     init, wctype, eps, maxiter, tol.coef,
                                     tol.score, warm, type_c,
                                     type_b, marginal, Z, weights, penalty)
  lambda <- lambda.list$lambda  # vector of length 2
  fit.results <- pgee.fit(X = X, Z = Z, yc = yc, yb = yb, N = N, m = 2,
                      wctype = wctype, family = "Mixed", lambda = lambda,
                      eps = eps, maxiter = maxiter, tol.coef = tol.coef,
                      tol.score = tol.score, init = init,
                      penalty = penalty, weights = weights, FDR = FDR,
                      fdr.corr = fdr.corr, fdr.type = fdr.type)
  # Add avgloss of optimal lambda to results
  fit.results$lambda.loss <- lambda.list$lambda.loss
  fit.results$LossMat <- lambda.list$LossMat
  fit.results
}

#' Cross validation for Penalized Generalized Estimating Equations
#'
#' Performs k-fold cross-validation for Penalized Generalized Estimating
#' Equations (PGEEs) over grid(s) of tuning parameters lambda. Linear and binary
#' logistic models are supported. In particular, can handle the case of
#' bivariate correlated mixed outcomes, in which each cluster consists of one
#' continuous outcome and one binary outcome.
#'
#' The function calls \code{pgee.fit} \code{K} times, each time leaving out
#' 1/\code{K} of the data. The cross-validation error is determined by the
#' arguments \code{type_c} and \code{type_b}. For family=="Mixed", the
#' cross-validation error is (by default) the sum of the continuous error and
#' the binary error.
#'
#' @param N Number of clusters.
#' @param m Cluster size. Assumed equal across all clusters. Should be set to 2
#'   for family=="Mixed".
#' @param X Design matrix. If family=="Mixed", then design matrix for continuous
#'   responses. For family!="Mixed", should have N*m rows. For family=="Mixed",
#'   should have N rows.
#' @param Z Design matrix for binary responses for family=="Mixed". Should not
#'   be provided for other family types. If not provided for family=="Mixed", is
#'   set equal to X. For family!="Mixed", should have N*m rows. For
#'   family=="Mixed", should have N rows.
#' @param y Response vector. Don't use this argument for family == "Mixed".
#'   Instead, use arguments yc and yb. Since the cluster size is assumed equal
#'   across clusters, the vector is assumed to have the form c(y_1,
#'   y_2,...,y_N), with y_i = c(y_i1,...,y_im).
#' @param yc Continuous response vector. Use only for family=="Mixed".
#' @param yb Binary (0/1) response vector. Use only for family=="Mixed".
#' @param K Number of folds.
#' @param grid1 For family!="Mixed", the grid of tuning parameters. For
#'   family=="Mixed", the grid of tuning parameters for coefficients
#'   corresponding to the continuous outcomes.
#' @param grid2 For family=="Mixed", the grid of tuning parameters for
#'   coefficients corresponding to the binary outcomes. Not used for
#'   family!="Mixed".
#' @param wctype Working correlation type; one of "Ind", "CS", or "AR1". For
#'   family=="Mixed", "CS" and "AR1" are equivalent.
#' @param family "Gaussian", "Binomial", or "Mixed" (use the last for bivariate
#'   mixed outcomes). Note that for "Binomial", currently only binary outcomes
#'   are supported.
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
#' @param warm Use warm starts?
#' @param weights Vector of cluster weights. All observations in a cluster are
#'   assumed to have the same weight.
#' @param type_c Loss function for continuous outcomes. "square" (square error
#'   loss, the default) or "absolute" (absolute error loss).
#' @param type_b Loss function for binary outcomes. "deviance" (binomial
#'   deviance, the default) or "classification" (prediction error).
#' @param marginal For the mixed outcomes case, set to 0 (the default) to
#'   account for both the continuous loss and the binary loss, set to 1 to only
#'   account for the continuous loss, and set to 2 to only account for the
#'   binary loss.
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
#'   \item{PenScore}{Vector of penalized score functions at the
#'   estimated regression coefficients. If the algorithm converges, then these
#'   should be close to zero.}
#'   \item{FDR}{Estimated FDR for family=="Mixed", if requested.}
#'   \item{lambda.loss}{Cross validation loss (error) for the
#'   optimal tuning parameter(s) lambda, averaged across folds.}
#'   \item{LossMat}{Matrix of cross validation losses. Rows denote tuning
#'   parameter values, columns denote folds.}
#' @examples
#' \dontrun{
#' # Gaussian
#' N <- 100
#' m <- 10
#' p <- 50
#' y <- rnorm(N * m)
#' # If you want standardize = TRUE, you must provide an intercept.
#' X <- cbind(1, matrix(rnorm(N * m * (p - 1)), N * m, p - 1))
#' gr1 <- seq(0.001, 0.1, length.out = 100)
#' fit <- cv.pgee(X = X, y = y, N = N, m = m, grid1 = gr1, wctype = "CS",
#'             family = "Gaussian")
#'
#' # Binary
#' y <- sample(0:1, N*m, replace = TRUE)
#' fit <- cv.pgee(X = X, y = y, N = N, m = m, grid1 = gr1, wctype = "CS",
#'             family = "Binomial")
#'
#' # Bivariate mixed outcomes
#' # Generate some data
#' Bc <- c(2.0, 3.0, 1.5, 2.0, rep(0,times=p-4))
#' Bb <- c(0.7, -0.7, -0.4, rep(0,times=p-3))
#' dat <- gen_mixed_data(Bc, Bc, N, 0.5)
#' # We require two grids of tuning parameters
#' gr2 <- seq(0.0001, 0.01, length.out = 100)
#' # Estimate regression coefficients and false discovery rate
#' fit <- cv.pgee(X = dat$X, Z = dat$Z, yc = dat$yc, yb = dat$yb, N = N, m = 2,
#'                wctype = "CS", family = "Mixed", grid1 = gr1, grid2 = gr2,
#'                FDR = TRUE)
#' }
#' @export
cv.pgee <- function(N, m, X, Z = NULL, y = NULL, yc = NULL, yb = NULL, K = 5,
                    grid1, grid2 = NULL, wctype = "Ind", family = "Gaussian",
                    eps = 1e-6, maxiter = 1000, tol.coef = 1e-03,
                    tol.score = 1e-03, init = NULL, standardize = TRUE,
                    penalty = "SCAD", warm = TRUE, weights = rep(1, N),
                    type_c = "square",  type_b = "deviance", marginal = 0,
                    FDR = FALSE, fdr.corr = NULL, fdr.type = "all") {
  # Checks
  valid_wc_types <- c("Ind", "CS", "AR1")
  valid_family_types <- c("Gaussian", "Binomial", "Mixed")
  if (!(wctype %in% valid_wc_types))
    stop("Error in cv.pgee():
         Invalid Working Correlation Matrix type specified")
  if (!(family %in% valid_family_types))
    stop("Error in cv.pgee(): Invalid family specified")

  # y checks
  if (family == "Mixed") {
    if (any(is.null(yc), is.null(yb))) {
      stop("Error in cv.pgee(): For family==\"Mixed\", use arguments yc and yb
           as response vectors.")
    }
    if (!is.null(y)) {
      warning("For family==\"Mixed\", argument y is ignored.
              Arguments yc and yb are used as response vectors.")
    }
    y <- c(rbind(yc, yb))
    } else {
      if (is.null(y)) {
        stop("Error in cv.pgee(): Please provide response vector using argument y.")
      }
      if (any(!is.null(yc), !is.null(yb))) {
        warning("For family \"Gaussian\" or \"Binomial\", arguments yc and yb
              are ignored. Argument y is used as the response vector.")
      }
    }

  # Force intercept to be first column
  if (!check_intercept(X))
    stop("Error in cv.pgee(): In design matrix X, specify the column of ones
          for the intercept in the first column only.")
  if (!check_intercept(Z))
    stop("Error in cv.pgee(): In design matrix Z, specify the column of ones
          for the intercept in the first column only.")

  # For mixed, if Z isn't provided, use X
  if (is.null(Z) & family == "Mixed") Z <- X

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

  # CV
  if (family == "Gaussian" | family == "Binomial") {
    cv.pgee.nonmixed(X = X, y = y, N = N, m = m, K = K, wctype = wctype,
                     family = family, grid_lambda = grid1, eps = eps,
                     maxiter = maxiter, tol.coef = tol.coef,
                     tol.score = tol.score, init = init, penalty = penalty,
                     warm = warm, type_c = type_c, type_b = type_b,
                     weights = weights)
  } else if (family == "Mixed") {
    # Check both grids provied
    if(is.null(grid2)) {
      stop("Error in cv.pgee(): For family == \"Mixed\",
            both the grid1 and grid2 arguments must be provided.")
    }
    cv.pgee.mixed(X = X, Z = Z, yc = yc, yb = yb, N = N, K = K,
                  wctype = wctype, grid1 = grid1, grid2 = grid2, eps = eps,
                  maxiter = maxiter, tol.coef = tol.coef,
                  tol.score = tol.score, init = init, penalty = penalty,
                  warm = warm, type_c = type_c, type_b = type_b,
                  marginal = marginal, weights = weights, FDR = FDR,
                  fdr.corr = fdr.corr, fdr.type = fdr.type)
  } else stop("Error in cv.fit(): Argument family must be either
              \"Gaussian\", \"Binomial\", or \"Mixed\". ")
}
