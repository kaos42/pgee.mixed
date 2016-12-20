# include <RcppArmadillo.h>

// v = (v1,...,vN) is the stacked NM length vector of variance fns, without phi
// Note that C++ matrix indexing starts from 0
// Note that phi is not required, as it gets cancelled - THIS IS NOT TRUE unless corr is estimated as in Liang and Zeger.
// Since we use biserial estimation for correlation, we must include phi
// This function is only to be used for responses of the same type (not mixed). Phi can then just be divided at the end.
// [[Rcpp::export]]
arma::vec CppScore(arma::mat X, arma::vec v, arma::mat Rhat, arma::vec y,
        arma::vec mu, int N, float phi, arma::vec w) {
    int m = Rhat.n_rows;
    int p = X.n_cols;
    arma::vec score = arma::zeros<arma::vec>(p);
    for(int i=1; i<N+1; i++) {
        arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
        arma::vec y_i = y.subvec((i-1)*m, (i*m)-1);
        arma::vec mu_i = mu.subvec((i-1)*m, (i*m)-1);
        arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);

        score += w(i-1) * X_i.t() * diagmat(sqrt(v_i)) * solve( Rhat,
                diagmat(sqrt(1/v_i)) * (y_i-mu_i) );
    }
    score /= arma::accu(w);
    score /= phi;
    return(score);
}

// v = (v1,...,vN) is the stacked NM length vector of variance fns, without phi
// Note that C++ matrix indexing starts from 0
// This function is used only for responses of same type (not mixed). Thus, phi can just be divided at the end.
// [[Rcpp::export]]
arma::mat CppHess(arma::mat X, arma::vec v, arma::mat Rhat, float phi, int N, arma::vec w) {
    int m = Rhat.n_rows;
    int p = X.n_cols;
    arma::mat H = arma::zeros<arma::mat>(p,p);
    for(int i=1; i<N+1; i++) {
        arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
        arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);

        H += w(i-1) * X_i.t() * diagmat(sqrt(v_i)) * solve( Rhat,
                diagmat(sqrt(v_i)) * X_i );
    }
    H /= arma::accu(w);
    H /= phi;
    return(H);
}

// Meat of sandwich estimator
// [[Rcpp::export]]
arma::mat CppM(arma::vec y, arma::mat X, arma::vec mu, arma::vec v,
      arma::mat Rhat, float phi, int N, arma::vec w) {
    int m = Rhat.n_rows;
    int p = X.n_cols;
    arma::mat M = arma::zeros<arma::mat>(p,p);
    for(int i=1; i<N+1; i++) {
        arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
        arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);
        arma::vec y_i = y.subvec((i-1)*m, (i*m)-1);
        arma::vec mu_i = mu.subvec((i-1)*m, (i*m)-1);
        arma::vec e_i = diagmat(sqrt(1/v_i)) * (y_i-mu_i);
        arma::mat L_i = X_i.t() * diagmat(sqrt(v_i)) * solve(Rhat, e_i);
        M += w(i-1) * L_i * L_i.t();
    }
    M /= (arma::accu(w)*arma::accu(w));
    return(M);
}

// The updation term in the NewtonRapson Beta update
// [[Rcpp::export]]
arma::vec CppIncBeta(Rcpp::NumericVector beta,
        Rcpp::NumericVector score,
        Rcpp::NumericMatrix Hess,
        Rcpp::NumericMatrix Emat, int n) {

    arma::vec b(beta.begin(), beta.size(), false);
    arma::vec s(score.begin(), score.size(), false);
    arma::mat H(Hess.begin(), Hess.nrow(), Hess.ncol(), false);
    arma::mat E(Emat.begin(), Emat.nrow(), Emat.ncol(), false);

    //arma::vec d = s - (n * (E * b));
    arma::vec d = s - ((E * b));
    int p = H.n_cols;
    //arma::mat A = H + (n * E) + 1e-6*arma::eye<arma::mat>(p,p);
    arma::mat A = H + (E) + 1e-6*arma::eye<arma::mat>(p,p);
    arma::vec inc = arma::solve(A, d);
    return(inc);
}

// [[Rcpp::export]]
float CppAlphaCS(arma::vec pres, arma::vec w, int m, int p, float phi) {

    // Since all clusters have equal size:
    int N = int(pres.size() / m);
    float alpha_hat = 0.0;
    // float denom = 0.5*m*(m-1)*N - p;
    float denom = 0.5*m*(m-1)*arma::accu(w);
    arma::vec pres_i = arma::zeros<arma::vec>(m);
    for(int i=0; i<N; ++i) {
        pres_i = pres.subvec(i*m, ((i+1)*m)-1);
        for(int j=0; j<m; ++j) {
            for(int k=j+1; k<m; ++k) {
                alpha_hat += w(i) * (pres_i(j) * pres_i(k));
            }
        }
    }
    return(alpha_hat / (denom * phi));
}

// [[Rcpp::export]]
float CppAlphaAR1(arma::vec pres, arma::vec w, int m, int p, float phi) {
    // Since all clusters have equal size:
    int N = int(pres.size() / m);
    float alpha_hat = 0.0;
    // float denom = N*(m-1) - p;
    float denom = arma::accu(w)*(m-1);
    arma::vec pres_i = arma::zeros<arma::vec>(m);
    for(int i=0; i<N; ++i) {
        pres_i = pres.subvec(i*m, ((i+1)*m)-1);
        for(int j=0; j<m-1; ++j) {
            alpha_hat += w(i) * pres_i(j) * pres_i(j+1);
        }
    }
    return(alpha_hat / (denom * phi));
}

// The penalized score to check solution
// [[Rcpp::export]]
arma::vec PenScore(Rcpp::NumericVector beta,
        Rcpp::NumericVector score,
        Rcpp::NumericMatrix Emat, int n) {

    arma::vec b(beta.begin(), beta.size(), false);
    arma::vec s(score.begin(), score.size(), false);
    arma::mat E(Emat.begin(), Emat.nrow(), Emat.ncol(), false);

    //arma::vec d = s - (n * (E * b));
    arma::vec d = s - ((E * b));
    return(d);
}
// check if two matrices are identical
// [[Rcpp::export]]
bool samemats(arma::mat X, arma::mat Y) {
    // Get dimensions
    int px = X.n_cols;
    int py = Y.n_cols;
    int nx = X.n_rows;
    int ny = Y.n_rows;
    if(px != py)
    {
      return(false);
    }
    if(nx != ny)
    {
      return(false);
    }
    // check each element
    for(int i=0; i<nx; ++i)
    {
      for(int j=0; j<px; ++j)
      {
        if(X(i,j) != Y(i,j))
        {
          return(false);
        }
      }
    }
    return(true);
}

// Function to compute W matrix required for FDR procedure
// X is in transformed form for mixed
// [[Rcpp::export]]
arma::mat CppW(arma::mat X, arma::vec v, arma::mat Rhat, int N) {
    int m = Rhat.n_rows;
    arma::mat I(m,m,arma::fill::eye);
    // for no correlation, W = X
    if(samemats(Rhat,I)) {
      return(X);
    }
    int p = X.n_cols;
    arma::mat W = arma::mat(N*m, p);
    for(int i=1; i<N+1; i++) {
        arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
        arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);
        W( arma::span((i-1)*m, (i*m)-1), arma::span::all ) =
                diagmat(sqrt(1/v_i)) * solve( Rhat, diagmat(sqrt(v_i))*X_i );
    }
    return(W);
}

// Hessian for bivariate mixed responses, with dispersion, with weights
// [[Rcpp::export]]
arma::mat CppHess2(arma::mat X, arma::vec v, float alpha, float phi, int N, arma::vec w) {
    int m = 2;
    int p = X.n_cols;
    // precompute the inverse of the correlation matrix
    arma::mat Rinv = arma::mat(2,2);
    Rinv(0,0) = 1.0;
    Rinv(0,1) = -alpha;
    Rinv(1,0) = -alpha;
    Rinv(1,1) = 1.0;
    Rinv = Rinv/(1-alpha*alpha);
    // P is the inverse of the square root of the diagonal dispersion matrix
    arma::mat P = arma::mat(2,2);
    P(0,0) = 1/sqrt(phi); // phi is sigma^2, the continuous dispersion
    P(0,1) = 0.0;
    P(1,0) = 0.0;
    P(1,1) = 1.0; // (inverse sqrt of) binary dispersion is unity
    arma::mat H = arma::zeros<arma::mat>(p,p);
    for(int i=1; i<N+1; i++) {
        arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
        arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);

        H += w(i-1) * X_i.t() * diagmat(sqrt(v_i)) * P * Rinv * P * diagmat(sqrt(v_i)) * X_i;
    }
    H /= arma::accu(w);
    return(H);
}

// Score for bivariate mixed responses, with dispersion, with weights
// [[Rcpp::export]]
arma::vec CppScore2(arma::mat X, arma::vec v, float alpha, arma::vec y,
                   arma::vec mu, int N, float phi, arma::vec w) {
  int m = 2;
  int p = X.n_cols;
  // precompute the inverse of the correlation matrix
  arma::mat Rinv = arma::mat(2,2);
  Rinv(0,0) = 1.0;
  Rinv(0,1) = -alpha;
  Rinv(1,0) = -alpha;
  Rinv(1,1) = 1.0;
  Rinv = Rinv/(1-alpha*alpha);
  // P is the inverse of the square root of the diagonal dispersion matrix
  arma::mat P = arma::mat(2,2);
  P(0,0) = 1/sqrt(phi); // phi is sigma^2, the continuous dispersion
  P(0,1) = 0.0;
  P(1,0) = 0.0;
  P(1,1) = 1.0; // (inverse sqrt of) binary dispersion is unity
  arma::vec score = arma::zeros<arma::vec>(p);
  for(int i=1; i<N+1; i++) {
    arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
    arma::vec y_i = y.subvec((i-1)*m, (i*m)-1);
    arma::vec mu_i = mu.subvec((i-1)*m, (i*m)-1);
    arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);

    score += w(i-1) * X_i.t() * diagmat(sqrt(v_i)) * P * Rinv * P * diagmat(sqrt(1/v_i)) * (y_i-mu_i);
  }
  score /= arma::accu(w);
  return(score);
}

// Meat of sandwich estimator, with weights
// [[Rcpp::export]]
arma::mat CppM2(arma::vec y, arma::mat X, arma::vec mu, arma::vec v,
               float alpha, float phi, int N, arma::vec w) {
  int m = 2;
  int p = X.n_cols;
  arma::mat M = arma::zeros<arma::mat>(p,p);
  // precompute the inverse of the correlation matrix
  arma::mat Rinv = arma::mat(2,2);
  Rinv(0,0) = 1.0;
  Rinv(0,1) = -alpha;
  Rinv(1,0) = -alpha;
  Rinv(1,1) = 1.0;
  Rinv = Rinv/(1-alpha*alpha);
  // P is the inverse of the square root of the diagonal dispersion matrix
  arma::mat P = arma::mat(2,2);
  P(0,0) = 1/sqrt(phi); // phi is sigma^2, the continuous dispersion
  P(0,1) = 0.0;
  P(1,0) = 0.0;
  P(1,1) = 1.0; // (inverse sqrt of) binary dispersion is unity
  for(int i=1; i<N+1; i++) {
    arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
    arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);
    arma::vec y_i = y.subvec((i-1)*m, (i*m)-1);
    arma::vec mu_i = mu.subvec((i-1)*m, (i*m)-1);
    arma::mat L_i = X_i.t() * diagmat(sqrt(v_i)) * P * Rinv * P * diagmat(sqrt(1/v_i)) * (y_i-mu_i);
    M += w(i-1) * L_i * L_i.t();
  }
  M /= (arma::accu(w)*arma::accu(w));
  return(M);
}

// Function to compute W matrix required for FDR procedure
// X is in transformed form for mixed
// allows weights
// [[Rcpp::export]]
arma::mat CppW2(arma::mat X, arma::vec v, arma::mat Rhat, int N, arma::vec w) {
  int m = Rhat.n_rows;
  arma::mat I(m,m,arma::fill::eye);
  // for no correlation, W = X
  if(samemats(Rhat,I)) {
    return(X);
  }
  int p = X.n_cols;
  arma::mat W = arma::mat(N*m, p);
  for(int i=1; i<N+1; i++) {
    arma::mat X_i = X( arma::span((i-1)*m, (i*m)-1), arma::span::all );
    arma::vec v_i = v.subvec((i-1)*m, (i*m)-1);
    W( arma::span((i-1)*m, (i*m)-1), arma::span::all ) = w(i-1) *
      diagmat(sqrt(1/v_i)) * solve( Rhat, diagmat(sqrt(v_i))*X_i );
  }
  return(W);
}
