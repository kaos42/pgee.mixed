# pgee.mixed
Penalized Generalized Estimating Equations for Bivariate Mixed Outcomes

Perform simultaneous estimation and variable selection for correlated 
bivariate mixed outcomes (one continuous outcome and one binary outcome 
per cluster) using penalized generalized estimating equations. In 
addition, clustered Gaussian and binary outcomes can also be modeled. 
The SCAD, MCP, and LASSO penalties are supported. Cross-validation can 
be performed to find the optimal regularization parameter(s).

## Installation from GitHub:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("kaos42/pgee.mixed")
```

### Installation note for OS X users:

This package uses Rcpp and RcppArmadillo. Mac OS X users may see the following
error in their console while trying to install the package:

```
ld: warning: directory not found for option '-L/usr/local/lib/gcc/x86_64-apple-darwin13.0.0/4.8.2'  
ld: library not found for -lgfortran
```

The solution is documented [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/). In a nutshell, type the following in a terminal:

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2  
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

## To do list

If there is sufficient interest in this package, the following features could
be added:

* Families other than Gaussian and binomial.  
* Working correlation structures other than independence, compound symmetry, and AR(1).  
* Specify a vector of cluster ids rather than force the equal cluster size structure.  
* Users can provide fixed working correlation and dispersion parameters.  
* Users can provide an index vector specifying which parameters are not to be penalized.  
* Add a ridge component to existing penalties.  
