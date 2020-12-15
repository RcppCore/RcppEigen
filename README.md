## RcppEigen: Rcpp Integration for the Eigen Templated Linear Algebra Library

[![Build Status](https://travis-ci.org/RcppCore/RcppEigen.svg)](https://travis-ci.org/RcppCore/RcppEigen)
[![CI](https://github.com/RcppCore/RcppEigen/workflows/ci/badge.svg)](https://github.com/RcppCore/RcppEigen/actions?query=workflow%3Aci)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-2.0.html)
[![License](https://img.shields.io/badge/license-MPL2-brightgreen.svg?style=flat)](https://www.mozilla.org/MPL/2.0/)
[![CRAN](https://www.r-pkg.org/badges/version/RcppEigen)](https://cran.r-project.org/package=RcppEigen)
[![Dependencies](https://tinyverse.netlify.com/badge/RcppEigen)](https://cran.r-project.org/package=RcppEigen)
[![Debian package](https://img.shields.io/debian/v/r-cran-rcppeigen/sid?color=brightgreen)](https://packages.debian.org/sid/r-cran-rcppeigen)
[![Last Commit](https://img.shields.io/github/last-commit/RcppCore/RcppEigen)](https://github.com/RcppCore/RcppEigen)  
[![Downloads](https://cranlogs.r-pkg.org/badges/RcppEigen?color=brightgreen)](https://www.r-pkg.org:443/pkg/RcppEigen)
[![CRAN use](https://jangorecki.gitlab.io/rdeps/RcppEigen/CRAN_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppEigen)
[![BioConductor use](https://jangorecki.gitlab.io/rdeps/RcppEigen/BioC_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppEigen)
[![StackOverflow](https://img.shields.io/badge/stackoverflow-rcpp-orange.svg)](https://stackoverflow.com/questions/tagged/rcpp)
[![JSS](https://img.shields.io/badge/JSS-10.18637%2Fjss.v052.i05-brightgreen)](https://dx.doi.org/10.18637/jss.v052.i05)


### Synopsis

[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a C++ template library for linear algebra:
matrices, vectors, numerical solvers and related algorithms.  It supports dense and sparse
matrices on integer, floating point and complex numbers, decompositions of such matrices,
and solutions of linear systems. Its performance on many algorithms is comparable with
some of the best implementations based on `Lapack` and level-3 `BLAS`.

RcppEigen provides an interface from R to and from [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) by
using the facilities offered by the [Rcpp](http://dirk.eddelbuettel.com/code/rcpp.html)
package for seamless R and C++ integration.

### Examples

A few examples are over at the [Rcpp Gallery](https://gallery.rcpp.org/tags/eigen/). A simple one is

```c++
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers

// [[Rcpp::export]]
VectorXd getEigenValues(Map<MatrixXd> M) {
    SelfAdjointEigenSolver<MatrixXd> es(M);
    return es.eigenvalues();
}
```

which can be turned into a function callable from R via a simple

```
sourceCpp("eigenExample.cpp")
```

due to the two Rcpp directives to use headers from the RcppEigen package, and to export
the `getEigenValues()` function -- but read [the full
post](https://gallery.rcpp.org/articles/eigen-eigenvalues/) for details.


### Status

The package is mature and under active development, following the
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) release cycle.

### Documentation

The package contains a pdf vignette which is a pre-print of the [paper by
Bates and Eddelbuettel](https://www.jstatsoft.org/article/view/v052i05)
in JSS (2013, v52i05).

### Authors

Douglas Bates, Dirk Eddelbuettel, Romain Francois, and Yixuan Qiu

### License

GPL (>= 2)
