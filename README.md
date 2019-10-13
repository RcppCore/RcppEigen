## RcppEigen

[![Build Status](https://travis-ci.org/RcppCore/RcppEigen.svg)](https://travis-ci.org/RcppCore/RcppEigen) 
[![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) 
[![License](http://img.shields.io/badge/license-MPL2-brightgreen.svg?style=flat)](http://www.mozilla.org/MPL/2.0/) 
[![CRAN](http://www.r-pkg.org/badges/version/RcppEigen)](https://cran.r-project.org/package=RcppEigen) 
[![Dependencies](https://tinyverse.netlify.com/badge/RcppEigen)](https://cran.r-project.org/package=RcppEigen)  
[![Downloads](http://cranlogs.r-pkg.org/badges/RcppEigen?color=brightgreen)](http://www.r-pkg.org/pkg/RcppEigen) 
[![CRAN use](https://jangorecki.gitlab.io/rdeps/RcppEigen/CRAN_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppEigen)
[![BioConductor use](https://jangorecki.gitlab.io/rdeps/RcppEigen/BioC_usage.svg?sanitize=true)](https://cran.r-project.org/package=RcppEigen)
[![StackOverflow](https://img.shields.io/badge/stackoverflow-rcpp-orange.svg)](https://stackoverflow.com/questions/tagged/rcpp)

### Overview

[Eigen](http://eigen.tuxfamily.org) is a C++ template library for linear
algebra: matrices, vectors, numerical solvers and related algorithms.  It
supports dense and sparse matrices on integer, floating point and complex
numbers, decompositions of such matrices, and solutions of linear
systems. Its performance on many algorithms is comparable with some of the
best implementations based on `Lapack` and level-3 `BLAS`.

The RcppEigen package includes the header files from the Eigen C++
template library (currently version 3.3.5). Thus users do not need to
install Eigen itself in order to use RcppEigen.

### Status

The package is mature and under active development, following the
[Eigen](http://eigen.tuxfamily.org) release cycle.

### Documentation

The package contains a pdf vignette which is a pre-print of the [paper by
Bates and Eddelbuettel](https://www.jstatsoft.org/article/view/v052i05)
in JSS (2013, v52i05).

### Authors

Douglas Bates, Dirk Eddelbuettel, Romain Francois, and Yixuan Qiu

### License

GPL (>= 2)
