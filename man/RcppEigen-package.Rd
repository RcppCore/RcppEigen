\name{RcppEigen-package}
\alias{RcppEigen-package}
\alias{RcppEigen}
\docType{package}
\title{
Rcpp/Eigen bridge
}
\description{
  The package eases the use of the Eigen C++ template library for linear
  algebra with Rcpp
}
\details{
  This package contains the header files for the Eigen C++ template
  library.  The typical usage is to install this package and list it in
  the \env{LinkingTo: } line in the \file{DESCRIPTION} file of
  other packages.  The C++ source code and the R source code in this
  package are for illustration only.

  As described at the \href{https://libeigen.gitlab.io/}{Eigen project home
  page} , Eigen is a a C++ template library for linear algebra: matrices,
  vectors, numerical solvers, and related algorithms.
}
\section{Threading}{
  The Eigen library can take advantage of OpenMP to execute computations in
  parallel via multi-threaded code. The number of cores uses can be set (or
  retrieved) explicitly via helper functions \code{EigenSetNbThreads()} and
  \code{EigenNbThreads()}. A default value is stored at package startup; it
  recognises R option value \code{Ncpus} and environment variable
  \code{OMP_THREAD_LIMITS}.  Additional helper functions
  \code{RcppEigen_throttle_cores()} and \code{RcppEigen_reset_cores()} are
  available to (temporarily) lower the number of cores uses and to reset to
  the package default value set at startup.
}
\seealso{
  \code{\link{RcppEigen_throttle_cores}}
}
\references{
  Douglas Bates and Dirk Eddelbuettel (2013). Fast and Elegant Numerical
  Linear Algebra Using the \pkg{RcppEigen} Package. \emph{Journal of
    Statistical Software}, \bold{52(5)}, 1-24.
  URL http://www.jstatsoft.org/v52/i05/.
}
\keyword{ package }
