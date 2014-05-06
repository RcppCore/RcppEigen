\name{NEWS}
\title{News for Package 'RcppEigen}
\newcommand{\cpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\section{Changes in RcppEigen version 0.3.2.1.2 (2014-05-05)}{
  \itemize{
    \item Applied upstream patch to not trigger g++ UBSAN warning on
    the cpu id comparison, with thanks to Gael Guennebaud for the patch
  }
}
    
\section{Changes in RcppEigen version 0.3.2.1.1 (2014-03-06)}{
  \itemize{
    \item Better \code{ifdef} on one directory entry feature, with
    thanks to Brian Ripley.
  }
}

\section{Changes in RcppEigen version 0.3.2.1.0 (2014-03-03)}{
  \itemize{
    \item Updated to version 3.2.0 of Eigen
  }
}

\section{Changes in RcppEigen version 0.3.2.0.3 (2014-03-01)}{
  \itemize{
    \item Updated and extended \code{RcppEigen.package.skeleton()} to
    use several examples via \cpkg{Rcpp} Attributes; also removed the
    deprecated \code{namespace} argument 
    \item Updated skeleton package example for \cpkg{Rcpp} 0.11.0 or
    later by removing needed for linking with user library
    \item Updated files \code{DESCRIPTION}, \code{NAMESPACE},
    \code{src/Makevars} and \code{src/Makevars.win} similarly
  }
}

\section{Changes in RcppEigen version 0.3.2.0.2 (2014-01-26)}{
  \itemize{
    \item Converted three unused unit test files to \cpkg{RUnit} and
    removed \code{Suggests:} of \cpkg{testthat}
    \item Add declaration to import a symbol from \cpkg{Rcpp} to
    \code{NAMESPACE} to ensure proper instantiation with the upcoming
    \cpkg{Rcpp} version
    \item Retire \code{SHLIB.maker} function
  }
}

\section{Changes in RcppEigen version 0.3.2.0.1 (2013-12-18)}{
  \itemize{
    \item New maintainer -- with a big thanks to Doug for all his work
    \item Applied two small patches to deal with non-g++ compilrs
    \item Clarifications concerning license and authorship of
    Eigen (as opposed to RcppEigen) code added to \code{DESCRIPTION} at
    the request of CRAN 
  }
}

\section{Changes in RcppEigen version 0.3.2.0 (2013-11-13)}{
  \itemize{
    \item Update to version 3.2.0 of Eigen
  }
}

\section{Changes in RcppEigen version 0.3.1.2.3 (2013-10-25)}{
  \itemize{
    \item Fix to RcppEigenCholmod.h to incorporate changes in the
    cholmod_factor struct. These changes are necessary if code
    compiled against RcppEigen that uses CHOLMOD factors is to be run
    with versions of the Matrix package >= 1.1-0
  }
}

\section{Changes in RcppEigen version 0.3.1.2 (2012-11-29)}{
  \itemize{
    \item Upgraded to Eigen 3.1.2
    \item Fixes to RcppEigenWrap.h and adjustment of tests accordingly.
    The changes allow RowMajor matrices to be wrapped (thanks to Gael
    Guennebaud) but cannot handle RowVector types. There will need to be
    more template metaprogramming done to redirect the case of RowVector,
    which cannot be changed to a ColMajor form.
    \item Because of changes in R, -DNDEBUG is automatic. One must
    override it with -UNDEBUG in the local ~/.R/Makevars to activate the
    debugging code. 
    \item New (unexported) functions CxxFlags() and RcppEigenCxxFlags()
    for use in Makefiles
    \item Fixes related to Rcpp 0.10.*
  }
}

\section{Changes in RcppEigen version 0.3.1 (2012-08-07)}{
  \itemize{
    \item Upgraded to Eigen 3.1.0
    \item Removed the "unsupported" Eigen module AutoDiff which defined a
    macro "sign" that conflicted with a function in the R API (which
    really should be visible as "Rf_sign", not sure why it shows up as
    "sign" and don't have time to investigate)
    \item Commented out several tests involving complex vectors and
    matrices.  Again there are compilation problems related to
    conflicting definitions in the std:: namespace and the R API and
    Eigen, which I don't have time to investigate.
  }
}

\section{Changes in RcppEigen version 0.2.0 (2012-03-12)}{
  \itemize{
    \item Upgraded the version of Eigen to 3.1.0-alpha2, in which the sparse
    matrix modules are now in the "supported" tree.
    \item Added several "unsupported" Eigen modules including
    \itemize{
      \item AutoDiff (a small automatic differentiation package adapted to
      vectors and matrices)
      \item IterativeSolvers (iterative linear and nonlinear solver algorithms)
      \item KroneckerProduct (as the name implies)
      \item MatrixFunctions (matrix cos, exp, log, sin, sinh, etc.)
      \item NonlinearOptimization (based on minpack but uses reverse
      communication - yay!)
      \item NumericalDiff (numerical differentiation of vector-valued or
      matrix-valued functions)
      \item Polynomials (polynomial representation and solution using a QR
      algorithm)
      \item Skyline (sparse skyline matrices useful in finite-element codes)
      \item SparseExtra (dynamic sparse matrices, now deprecated, and Matrix
      Market I/O functions)
      \item Splines (multidimensional spline representations and spline
      interpolation)
    }
    \item At present all these modules, including the MatrixFunctions
    module, are included with RcppEigen.h but that may change if too
    many people get unexpected results from A.exp()
    \item The ability to wrap RowMajor sparse matrices and to use
    as<Eigen::ArrayXd> etc.
    \item Migrated some tests to the testthat package.  Currently there is
    some difficulty with combining testthat, inline and R CMD check.
  }
}
 
  
