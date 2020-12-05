\name{fastLm}
\alias{fastLm}
\alias{fastLm.default}
\alias{fastLm.formula}
\alias{fastLmPure}
\concept{regression}
\title{Bare-bones linear model fitting function}
\description{
  \code{fastLm} estimates the linear model using one of several methods
  implemented using the \code{Eigen} linear algebra library.
}
\usage{
fastLmPure(X, y, method = 0L)
fastLm(X, \dots)
\method{fastLm}{default}(X, y, method = 0L, \dots)
\method{fastLm}{formula}(formula, data = list(), method = 0L, \dots)
}
\arguments{

  \item{y}{the response vector}

  \item{X}{a model matrix}

  \item{formula}{an object of class \code{"\link{formula}"} (or one that
    can be coerced to that class): a symbolic description of the
    model to be fitted.  The details of model specification are given
    in the \sQuote{Details} section of the documentation for \code{\link{lm}}.}

  \item{data}{an optional data frame, list or environment (or object
    coercible by \code{\link{as.data.frame}} to a data frame) containing
    the variables in the model.  If not found in \code{data}, the
    variables are taken from \code{environment(formula)},
    typically the environment from which \code{lm} is called.}

  \item{method}{an integer scalar with value 0 for the column-pivoted QR
    decomposition, 1 for the unpivoted QR decomposition, 2 for the LLT
    Cholesky, 3 for the LDLT Cholesky, 4 for the Jacobi singular value
    decomposition (SVD) and 5 for a method based on the
    eigenvalue-eigenvector decomposition of
    \eqn{\mathbf{X}^\prime\mathbf{X}}{X'X}.  Default is zero.}

  \item{\dots}{not used}
}
\details{
  Linear models should be estimated using the \code{\link{lm}} function. In
  some cases, \code{\link{lm.fit}} may be appropriate.

  The \code{fastLmPure} function provides a reference use case of the \code{Eigen}
  C++ template library via the wrapper functions in the \pkg{RcppEigen} package.

  The \code{fastLm} function provides a more standard implementation of
  a linear model fit, offering both a default and a formula interface as
  well as \code{print}, \code{summary} and \code{predict} methods.

  Internally the \code{fastLm} function, by default, uses a QR
  decomposition with column pivots, which is a rank-revealing
  decomposition, so that it can handle rank-deficient cases
  effectively.  Other methods for determining least squares solutions
  are available according to the value of the \code{method} argument.

  An example of the type of situation requiring extra care in checking
  for rank deficiency is a two-way layout with missing cells (see the
  examples section).  These cases require a special pivoting scheme of
  \dQuote{pivot only on (apparent) rank deficiency} which is not part of
  conventional linear algebra software.
}
\value{
  \code{fastLmPure} returns a list with several components:
  \item{coefficients}{a vector of coefficients}
  \item{se}{a vector of the standard errors of the coefficient estimates}
  \item{rank}{a scalar denoting the computed rank of the model matrix}
  \item{df.residual}{a scalar denoting the degrees of freedom in the model}
  \item{residuals}{the vector of residuals}
  \item{s}{a numeric scalar - the root mean square for residuals}
  \item{fitted.values}{the vector of fitted value}

  \code{fastLm} returns a richer object which also includes the
     call argument similar to the \code{\link{lm}} or
  \code{\link[MASS]{rlm}} functions..
}
\seealso{\code{\link{lm}}, \code{\link{lm.fit}}}
\references{
  Douglas Bates and Dirk Eddelbuettel (2013). Fast and Elegant Numerical
  Linear Algebra Using the \pkg{RcppEigen} Package. \emph{Journal of
    Statistical Software}, \bold{52(5)}, 1-24.
  URL http://www.jstatsoft.org/v52/i05/.
}
\author{
  Eigen is described at \url{http://eigen.tuxfamily.org/index.php?title=Main_Page}.
  RcppEigen is written by Douglas Bates, Dirk Eddelbuettel and Romain Francois.
}
\examples{
  data(trees, package="datasets")
  mm <- cbind(1, log(trees$Girth))   # model matrix
  y  <- log(trees$Volume)            # response

  ## bare-bones direct interface
  flm <- fastLmPure(mm, y)
  print(flm)

  ## standard R interface for formula or data returning object of class fastLm
  flmmod <- fastLm( log(Volume) ~ log(Girth), data=trees)
  summary(flmmod)

  ## case where non-rank-revealing methods break down
  dd <- data.frame(f1 = gl(4, 6, labels = LETTERS[1:4]),
                   f2 = gl(3, 2, labels = letters[1:3]))[-(7:8), ]
  xtabs(~ f2 + f1, dd)     # one missing cell
  mm <- model.matrix(~ f1 * f2, dd)
  kappa(mm)                # large, indicating rank deficiency
  set.seed(1)
  dd$y <- mm \%*\% seq_len(ncol(mm)) + rnorm(nrow(mm), sd = 0.1)
  summary(lm(y ~ f1 * f2, dd))     # detects rank deficiency
  try(summary(fastLm(y ~ f1 * f2, dd))) # also detects rank deficiency
}
\keyword{regression}
