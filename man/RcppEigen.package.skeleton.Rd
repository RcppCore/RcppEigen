\name{RcppEigen.package.skeleton}
\alias{RcppEigen.package.skeleton}
\title{
Create a skeleton for a new package that intends to use RcppEigen
}
\description{
  \code{RcppEigen.package.skeleton} automates the creation of 
  a new source package that intends to use features of RcppEigen.
  
  It is based on the \link[utils]{package.skeleton} function
  which it executes first.
}
\usage{
RcppEigen.package.skeleton(name = "anRpackage", list = character(), 
	environment = .GlobalEnv, path = ".", force = FALSE, 
	code_files = character(), example_code = TRUE)
}
\arguments{
  \item{name}{See \link[utils]{package.skeleton}}
  \item{list}{See \link[utils]{package.skeleton}}
  \item{environment}{See \link[utils]{package.skeleton}}
  \item{path}{See \link[utils]{package.skeleton}}
  \item{force}{See \link[utils]{package.skeleton}}
  \item{code_files}{See \link[utils]{package.skeleton}}
  \item{example_code}{If TRUE, example C++ code using RcppEigen is added to the package}
}
\details{
  In addition to \link[utils]{package.skeleton} : 
  
  The \samp{DESCRIPTION} file gains a Depends line requesting that 
  the package depends on Rcpp and RcppEigen and 
  a LinkingTo line so that the package finds Rcpp and RcppEigen header files.
	
  The \samp{NAMESPACE} gains a \code{useDynLib} directive.
	
  The \samp{src} directory is created if it does not exists and 
  a \samp{Makevars} file is added setting the environment variable
  \samp{PKG_LIBS} to accomodate the necessary flags
  to link with the Rcpp library. 
	
  If the \code{example_code} argument is set to \code{TRUE}, 
  example files \samp{rcppeigen_hello_world.h} and \samp{rcppeigen_hello_world.cpp}
  are also created in the \samp{src}. An R file \samp{rcppeigen_hello_world.R} is 
  expanded in the \samp{R} directory, the \code{rcppeigen_hello_world} function
  defined in this files makes use of the C++ function \samp{rcppeigen_hello_world}
  defined in the C++ file. These files are given as an example and should 
  eventually by removed from the generated package.
}
\value{
  Nothing, used for its side effects
}
\seealso{
  \link[utils]{package.skeleton}
}
\references{
  Read the \emph{Writing R Extensions} manual for more details.

  Once you have created a \emph{source} package you need to install it:
  see the \emph{R Installation and Administration} manual,
  \code{\link{INSTALL}} and \code{\link{install.packages}}.
}
\examples{
\dontrun{
  RcppEigen.package.skeleton("foobar")
}
}
\keyword{ programming }

