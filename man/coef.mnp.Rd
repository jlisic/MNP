\name{coef.mnp}

\alias{coef.mnp}
\alias{coefficients.mnp}

\title{Extract Multinomial Probit Model Coefficients}

\description{
  \code{coef.mnp} is a function which extracts multinomial probit model
  coefficients from ojbects returned by
  \code{mnp}. \code{coefficients.mnp} is an alias for it. \code{coef}
  method for class \code{mnp}. 
}

\usage{
  \method{coef}{mnp}(object, subset = NULL, ...)
}

\arguments{
  \item{object}{An output object from \code{mnp}.}
  \item{subset}{A scalar or a numerical vector specifying the row
    number(s) of \code{param} in the output object from \code{mnp}. If
    specified, the posterior draws of coefficients for those rows are
    extracted. The default is \code{NULL} where all the posterior draws
    are extracted. 
  }
  \item{...}{further arguments passed to or from other methods.}
}

\value{
  \code{coef.mnp} returns a matrix (when a numerical vector or
  \code{NULL} is specified for \code{subset} argument) or a vector (when
  a scalar is specified for \code{subset} arugment) of multinomila
  probit model coefficients.
}

\seealso{\code{mnp}, \code{cov.mnp}; MNP home page at
  \url{http://imai.princeton.edu/research/MNP.html}}

\author{
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@Princeton.Edu}
}

\keyword{methods}
