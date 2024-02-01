% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{print.minimaxApprox}
\alias{print.minimaxApprox}
\encoding{UTF-8}

\title{Print method for a \code{"minimaxApprox object"}}

\description{
Provides a more human-readable output of a \code{"minimaxApprox"} object.
}

\usage{
\method{print}{minimaxApprox}(x, digits = 6L, \dots)
}

\arguments{
    \item{x}{An object inheriting from \link{class} \code{"minimaxApprox"}.}
    \item{digits}{integer; Number of digits to which to round the ratio.}
    \item{\dots}{Further arguments to \code{print}.}
}

\details{
To print the raw \code{"minimaxApprox"} object use \code{print.default}.
}

\value{No return value; called for side effects.}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\seealso{
\code{\link{minimaxApprox}}
}

\examples{
PP <- minimaxApprox(sin, 0, 1, 8)
PP
print(PP, digits = 2L)
print.default(PP)
}

\keyword{print}
\keyword{methods}
\keyword{NumericalMathematics}
