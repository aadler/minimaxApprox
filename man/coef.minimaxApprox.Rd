% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{coef.minimaxApprox}
\alias{coef.minimaxApprox}
\encoding{UTF-8}

\title{Extract coefficients from a \code{"minimaxApprox"} object}

\description{
Extracts the numerator and denominator vectors from a \code{"minimaxApprox"}
object. For objects with both Chebyshev and monomial coefficients, it will
extract both.
}

\usage{
\method{coef}{minimaxApprox}(object, \dots)
}

\arguments{
    \item{object}{An object inheriting from \link{class} \code{"minimaxApprox"}.
}
    \item{\dots}{Other arguments.}
}

\value{Coefficients extracted from the \code{"minimaxApprox"} object. A
\link{list} containing:

    \item{a}{The polynomial coefficients or the rational numerator
    coefficients.}
    \item{b}{The rational denominator coefficients. Missing for polynomial
    approximation.}
    \item{aMono}{The polynomial coefficients or the rational numerator
    coefficients for the monomial basis when the approximation was done using
    Chebyshev polynomials. Missing if only the monomial basis was used.}
    \item{bMono}{The rational denominator coefficients for the monomial basis
    when the approximation was done using Chebyshev polynomials. Missing if
    either only the monomial basis was used or for polynomial approximation.}
}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\seealso{
\code{\link{minimaxApprox}}
}

\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
coef(PP)
identical(unlist(coef(PP), use.names = FALSE), c(PP$a, PP$aMono))

RR <- minimaxApprox(exp, 0, 1, c(2, 3), basis = "m")
coef(RR)
identical(coef(RR), list(a = RR$a, b = RR$b))
}

\keyword{methods}
\keyword{NumericalMathematics}
