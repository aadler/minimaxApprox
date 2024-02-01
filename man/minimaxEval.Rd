% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{minimaxEval}
\alias{minimaxEval}
\encoding{UTF-8}

\title{Evaluate Minimax Approximation}
\description{
Evaluates the rational or polynomial approximation stored in \code{mmA} at
\code{x}.
}
\usage{
minimaxEval(x, mmA, basis = "Chebyshev")
}
\arguments{
  \item{x}{a numeric vector}
  \item{mmA}{a \code{"minimaxApprox"} return object}
  \item{basis}{character; Which polynomial basis to use in to evaluate the
  function; see \code{\link{minimaxApprox}} for more details. If Chebyshev is
  requested but the analysis used only monomials, the calculation will proceed
  using the monomials with a message. The default is "Chebyshev", and the
  parameter is case-insensitive and may be abbreviated.}
}
\details{
This is a convenience function to evaluate the approximation at \code{x}.
}
\value{
A vector of the same length as \code{x} containing the approximated values.
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}
\seealso{\code{\link{minimaxApprox}}, \code{\link{minimaxErr}}}
\examples{
# Show results
x <- seq(0, 0.5, length.out = 11L)
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
apErr <- abs(exp(x) - minimaxEval(x, mmA))
all.equal(max(apErr), mmA$ExpErr)

# Plot results
curve(exp, 0.0, 0.5, lwd = 2)
curve(minimaxEval(x, mmA), 0.0, 0.5, add = TRUE, col = "red", lty = 2L, lwd = 2)
}

\keyword{NumericalMathematics}
