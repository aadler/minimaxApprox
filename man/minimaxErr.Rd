\name{minimaxErr}
\alias{minimaxErr}
\title{
Evaluate the Minimax Approximation Error
}
\description{
Evaluates the difference between the function and the minimax approximation at
\code{x}.
}
\usage{
minimaxErr(x, mmA)
}
\arguments{
  \item{x}{a numeric vector}
  \item{mmA}{a \code{"minimaxApprox"} return object}
}
\details{
This is a convenience function to evaluate the approximation error at \code{x}.
It will use the same polynomial basis as was used in the approximation; see
\code{\link{minimaxApprox}} for more details.
}
\value{
A vector of the same length as \code{x} containing the approximation error
values.
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}
\seealso{\code{\link{minimaxApprox}}, \code{\link{minimaxEval}}}
\examples{
# Show results
x <- seq(0, 0.5, length.out = 11L)
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
err <- minimaxEval(x, mmA) - exp(x)
all.equal(err,  minimaxErr(x, mmA))

# Plot results
x <- seq(0, 0.5, length.out = 1001L)
plot(x, minimaxErr(x, mmA), type = "l")
}

\keyword{NumericalMathematics}
