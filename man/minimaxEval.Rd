\name{minimaxEval}
\alias{minimaxEval}
\encoding{UTF-8}

\title{Evaluate Minimax Approximation}
\description{
The function evaluates the rational or polynomial approximation at x.
}
\usage{
minimaxEval(x, mmA)
}
\arguments{
  \item{x}{a numeric vector}
  \item{mmA}{a \code{"minimaxApprox"} return object}
}
\details{
This is a convenience function to evalute the approximation at supplied \code{x}.
}
\value{
A vector of the same length as \code{x} containing the approximation values.
}
\author{Avraham Adler \email{Avraham.Adler@gmail.com}}
\seealso{\code{\link{minimaxApprox}}}
\examples{
x <- seq(0.1, 0.4, 0.025)
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
exp(x)
minimaxEval(x, mmA)
exp(x) - minimaxEval(x, mmA)
all(exp(x) - minimaxEval(x, mmA) <= mmA$EE)
}

\keyword{NumericalMathematics}
