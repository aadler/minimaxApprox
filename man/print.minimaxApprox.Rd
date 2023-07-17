\name{print.minimaxApprox}
\alias{print.minimaxApprox}

\title{Print method for a \dQuote{minimaxApprox} object}

\description{
Provides a more human-readable output of a \dQuote{minimaxApprox} object.
}

\usage{
\method{print}{minimaxApprox}(x, \dots)
}

\arguments{
  \item{x}{An object inheriting from \link{class} \dQuote{minimaxApprox}}
  \item{\dots}{Further arguments to \code{print}}
}

\details{
To print the raw \dQuote{minimaxApprox} object use \code{print.default}.
}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\seealso{
\code{\link{plot.minimaxApprox}}, \code{\link{coeff.minimaxApprox}}
}

\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
PP
print.default(PP)
}

\keyword{print}
\keyword{methods}
\keyword{NumericalMathematics}
