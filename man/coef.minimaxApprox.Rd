\name{coef.minimaxApprox}
\alias{coef.minimaxApprox}
\title{Extract coefficiencts from a \dQuote{minimaxApprox} object}

\description{
Extracts the numerator and denominator vectors from a \dQuote{minimaxApprox}
object.
}

\usage{
\method{coef}{minimaxApprox}(object, \dots)
}

\arguments{
  \item{object}{An object inheriting from \link{class} \dQuote{minimaxApprox}.
}
  \item{\dots}{Other arguments.}
}

\value{Coefficients extracted from the \dQuote{minimaxApprox} object. A
\link{list} containing:

  \item{a}{The polynomial coefficents or the rational numerator coefficents}
  \item{b}{The rational denominator coefficents. Missing for polynomial
  approximation}
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
identical(unlist(coef(PP), use.names = FALSE), PP$a)
RR <- minimaxApprox(exp, 0, 1, c(2, 3))
coef(RR)
identical(coef(RR), list(a = RR$a, b = RR$b))
}

\keyword{methods}
\keyword{NumericalMathematics}
