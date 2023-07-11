\name{minimaxApprox}
\alias{minimaxApprox}
\alias{print.minimaxApprox}
\alias{plot.minimaxApprox}
\alias{coef.minimaxApprox}

\title{Minimax Approximation of Functions}
\description{
Calculates minimax approximations to functions. Polynomial approximations use
the second Remez algorithm. Rational approximation uses a form of the
Fraser-Hart algorithm.
}
\usage{
minimaxApprox(fn, lower, upper, degree, errType = "abs", xi = NULL, opts = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
  \strong{function}: A univariate function having \strong{x} as its first value
  }
  \item{lower}{\strong{numeric}: The lower bound of the approximation interval}
  \item{upper}{\strong{numeric}: The upper bound of the approximation interval}
  \item{degree}{
  Either a single \strong{integer} representing the requested degree for
  polynomial approximation or an \strong{integer} vector of length 2
  representing the requested degrees of the numerator and denominator for
  rational approximation.
}
  \item{errType}{\strong{character}: The type of error to be minimized over the
  interval. Either '\strong{abs}'olute or '\strong{rel}'ative.
}
  \item{xi}{
  \strong{numeric}: OPTIONAL: A vector of initial points of the correct length,
  \eqn{\sum(\code{degree}) + 2}{sum(degree) + 2}, for use in rational
  approximation. Polynomial approximation will ignore \code{xi} with a message
  as it always uses the appropriate Chebyeshev nodes as a starting point. If
  left out, rational approximation will also use the Chebyeshev nodes.
}
  \item{opts}{
%%     ~~Describe \code{opts} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
RR <- minimaxApprox(exp, 0, 1, c(2, 3))
}

\keyword{NumericalMathematics}
