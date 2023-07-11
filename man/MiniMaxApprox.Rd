\name{minimaxApprox}
\alias{minimaxApprox}
\alias{print.minimaxApprox}
\alias{plot.minimaxApprox}
\alias{coef.minimaxApprox}

\title{Minimax Approximation of Functions}
\description{
Calculates minimax approximations to functions. Polynomial approximation uses
the Remez algorithm. Rational approximation uses a form of the Fraser-Hart
version of the Remez algorithm.
}
\usage{
minimaxApprox(fn, lower, upper, degree, errType = "abs", xi = NULL, opts = list())

\method{print}{minmaxApprox}(x, \dots)
\method{plot}{minmaxApprox}(x, ylim = NULL, \dots)
\method{coef}{minmaxApprox}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{function; A univariate vectorized function having \strong{x} as its
  first value.
  }
  \item{lower}{numeric; The lower bound of the approximation interval.}
  \item{upper}{numeric; The upper bound of the approximation interval.}
  \item{degree}{integer; Either a single value representing the requested degree
  for polynomial approximation or a vector of length 2 representing the
  requested degrees of the numerator and denominator for rational approximation.
}
  \item{errType}{character; The type of error to be minimized over the interval.
  Either '\strong{abs}'olute (default) or '\strong{rel}'ative.
}
  \item{xi}{numeric; For rational approximation, a vector of initial points of
  the correct length---\eqn{\sum(\code{degree}) + 2}{sum(degree) + 2}. If
  missing, approximation will use the appropriate Chebyeshev nodes. Polynomial
  approximation \emph{always} uses Chebyeshev nodes and will ignore \code{xi}
  with a message.
}
  \item{opts}{\link{list}; Configuration options including:
  \itemize{
  \item \code{maxiter}: integer; The maximum number of iterations to attempt
  convergence. Defaults to 100.
  \item \code{miniter}:  integer; The minimum number of iterations before
  allowing convergence. Defaults to 10.
  \item \code{conviter}: integer; For rational approximation only, the number of
  successive iteration with the same results allowed before assuming no further
  convergence is possible. Defaults to 10L. Will overwrite \code{maxiter} if
  passed and larger.
  \item \code{showProgress}: logical; If \code{TRUE} will print error values at
  each iteration.
  \item \code{convRatio}: numeric; The convergence ratio tolerance. Defaults to
  \code{1+1e-9}. See \strong{Details}.
  \item \code{tol}: numeric; The absolute difference tolerance. Defaults to
  \code{1e-14}. See \strong{Details}.
  }
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
PP
RR <- minimaxApprox(exp, 0, 1, c(2, 3))
RR
}

\keyword{NumericalMathematics}
