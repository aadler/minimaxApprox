\name{minimaxApprox}
\encoding{UTF-8}
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

\method{print}{minimaxApprox}(x, \dots)
\method{plot}{minimaxApprox}(x, y, \dots)
\method{coef}{minimaxApprox}(object, \dots)
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
  \item{x}{For \code{print} and \code{plot}, a \code{minimaxApprox} object.}
  \item{y}{
  Ignored by \code{plot}; in call as required by
  \href{https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods}{Writing R Extensions:chapter 7}.
}
  \item{object}{For \code{coef}, a \code{minimaxApprox} object.}
  \item{\dots}{
  Further arguments to be passed to other methods. Specifically to pass
  \code{ylim} to plot to allow zooming in or out.
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
Remez, E. I. (1962) \emph{General computational methods of Chebyshev
approximation: The problems with linear real parameters}. US Atomic Energy
Commission, Division of Technical Information. AEC-tr-4491

Fraser W. and Hart J. F. (1962) \dQuote{On the computation of rational
approximations to continuous functions}, \emph{Communications of the ACM},
\bold{5}(7), 401--403, \doi{10.1145/368273.368578}

Cody, W. J. (1968) \dQuote{Rational Chebyshev approximation using linear
equations}, \emph{Numerische Mathematik}, \bold{12}, 242--251,
\doi{10.1007/BF02162506}

}
\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning}{
As of this version, the algorithms are implemented using machine double
precision, which means that the approximations are at best slightly worse. The
value of the approximations calculated by this package are for increasing
efficiency in coding and not, as of yet, for highly precise estimates.
}

\seealso{
\code{\link{minimaxEval}} for a conveniece function to calculate approximation
values and \code{\link[Pade:Pade]{Pade}} for a function to calculate
\enc{Padé}{Pade} coefficients given suitable Taylor series coefficients.
}
\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
PP
RR <- minimaxApprox(exp, 0, 1, c(2, 3))
RR
}

\keyword{NumericalMathematics}
