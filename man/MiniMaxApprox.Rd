\name{MiniMaxApprox}
\alias{MiniMaxApprox}
\alias{print.MiniMaxApprox}
\alias{plot.MiniMaxApprox}
\alias{coef.MiniMaxApprox}

\title{Minimax Approximation of Functions}
\description{
Calculates minimax approximations to functions. Polynomial approximations use
the second Remez algorithm. Rational approximation uses a form of the
Fraser-Hart algorithm.
}
\usage{
MiniMaxApprox(fn, lower, upper, degree, errType = "abs", xi = NULL, opts = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fn}{
%%     ~~Describe \code{fn} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{degree}{
%%     ~~Describe \code{degree} here~~
}
  \item{errType}{
%%     ~~Describe \code{errType} here~~
}
  \item{xi}{
%%     ~~Describe \code{xi} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
PP <- MiniMaxApprox(exp, 0, 1, 5)
RR <- MiniMaxApprox(exp, 0, 1, c(2, 3))
}

\keyword{NumericalMathematics}
