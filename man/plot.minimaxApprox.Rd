% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{plot.minimaxApprox}
\alias{plot.minimaxApprox}
\encoding{UTF-8}

\title{Plot errors from a \code{"minimaxApprox"} object}
\description{
Produces a plot of the error of the \code{"minimaxApprox"} object, highlighting
the basis points and error bounds.
}
\usage{\method{plot}{minimaxApprox}(x, y, \dots)}
\arguments{
    \item{x}{An object inheriting from \link{class} \code{"minimaxApprox"}.}
    \item{y}{Ignored. In call as required by \R in
    \href{https://cran.r-project.org/doc/manuals/R-exts.html#Generic-functions-and-methods}{Writing R Extensions:chapter 7}.}
    \item{\dots}{Further arguments to plot. Specifically to pass \code{ylim} to
    allow zooming in or out.}
}
\value{A plot of the approximation error.}
\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}
\seealso{
\code{\link{minimaxApprox}}
}
\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
plot(PP)
}

\keyword{hplot}
\keyword{methods}
\keyword{NumericalMathematics}
