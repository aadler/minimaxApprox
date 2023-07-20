% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{minimaxApprox}
\alias{minimaxApprox}
\encoding{UTF-8}

\title{Minimax Approximation of Functions}

\description{
Calculates minimax approximations to functions. Polynomial approximation uses
the Remez (1962) algorithm. Rational approximation uses the Cody-Fraser-Hart
(Cody et al., 1968) version of the algorithm.
}

\usage{
minimaxApprox(fn, lower, upper, degree, relErr = FALSE, xi = NULL,
              opts = list())
}

\arguments{
    \item{fn}{function; A vectorized univariate function having \strong{x} as
    its first value.}
    \item{lower}{numeric; The lower bound of the approximation interval.}
    \item{upper}{numeric; The upper bound of the approximation interval.}
    \item{degree}{integer; Either a single value representing the requested
    degree for polynomial approximation or a vector of length 2 representing the
    requested degrees of the numerator and denominator for rational
    approximation.}
    \item{relErr}{logical; If \code{TRUE}, calculate the minimax approximation
    using \emph{relative} error. The default is \code{FALSE} which uses
    \emph{absolute} error.}
    \item{xi}{numeric; For rational approximation, a vector of initial points of
    the correct length---\eqn{\sum(\code{degree}) + 2}{sum(degree) + 2}. If
    missing, approximation will use the appropriate Chebyshev nodes. Polynomial
    approximation \strong{always} uses Chebyshev nodes and will ignore \code{xi}
    with a message.}
    \item{opts}{\link{list}; Configuration options including:
        \itemize{
            \item \code{maxiter}: integer; The maximum number of iterations to
            attempt convergence. Defaults to 100.
            \item \code{miniter}:  integer; The minimum number of iterations
            before allowing convergence. Defaults to 10.
            \item \code{conviter}: integer; The number of successive iterations
            with the same results allowed before assuming no further convergence
            is possible. Defaults to 10. Will overwrite \code{maxiter} if
            \code{conviter} is explicitly passed and is larger than
            \code{maxiter}.
            \item \code{showProgress}: logical; If \code{TRUE} will print error
            values at each iteration.
             \item \code{convRatio}: numeric; The convergence ratio tolerance.
             Defaults to \code{1+1e-9}. See \strong{Details}.
             \item \code{tol}: numeric; The absolute difference tolerance.
             Defaults to \code{1e-14}. See \strong{Details}.
        }
    }
}

\details{
The function implements the Remez algorithm using linear approximation, chiefly
as described by Cody (1968). Convergence is considered achieved when all three
of the following criteria are met:
    \enumerate{
        \item The observed error magnitudes are within tolerance of the expected
         error (Distance Test).
        \item The observed error magnitudes are within tolerance of each other
        (Magnitude Test).
        \item The observed error signs oscillate (Oscillating Test).
    }
\dQuote{Within tolerance} can be met in one of two ways:
    \enumerate{
        \item The difference between the absolute magnitudes is less than or
        equal to \code{tol}.
        \item The ratio between the larger and smaller is less than or equal to
        \code{convRatio}.
    }

For efficiency, the Distance Test is taken between the absolute value of the
largest observed error and the absolute value of the expected error. Similarly,
the Magnitude test is taken between the absolute value of the largest observed
error and the absolute value of the smallest observed error. Both the Magnitude
Test and the Distance Test can be passed by \strong{either} being within
\code{tol} or \code{convRatio} as described above.
}

\value{
\code{minimaxApprox} returns an object of \link{class} \code{"minimaxApprox"}
which inherits from the class \link{list}.

The generic accessor function \code{coef} will extract the numerator and
denominator vectors. There are also default \code{print} and \code{plot}
methods.

An object of class \code{"minimaxApprox"} is a list containing the following
components:

    \item{a}{The polynomial coefficients or the rational numerator
    coefficients.}
    \item{b}{The rational denominator coefficients. Missing for polynomial
    approximation.}
    \item{EE}{The absolute value of the expected error as calculated by the
    Remez algorithms.}
    \item{OE}{The absolute value of largest observed error between the function
    and the approximation at the extremal basis points.}
    \item{iterations}{The number of iterations of the algorithm. This does not
    include any iterations required to converge the error value in rational
    approximation.}
    \item{x}{The basis points at which the minimax error was achieved.}
    \item{Warning}{A logical flag indicating if any warnings were thrown.}
}

\references{
Remez, E. I. (1962) \emph{General computational methods of Chebyshev
approximation: The problems with linear real parameters}. US Atomic Energy
Commission, Division of Technical Information. AEC-tr-4491

Fraser W. and Hart J. F. (1962) \dQuote{On the computation of rational
approximations to continuous functions}, \emph{Communications of the ACM},
\bold{5}(7), 401--403, \doi{10.1145/368273.368578}

Cody, W. J. and Fraser W. and Hart J. F. (1968) \dQuote{Rational Chebyshev approximation using linear equations}, \emph{Numerische Mathematik}, \bold{12},
242--251, \doi{10.1007/BF02162506}
}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\note{
At present, the algorithms are implemented using machine double precision, which
means that the approximations are at best slightly worse. The value of the
approximations calculated by this package are for increasing efficiency in
coding and not---as of yet---for highly precise estimates.
}

\seealso{
\code{\link{minimaxEval}} for a convenience function to calculate approximation
values and \code{\link[Pade:Pade]{Pade}} for a function to calculate
\enc{Padé}{Pade} coefficients given suitable Taylor series coefficients.
}

\examples{
PP <- minimaxApprox(exp, 0, 1, 5)
PP
RR <- minimaxApprox(exp, 0, 1, c(2, 3))
RR
}

\keyword{optimize}
\keyword{NumericalMathematics}
