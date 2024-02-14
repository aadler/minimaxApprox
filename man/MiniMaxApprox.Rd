% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{minimaxApprox}
\alias{minimaxApprox}
\encoding{UTF-8}

\title{Minimax Approximation of Functions}

\description{
Calculates minimax approximations to functions. Polynomial approximation uses
the Remez (1962) algorithm. Rational approximation uses the Cody-Fraser-Hart
(Cody et al., 1968) version of the algorithm. When using monomials as the
polynomial basis, the Compensated Horner Scheme of Langlois et al. (2006) is
used.
}

\usage{
minimaxApprox(fn, lower, upper, degree, relErr = FALSE, basis ="Chebyshev",
              xi = NULL, opts = list())
}

\arguments{
    \item{fn}{function; A vectorized univariate function having \code{x} as its
    first argument. This could be a built-in \R function, a predefined function,
    or an anonymous function defined in the call; see \strong{Examples}.}
    \item{lower}{numeric; The lower bound of the approximation interval.}
    \item{upper}{numeric; The upper bound of the approximation interval.}
    \item{degree}{integer; Either a single value representing the requested
    degree for polynomial approximation or a vector of length 2 representing the
    requested degrees of the numerator and denominator for rational
    approximation.}
    \item{relErr}{logical; If \code{TRUE}, calculate the minimax approximation
    using \emph{relative} error. The default is \code{FALSE} which uses
    \emph{absolute} error.}
    \item{basis}{character; Which polynomial basis to use in the analysis.
    "Monomial" uses the standard \eqn{x^k} basis. "Chebyshev" uses the Chebyshev
    polynomials of the first kind, \eqn{T_k}. The default is "Chebyshev", and
    the parameter is case-insensitive and may be abbreviated.}
    \item{xi}{numeric; For rational approximation, a vector of initial points of
    the correct length---\eqn{\sum(\code{degree}) + 2}{sum(degree) + 2}. If
    missing, the approximation will use the appropriate Chebyshev nodes.
    Polynomial approximation \strong{always} uses Chebyshev nodes and will
    ignore \code{xi} with a message.}
    \item{opts}{\link{list}; Configuration options including:
        \itemize{
            \item \code{maxiter}: integer; The maximum number of iterations to
            attempt convergence. Defaults to 100.
            \item \code{miniter}:  integer; The minimum number of iterations
            before allowing convergence. Defaults to 10.
            \item \code{conviter}: integer; The number of successive iterations
            with the same results allowed before assuming no further convergence
            is possible. Defaults to 30. Will overwrite \code{maxiter} if
            \code{conviter} is explicitly passed and is larger than
            \code{maxiter}.
            \item \code{showProgress}: logical; If \code{TRUE} will print error
            values at each iteration.
            \item \code{convrat}: numeric; The convergence ratio tolerance.
            Defaults to \eqn{1 + 1 \times 10^{-9}}{1+1e-9}. See \strong{Details}.
            \item \code{tol}: numeric; The absolute difference tolerance.
            Defaults to \eqn{1 \times 10^{-14}}{1e-14}. See \strong{Details}.
            \item \code{tailtol}: numeric; The tolerance of the coefficient of
            the largest power of \code{x} to be ignored when performing the
            polynomial approximation a second time. Defaults to the smaller of
            \eqn{1 \times 10^{-10}}{1e-10} or
            \eqn{\frac{\code{upper} - \code{lower}}{10^6}}{(\code{upper} -
            \code{upper}) / 1e6}. Set to \code{NULL} to skip the
            \code{degree + 1} check completely. See \strong{Details}.
            \item \code{ztol}: numeric; The tolerance for each polynomial or
            rational numerator or denominator coefficient's contribution to
            \strong{not} to be set to 0. Similar to polynomial \code{tailtol}
            but applied at each step of the algorithm. Defaults to
            \code{NULL} which leaves all coefficients as they are regardless
            of magnitude. See \strong{Details}.
        }
    }
}

\details{
\subsection{Convergence}{
The function implements the Remez algorithm using linear approximation, chiefly
as described by Cody et al. (1968). Convergence is considered achieved when all
three of the following criteria are met:
    \enumerate{
        \item The observed error magnitudes are within tolerance of the expected
         error---the \strong{Distance Test}.
        \item The observed error magnitudes are within tolerance of each
        other---the \strong{Magnitude Test}.
        \item The observed error signs oscillate---the
        \strong{Oscillation Test}.
    }
\dQuote{Within tolerance} can be met in one of two ways:
    \enumerate{
        \item \strong{Difference}: The difference between the absolute
        magnitudes is less than or equal to \code{tol}.
        \item \strong{Ratio}: The ratio between the absolute magnitudes of the
        larger and smaller is less than or equal to \code{convrat}.
    }

For efficiency, the \strong{Distance Test} is taken between the absolute value
of the largest observed error and the absolute value of the expected error.
Similarly, the \strong{Magnitude Test} is taken between the absolute value of
the largest observed error and the absolute value of the smallest observed
error. Both tests can be passed by \strong{either} being within \code{tol} or
\code{convrat} as described above. However, when the \strong{Difference} test
returns values less than machine precision, it is ignored in favor of the
\strong{Ratio} test.
}
\subsection{Polynomial Evaluation}{
Monomial polynomials are evaluated using the Compensated Horner Scheme of
Langlois et al. (2006) to enhance both stability and precision. Chebyshev
polynomials are evaluated normally. There may be cases where the algorithm will
fail using the monomial basis but succeed using Chebyshev polynomials. The
default is to use the Chebyshev polynomials.
}
\subsection{Polynomial Algorithm \dQuote{Singular Error} Response}{
When too high of a degree is requested for the tolerance of the algorithm, it
often fails with a singular matrix error. In this case, for the
\emph{polynomial} version, the algorithm will try looking for an approximation
of degree \code{n + 1}. If it finds one, \strong{and} the contribution of that
coefficient to the approximation is \eqn{\le} \code{tailtol}, it will ignore
that coefficient and return the resulting degree \code{n} polynomial, as the
largest coefficient is effectively 0. The contribution is measured by
multiplying that coefficient by the endpoint with the larger absolute magnitude
raised to the \code{n + 1} power. This is done to prevent errors in cases where
a very small coefficient is found on a range with very large absolute values and
the resulting contribution to the approximation is \strong{not}
\emph{de minimis}. Setting \code{tailtol} to \code{NULL} will skip the
\code{n + 1} test completely.
}
\subsection{Close-to-Zero Tolerance}{
For each step of the algorithms' iterations, the contribution of the found
coefficient to the total sum (as measured in the above section) is compared to
the \code{ztol} option. When less than or equal to \code{ztol}, that coefficient
is set to 0. Setting \code{ztol} to \code{NULL} skips the test completely. For
intervals near or containing zero, setting this option to anything other than
\code{NULL} may result in either non-convergence or poor results. It is
recommended to keep it as \code{NULL}.
}
}
\value{
\code{minimaxApprox} returns an object of \link{class} \code{"minimaxApprox"}
which inherits from the class \link{list}.

The generic accessor function \code{coef} will extract the numerator and
denominator vectors. There are also default \code{print} and \code{plot}
methods.

An object of class \code{"minimaxApprox"} is a list containing the following
components:

    \item{a}{The polynomial or rational numerator coefficients. When using
    Chebyshev polynomials, these are the coefficients for \eqn{T_k}.}
    \item{b}{The rational denominator coefficients. When using Chebyshev
    polynomials, these are the coefficients for \eqn{T_k}. Missing for
    polynomial approximation.}
    \item{aMono}{When using Chebyshev polynomials, these are the polynomial or
    rational numerator coefficients for monomial expansion in \eqn{x^k}. Missing
    for monomial-based approximation.}
    \item{bMono}{When using Chebyshev polynomials, these are the rational
    denominator coefficients for monomial expansion in \eqn{x^k}. Missing for
    both polynomial and monomial-based rational approximation.}
    \item{ExpErr}{The absolute value of the expected error as calculated by the
    Remez algorithms.}
    \item{ObsErr}{The absolute value of largest observed error between the
    function and the approximation at the extremal points.}
    \item{iterations}{The number of iterations of the algorithm. This does not
    include any iterations required to converge the error value in rational
    approximation.}
    \item{Extrema}{The extrema at which the minimax error was achieved.}
    \item{Warning}{A logical flag indicating if any warnings were thrown.}

The object also contains the following attributes:

    \item{type}{"Rational" or "Polynomial".}
    \item{basis}{"Monomial" or "Chebyshev".}
    \item{func}{The function being approximated.}
    \item{range}{The range on which the function is being approximated.}
    \item{relErr}{A logical indicating that relative error was used. If
    \code{FALSE}, then absolute error was used.}
    \item{tol}{The tolerance used for the \strong{Distance Test}.}
    \item{convrat}{The tolerance used for the \strong{Magnitude Test}.}
}

\references{
Remez, E. I. (1962) \emph{General computational methods of Chebyshev
approximation: The problems with linear real parameters}. US Atomic Energy
Commission, Division of Technical Information. AEC-tr-4491

Fraser W. and Hart J. F. (1962) \dQuote{On the computation of rational
approximations to continuous functions}, \emph{Communications of the ACM},
\bold{5}(7), 401--403, \doi{10.1145/368273.368578}

Cody, W. J. and Fraser W. and Hart J. F. (1968) \dQuote{Rational Chebyshev
approximation using linear equations}, \emph{Numerische Mathematik}, \bold{12},
242--251, \doi{10.1007/BF02162506}

Langlois, P. and Graillat, S. and Louvet, N. (2006) \dQuote{Compensated Horner
Scheme}, in \emph{Algebraic and Numerical Algorithms and Computer-assisted
Proofs.} Dagstuhl Seminar Proceedings, \bold{5391},
\doi{10.4230/DagSemProc.05391.3}
}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\note{
At present, the algorithms are implemented using machine double precision, which
means that the approximations are at best slightly worse. Using long doubles for
internal calculations for both monomials and Chebyshev polynomials did not
result in any material gains. Research proceeds on more precise, stable, and
efficient implementations. So long as the package remains in an experimental
state---noted by a 0 major version---the API may change at any time.
}

\seealso{\code{\link{minimaxEval}}, \code{\link{minimaxErr}}}

\examples{
minimaxApprox(exp, 0, 1, 5)                              # Built-in & polynomial

fn <- function(x) sin(x) ^ 2 + cosh(x)                   # Pre-defined
minimaxApprox(fn, 0, 1, c(2, 3), basis = "m")            # Rational

minimaxApprox(function(x) x ^ 3 / sin(x), 0.7, 1.6, 6L)  # Anonymous

fn <- function(x) besselJ(x, nu = 0)                     # More than one input
b0 <- 0.893576966279167522                               # Zero of besselY
minimaxApprox(fn, 0, b0, c(3L, 3L))                      # Cf. DLMF 3.11.19
}

\keyword{optimize}
\keyword{NumericalMathematics}
