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
    \code{"Monomial"} uses the standard \eqn{x^k} basis. \code{"Chebyshev"} uses
    the Chebyshev polynomials of the first kind, \eqn{T_k}. \code{"Barycentric"}
    for polynomials uses the barycentric-Remez algorithm of Pach\enc{ó}{o}n and
    Trefethen (2009), which represents the trial polynomial by its values at the
    reference rather than by coefficients in a fixed basis, avoiding the
    ill-conditioned linear solve entirely (see \strong{Barycentric Basis}
    below). For rational approximation, \code{"Barycentric"} uses the adaptive
    barycentric algorithm of Filip, Nakatsukasa, Trefethen, and Beckermann
    (2018) and supports \strong{absolute} error only; requesting it with
    \code{relErr = TRUE} is an error. The default is \code{"Chebyshev"}, and the
    parameter is case-insensitive and may be abbreviated.}
    \item{xi}{numeric; For rational approximation, a vector of initial points of
    the correct length---\eqn{\sum(\code{degree}) + 2}{sum(degree) + 2}. If
    missing, the approximation will use the appropriate Chebyshev nodes,
    except that rational approximation in the barycentric basis initializes
    from the error extrema of an AAA approximant (see \strong{Barycentric
    Basis} below); a supplied \code{xi} of the correct length overrides that
    initialization. Polynomial approximation \strong{always} uses Chebyshev
    nodes and will ignore \code{xi} with a message.}
    \item{opts}{\link{list}; Configuration options including:
        \itemize{
            \item \code{maxiter}: integer; The maximum number of iterations to
            attempt convergence. Defaults to 100.
            \item \code{miniter}:  integer; The minimum number of iterations
            before allowing convergence. Defaults to 10.
            \item \code{conviter}: integer; The number of successive iterations
            with the same results allowed before assuming no further convergence
            is possible. Defaults to 30. Will overwrite \code{maxiter} and
            \code{miniter} if \code{conviter} is explicitly passed and is larger
            than either one.
            \item \code{showProgress}: logical; If \code{TRUE} will print error
            values at each iteration.
            \item \code{convrat}: numeric; The convergence ratio tolerance.
            Defaults to \eqn{1 + 1 \times 10^{-9}}{1+1e-9}. See
            \strong{Details}.
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

When the error values remain within tolerance of each other over \code{conviter}
iterations, the algorithm will stop, as it is expected that no further precision
will be gained by continued iterations.
}
\subsection{Polynomial Evaluation}{
Monomial polynomials are evaluated using the Compensated Horner Scheme of
Langlois et al. (2006) to enhance both stability and precision. Chebyshev
polynomials are evaluated on an internally affine-mapped domain (see
\strong{Mapped Chebyshev Basis} below). There may be cases where the algorithm
will fail using the monomial basis but succeed using Chebyshev polynomials and
vice versa. The default is to use the Chebyshev polynomials.
}
\subsection{Mapped Chebyshev Basis}{
Chebyshev polynomials \eqn{T_k} are only well-conditioned on \eqn{[-1, 1]}; off
that interval, conditioning degrades severely and grows with both distance from
\eqn{[-1, 1]} and degree (for example, degree 10 on \eqn{[5, 6]} has a condition
number around \eqn{10^{24}} unmapped versus about \eqn{1.4} mapped). To avoid
this, the Chebyshev basis internally maps the requested range \eqn{[l, u]} to
\eqn{[-1, 1]} via \eqn{z=\frac{2x - (l + u)}{u-l}}{z = (2x - (l + u)) / (u - l)}
and evaluates \eqn{T_k} on \eqn{z}. The returned \code{a} (and \code{b})
are therefore coefficients of \eqn{T_k(z)}, \strong{not} \eqn{T_k(x)}, whenever
the fitted range is not \eqn{[-1, 1]}; on \eqn{[-1, 1]} itself the map is the
identity and \code{a}/ \code{b} are unchanged from prior releases.

\code{aMono}/\code{bMono} remain coefficients in the natural, raw \eqn{x} in
every case, and are the stable way to consume a Chebyshev-basis result off
\eqn{[-1, 1]} without needing to know about the internal map. They are obtained
by converting the mapped-Chebyshev coefficients to monomials in \eqn{z} as
before and then composing that polynomial with the affine map above to recover
monomial coefficients in \eqn{x}.

This composition step is itself a monomial-basis operation and inherits
monomial ill-conditioning at high degree, the same way \code{aMono} always has.
On ranges far from \eqn{[-1, 1]} at higher degrees, the composition's own error
can grow to dominate the underlying fit's accuracy \eqn{-} for example, a
degree-15 fit on \eqn{[10, 20]} can have a well-conditioned Chebyshev-basis
error near \eqn{10^{-7}} while the composed \code{aMono} coefficients disagree
with the Chebyshev-basis evaluation by closer to \eqn{10^{-6}} on the same
range. When degree is high and the range is both wide and far from
\eqn{[-1, 1]}, prefer evaluating with the Chebyshev basis directly (the default
for \code{\link{minimaxEval}}/\code{\link{minimaxErr}}) over consuming
\code{aMono} directly, if accuracy at that level matters.
}
\subsection{Barycentric Basis}{
\subsection{Polynomial Approximation}{
Selecting \code{basis = "Barycentric"} for \strong{polynomial} approximation
uses the barycentric-Remez algorithm of Pach\enc{ó}{o}n and Trefethen (2009).
Rather than solving an \eqn{(n + 2) \times (n + 2)} linear system in a fixed
basis at each iteration---the source of the conditioning failures that cap
the achievable degree in the classical bases---it represents the trial
polynomial by its \emph{values} at the reference points, with barycentric
weights, and computes the leveled error \eqn{h} in closed form. Evaluation uses
the second barycentric formula, which is forward-stable for well-distributed
nodes such as these (Higham, 2004). Weights are
assembled with capacity scaling so that wide or far-off-center intervals and
high degrees do not over- or underflow. As a result the barycentric basis
typically converges at degrees where the classical bases fail with a
singular-matrix error, and is well-conditioned off \eqn{[-1, 1]} without a
separate mapping step.

Functions that are exactly representable, or resolved to machine precision at
the requested degree (the cases described in \strong{Machine-Precision-Resolved
Functions} below), are handled \emph{natively}: the closed-form leveled error is
\eqn{0} and the trial polynomial interpolates \code{fn}, so no singular-solve
rescue is needed. Relative error is supported via the analogous closed form;
however, if \code{fn} is \strong{exactly} \eqn{0} at a reference point
(typically an endpoint of \eqn{[l, u]}, which is always a reference node), the
relative error there is \eqn{0/0} and the relative minimax does not exist---this
is reported as a clear error rather than a silent or crashed result. Use
absolute error, or an interval whose endpoints are not zeros of \code{fn}.

The returned object stores the barycentric representation (reference nodes,
weights, and nodal values) as its primary form; \code{a}, \code{aMono}, and the
convenience methods are derived from it. Because the coefficient recovery is the
one step that can lose accuracy at high degree, a conversion residual is
reported (see \code{convResid} in \strong{Value}); for functions with poles near
the approximation interval this residual can be noticeably larger at some
intermediate degrees (for example, a degree-20 fit of \eqn{1/(1 + x^2)} on
\eqn{[-1, 1]} shows a residual near \eqn{10^{-8}} while degrees 10 and 50 are
near \eqn{10^{-15}}). The \code{tailtol} degree-\code{n + 1} retry and the
\code{ztol} coefficient-zeroing test do not apply to the barycentric basis,
which has no coefficient solve to restart or prune.
}
\subsection{Rational Approximation}{
\strong{Rational approximation} in the barycentric basis implements the
adaptive-barycentric minimax algorithm of Filip, Nakatsukasa, Trefethen, and
Beckermann (2018): the trial rational is represented as a quotient of
barycentric forms over \emph{support points} chosen as every other reference
point, and the leveled error at each iteration is an eigenvalue of a small
\emph{symmetric} matrix, so there is again no ill-conditioned linear solve.
The initial reference is taken from the error extrema of an AAA approximant
(Nakatsukasa, S\enc{è}{e}te, and Trefethen, 2018) computed on a dense sample
grid, falling back to Chebyshev points when the AAA error curve does not
supply enough alternations. The returned \code{a}/\code{b} (and
\code{aMono}/\code{bMono}) are recovered from the barycentric representation
with the denominator normalized to a constant term of \eqn{1}, matching the
classical rational path, and \code{convResid} is a named length-2 vector
reporting the relative conversion residual of the numerator (\code{"a"}) and
denominator (\code{"b"}) separately.

The rational barycentric implementation is a deliberate subset of the
published algorithm, with the following limitations. Only \strong{absolute}
error is supported; \code{relErr = TRUE} is rejected with an error.
Degenerate or \emph{defective} problems---for example, requesting even
numerator and denominator degrees for an even function, where the true best
approximation has fewer effective parameters than the requested type---are
\strong{detected and reported as errors} rather than repaired by degree
reduction as in the full published algorithm; the error message identifies
the detected cause (a rank-deficient basis matrix, no pole-free solution of
the leveled-error equations, or a collapsed barycentric weight). A trial
denominator with a zero inside the approximation interval is likewise
reported as an error, mirroring the classical rational path. Finally, for
\emph{non-normal} requests (degrees for which the true best approximation
equioscillates on more than \eqn{m + n + 2} points, as happens for even or odd
functions at generic degrees) the fixed-size exchange can converge to a leveled
fixed point that is not the global minimax; the returned approximation is then
valid and its expected error is a lower bound on the true minimax error, and a
warning to that effect is raised whenever the approximation's dense-grid error
exceeds its leveled error. The same reference-local caveat applies, without the
automatic check, to results within a few orders of magnitude of machine
precision on any basis.
}
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
\subsection{Machine-Precision-Resolved Functions}{
When the requested degree is high enough that \code{fn} is resolved to (or
extremely near) machine double precision---including functions that are
\emph{exactly} representable, such as a quadratic requested at degree 2---there
may be no minimax problem left to solve: any polynomial passing through a
correctly-chosen reference already matches \code{fn} to the limits of floating
point, so the Remez iteration has nothing left to optimize. This is detected
either as a singular linear solve at both degree \code{n} and degree
\code{n + 1} (see above), or, on some platforms, as the iteration running to
\code{maxiter} or stalling \dQuote{unchanging} while the observed error is
already at the machine-precision floor. In that situation, \code{minimaxApprox}
returns the degree-\code{n} polynomial that \emph{interpolates} \code{fn}
through the Chebyshev reference, rather than a Remez-optimized polynomial,
together with a warning stating that the result is not technically a Remez
result. This interpolant equals the true minimax polynomial when \code{fn} is
exactly representable, and is indistinguishable from it (to within floating
point) otherwise; \code{ExpErr} and \code{ObsErr} are equal in this case, since
there is no separate leveled error to report. If the observed error is
\strong{not} at the machine floor, the singularity or stall has a genuine cause
unrelated to representability, and the original error or convergence warning is
raised as before---this fallback never masks real non-convergence.
}
\subsection{Close-to-Zero Tolerance}{
For each step of the algorithms' iterations, the contribution of the found
coefficient to the total sum (as measured in the above section) is compared to
the \code{ztol} option. When less than or equal to \code{ztol}, that coefficient
is set to 0. Setting \code{ztol} to \code{NULL} skips the test completely. For
intervals near or containing zero, setting this option to anything other than
\code{NULL} may result in either non-convergence or poor results. It is
recommended to keep it as \code{NULL}, although there are edge cases where
it may allow convergence where a standard call may fail.
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
    Chebyshev polynomials, these are the coefficients for \eqn{T_k} of the
    internally mapped variable \eqn{z} (see \strong{Mapped Chebyshev Basis}
    above) whenever the fitted range is not \eqn{[-1, 1]}; on \eqn{[-1, 1]}
    itself, \eqn{z = x} and these are unchanged from prior releases. When
    using the barycentric basis, these are the mapped-Chebyshev coefficients
    recovered from the barycentric representation (see \strong{Barycentric
    Basis} above). When using monomials, these are the coefficients for
    \eqn{x^k}.}
    \item{b}{The rational denominator coefficients. Same semantics as
    \code{a} above. Missing for polynomial approximation.}
    \item{aMono}{When using Chebyshev or barycentric polynomials, these are the
    polynomial or rational numerator coefficients for monomial expansion in
    \eqn{x^k} (natural, raw \eqn{x}), regardless of the internal map. See
    \strong{Mapped Chebyshev Basis} above for an accuracy caveat at high degree
    on wide, far-off-\eqn{[-1,1]} ranges; the barycentric basis reports the size
    of this conversion loss directly in \code{convResid}. Missing for
    monomial-based approximation.}
    \item{bMono}{When using Chebyshev polynomials, these are the rational
    denominator coefficients for monomial expansion in \eqn{x^k} (natural, raw
    \eqn{x}; same caveat as \code{aMono}). Missing for both polynomial and
    monomial-based rational approximation.}
    \item{ExpErr}{The absolute value of the expected error as calculated by the
    Remez algorithms. For a machine-precision-resolved result (see
    \strong{Machine-Precision-Resolved Functions} above), there is no separate
    Remez-computed expected error; \code{ExpErr} equals \code{ObsErr} in that
    case.}
    \item{ObsErr}{The absolute value of largest observed error between the
    function and the approximation at the extremal points.}
    \item{iterations}{The number of iterations of the algorithm. This does not
    include any iterations required to converge the error value in rational
    approximation.}
    \item{Extrema}{The extrema at which the minimax error was achieved.}
    \item{Warning}{A logical flag indicating if any warnings were thrown.}
    \item{convResid}{Barycentric basis only. The conversion residual: the
    maximum absolute difference between the (accurate) barycentric evaluation of
    the fitted polynomial and its reconstruction from the recovered coefficients,
    measured on a dense probe grid. A diagnostic of how much accuracy the
    coefficient recovery lost; see \strong{Barycentric Basis} above. For
    rational approximation it is a named length-2 vector giving the
    \emph{relative} residual of the numerator (\code{"a"}) and denominator
    (\code{"b"}) polynomials separately.}

The object also contains the following attributes:

    \item{type}{"Rational" or "Polynomial".}
    \item{basis}{"Monomial", "Chebyshev", or "Barycentric".}
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

Pach\enc{ó}{o}n, R. and Trefethen, L. N. (2009) \dQuote{Barycentric-Remez
algorithms for best polynomial approximation in the chebfun system}, \emph{BIT
Numerical Mathematics}, \bold{49}(4), 721--741, \doi{10.1007/s10543-009-0240-1}

Berrut, J.-P. and Trefethen, L. N. (2004) \dQuote{Barycentric Lagrange
interpolation}, \emph{SIAM Review}, \bold{46}(3), 501--517,
\doi{10.1137/S0036144502417715}

Higham, N. J. (2004) \dQuote{The numerical stability of barycentric Lagrange
interpolation}, \emph{IMA Journal of Numerical Analysis}, \bold{24}(4),
547--556, \doi{10.1093/imanum/24.4.547}

Filip, S.-I. and Nakatsukasa, Y. and Trefethen, L. N. and Beckermann, B.
(2018) \dQuote{Rational minimax approximation via adaptive barycentric
representations}, \emph{SIAM Journal on Scientific Computing}, \bold{40}(4),
A2427--A2455, \doi{10.1137/17M1132409}

Nakatsukasa, Y. and S\enc{è}{e}te, O. and Trefethen, L. N. (2018) \dQuote{The
AAA algorithm for rational approximation}, \emph{SIAM Journal on Scientific
Computing}, \bold{40}(3), A1494--A1522, \doi{10.1137/16M1106122}
}

\author{
Avraham Adler \email{Avraham.Adler@gmail.com}
}

\note{
At present, the algorithms are implemented using machine double precision, which
means that the approximations are at best slightly worse. Research proceeds on
more precise, stable, and efficient implementations. So long as the package
remains in an experimental state---noted by a 0 major version---the API may
change at any time.
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

minimaxApprox(exp, -1, 1, 14, basis = "b")              # Barycentric: converges
                                                        # at degrees where the
                                                        # classical bases fail
}

\keyword{optimize}
\keyword{NumericalMathematics}
