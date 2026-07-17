% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{NEWS}
\alias{NEWS}

\title{NEWS file for the minimaxApprox package}

\section{Developmental Note: Version 0.0.0+}{
  \itemize{
    \item So long as the package remains in its experimental development
    state---noted by a 0 major version---the API may change without incrementing
    the major version. Please read these development notes carefully. Breaking
    changes will be prefaced by \dQuote{\strong{Breaking:}}.
  }
}

\section{Version 0.6.0 (2026-10-01)}{
  \subsection{Fixed}{
    \itemize{
      \item Fixed the Remez exchange silently converging to a
      \emph{reference-local} fixed point that is not the global minimax,
      returning a materially suboptimal approximation with no warning (for
      example, \code{atan} on \eqn{[0, 3]} at degree 3 reported an expected
      error 1.8 times smaller than the returned approximation's true maximum
      error). Two structural defects were responsible: the root search only
      ever examined the intervals \emph{between} consecutive current reference
      points, so a sign change between an endpoint and the outermost reference
      point was invisible; and the maximize/minimize schedule was seeded at
      the lower endpoint and flipped mechanically, so an error curve with one
      more oscillation than the schedule assumed desynced the exchange into
      returning roots instead of extrema. Both are fixed by the exchange
      redesign described under \dQuote{Changed}.

      \item Fixed a crash (\dQuote{segfault from C stack overflow}) evaluating a
      Chebyshev-basis approximation at several million or more points at once.

      \item Added a clear error, instead of undefined behavior, when an internal
      Chebyshev evaluation would require more than 2^31 matrix cells.

      \item Fixed \code{lower > upper} (or \code{lower == upper}) being silently
      accepted and producing a badly suboptimal approximation with no warning;
      this now raises a clear error suggesting the arguments be swapped.

      \item Fixed the Remez-iteration stagnation check (\code{isUnchanging})
      reporting convergence as \dQuote{unchanging} when errors were still
      improving rapidly between iterations. The ratio test was checked in only
      one direction. Now, it requires the ratio to be close to 1 in either
      direction before treating the solution as stalled.

      \item Fixed \code{opts$ztol} using the monomial-basis coefficient scale
      even when the Chebyshev basis was requested, which could incorrectly zero
      out genuine small-but-relevant Chebyshev coefficients on ranges
      intersecting \eqn{[-1, 1]}.

      \item Fixed the degree-\eqn{n+1} restart's \dQuote{uppermost coefficient
      is effectively 0} check silently accepting any negative top coefficient
      regardless of magnitude (a missing \code{abs()}), which could return a
      polynomial one degree lower than warranted.

      \item Fixed severe numerical conditioning loss when fitting the
      Chebyshev basis on ranges far from \eqn{[-1, 1]} (e.g. degree 10 on
      \eqn{[5, 6]} had a condition number around \eqn{10^{24}}), which could
      cap the achievable accuracy or force convergence warnings well short of
      the true minimax error (that \eqn{[5, 6]} case: observed error improved
      from \eqn{\sim 1.8 \times 10^{-10}} to \eqn{\sim 3 \times 10^{-12}},
      matching the true value). See the \strong{Breaking} change below.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item The Remez exchange was redesigned to be driven by the error curve
      itself rather than by the current reference's structure. Roots of the
      error equation are now located over the \emph{whole} approximation
      interval from an oversampled grid; each sign-region's search direction
      is taken from the error's own sign there instead of a mechanically
      alternated schedule; the interval endpoints are always retained as
      exchange candidates; and when the curve offers more alternating extrema
      than the reference holds, the new reference is the alternating window
      that contains the global extremum and maximizes its smallest error (the
      classical conditions for the leveled error to increase strictly toward
      the true minimax). Converged results for healthy inputs agree with
      0.5.0 at tolerance level but are \strong{not bitwise identical}, and
      inputs that previously exited through a singular-matrix error or a
      stall may now exit through a different---equally loud or strictly
      better---path, since the iteration trajectory differs.

      \item \strong{Breaking:} Chebyshev-basis approximation now internally
      maps the requested range \eqn{[l, u]} to \eqn{[-1, 1]} before evaluating
      \eqn{T_k}, fixing the conditioning issue noted above. \code{a} (and
      \code{b}) are therefore coefficients of \eqn{T_k} of the internally
      mapped variable, not \eqn{T_k(x)}, whenever the fitted range is not
      \eqn{[-1, 1]}; results ON \eqn{[-1, 1]} are unchanged (the map is the
      identity there, verified bitwise-identical). \code{aMono}/\code{bMono}
      remain coefficients in the natural, raw \eqn{x} in every case, and are
      the stable way to consume a Chebyshev-basis result without needing to
      know about the internal map; see \code{\link{minimaxApprox}} for a
      documented accuracy caveat on this conversion at high degree on wide,
      far-off-\eqn{[-1,1]} ranges.

      \item Added upfront validation of \code{fn}, \code{lower}, \code{upper},
      \code{degree}, and \code{opts} in \code{minimaxApprox()}, replacing
      cryptic low-level failures (e.g. \dQuote{missing value where TRUE/FALSE
      needed}, \dQuote{invalid argument type}) with clear, argument-specific
      error messages. \code{fn} must be a function whose first argument is
      \code{x} (primitives such as \code{exp} and \code{sin} are still
      accepted); \code{lower}/\code{upper} must be finite non-missing numeric
      scalars; \code{degree} must be finite and non-missing;
      \code{opts$maxiter}/\code{miniter}/\code{conviter} must be single
      positive integers and \code{opts$tol}/\code{convrat}/\code{tailtol}/
      \code{ztol} must be finite numeric scalars when supplied
      (\code{tailtol}/\code{ztol} may still be \code{NULL} by design).

      \item Corrected wording in the degree-\code{n+1} restart messages, which
      described the returned polynomial's length as its degree
      (a degree-\eqn{n} polynomial has \eqn{n+1} coefficients).

      \item \code{isOscil}, used internally to test for the expected
      alternating-sign error pattern, no longer propagates \code{NA} when passed
      a \code{NaN}/\code{NA} error value; this previously could error deep
      inside the convergence check instead of failing cleanly.

      \item The internal zero-basis perturbation (relative-error mode, used when
      the approximated function is exactly 0 at a candidate extremum) now
      escalates to a magnitude-scaled step when the fixed \code{1e-12} nudge
      would be a no-op (points with \eqn{|x| \gtrsim 4.5 \times 10^3}), while
      remaining identical to the previous behavior at ordinary magnitudes.

      \item Functions that the requested degree resolves to (or near) machine
      floating-point precision---including functions that are exactly
      polynomials, such as \code{minimaxApprox(function(x) x^2, 0, 1, 2)}, and
      well-behaved cases such as \code{minimaxApprox(exp, -1, 1, 14)}---no
      longer hard-error or emit misleading non-convergence warnings. This case
      can surface two ways, both now handled: (1) both the degree-\eqn{n} and
      degree-\eqn{n+1} Remez solves fail singular, previously
      \dQuote{neither converged}; or (2) on some \acronym{BLAS}/\acronym{LAPACK}
      platforms the same case does not report singular but instead runs to
      \code{maxiter} (or stalls \dQuote{unchanging}) while genuinely wandering
      at the machine-precision floor. In either case, when there is no minimax
      problem left to solve, the degree-\eqn{n} interpolant through the
      Chebyshev reference is returned instead, together with a warning that it
      is not technically a Remez result but that its error is below the level
      at which a distinct minimax solution can be discerned. Genuinely
      non-representable cases (where no polynomial comes near the function)
      still raise the original error or convergence warning, unaffected.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added a reference-local certificate check to \strong{all}
      approximation paths (previously rational barycentric only): at a
      converged exit, the returned approximation's maximum error on a dense
      grid is compared against its leveled error, and a warning is raised when
      the former materially exceeds the latter. In that case the reported
      expected error is a valid \emph{lower bound} on the true minimax error
      (de la Vall\enc{é}{e}e Poussin) rather than a certificate of it, which
      can occur for degrees that are non-normal for the function (e.g. even or
      odd functions at generic degrees) or when the exchange stalls at a
      reference-local fixed point. The check is skipped where its ratio would
      be noise: results at the machine-precision floor and relative-error
      results with an exact zero of \code{fn} on the grid.

      \item Added a \code{"Barycentric"} option to the \code{basis} argument,
      implementing the barycentric-Remez algorithm for best polynomial
      approximation of Pach\enc{ó}{o}n and Trefethen (2009). It represents the
      trial polynomial by its values at the reference points rather than by
      coefficients in a fixed basis, avoiding the ill-conditioned linear solve
      that limits the achievable degree in the classical bases. In consequence
      it converges at degrees and on ranges where the Monomial and Chebyshev
      bases fail with a singular-matrix error---for example,
      \code{minimaxApprox(exp, -1, 1, 14)} (and degree 50), and \code{exp} on
      \eqn{[5, 6]} at degree 10---and, using capacity scaling, remains stable on
      very wide or far-off-center intervals. Exactly-representable and
      machine-precision-resolved functions are handled natively (leveled error
      \eqn{0}, no rescue needed). Relative error is supported; a function that is
      \strong{exactly} zero at a reference point (typically an interval
      endpoint), where the relative minimax does not exist, now raises a clear
      error. The returned object carries the barycentric representation and a
      coefficient-conversion residual diagnostic (\code{convResid}). See
      \code{\link{minimaxApprox}}.

      \item Added \strong{rational} approximation support to the
      \code{"Barycentric"} basis, implementing (a subset of) the adaptive
      barycentric minimax algorithm of Filip, Nakatsukasa, Trefethen, and
      Beckermann (2018), with the initial reference taken from the error
      extrema of an AAA approximant (Nakatsukasa, S\enc{è}{e}te, and
      Trefethen, 2018). The leveled error is computed as an eigenvalue of a
      small symmetric matrix, so, as for the polynomial case, there is no
      ill-conditioned linear solve; the classical smooth-function results are
      reproduced to convergence tolerance while cases such as
      \eqn{|x|} converge to the certified minimax where reference-based
      classical iterations can struggle. Absolute error only in this release
      (\code{relErr = TRUE} is an error); degenerate or defective requests
      (e.g. even degrees for an even function) are detected and reported with
      a clear error rather than repaired by degree reduction; a trial
      denominator with a zero in the interval is an error mirroring the
      classical rational path. The returned object carries \code{a}/\code{b}
      normalized as in the classical path, \code{aMono}/\code{bMono}, the
      barycentric representation, and a named two-component \code{convResid}.
      See \code{\link{minimaxApprox}}.

      \item Added the package's first vignette, \dQuote{Using minimaxApprox},
      covering the minimax criterion, the three bases and how to choose among
      them, rational approximation and its barycentric-basis limitations,
      relative error, the \code{Warning}/\code{convResid} diagnostics, and the
      handful of \code{opts} most users need.
    }
  }
}

\section{Version 0.5.0 (2026-01-08)}{
  \subsection{Fixed}{
    \itemize{
      \item Fixed error in one package metadata unit test.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item \strong{Breaking:} Converted \acronym{C} code for Compensated Horner
      method to use \code{long double} internally to eke out a bit more
      precision. One unit test had to be updated as what was formerly an error
      now converged.
      \item Updated Github-related files and workflows.
      \item Updated \file{WORDLIST}.
      \item Reduced overwriting of base \R function names.
      \item Simplified code and nesting of various functions.
      \item Package metadata tests now use \code{vapply} for type safety.
    }
  }
}

\section{Version 0.4.3 (2024-06-20)}{
  \subsection{Changed}{
    \itemize{
      \item Efficiency tweaks to C code.
      \item Use native \acronym{CRAN} \acronym{DOI} and update \file{README}.
      \item Expanded package metadata testing.
    }
  }
}

\section{Version 0.4.2 (2024-05-20)}{
  \subsection{Changed}{
    \itemize{
      \item If an explicitly passed \code{conviter} is larger than
      \code{miniter}, then \code{miniter} is set to \code{conviter}.
      \item Documentation made more consistent.
      \item C code tweaked for clarity and efficiency.
      \item Updated unit tests.
      \item Set up Github workflow for testing with \pkg{rhub}.
    }
  }
}

\section{Version 0.4.1 (2024-02-28)}{
  \subsection{Changed}{
    \itemize{
      \item One unit test adjusted to pass on Professor Ripley's OpenBLAS test
      machine which has different tolerances than the standard twelve on
      \acronym{CRAN}, six on Github, three on \R Win Builder, and the one on \R
      Mac Builder.
    }
  }
}

\section{Version 0.4.0 (2024-02-26)}{
  \subsection{Added}{
    \itemize{
      \item \strong{Breaking:} Analysis now uses the Chebyshev polynomials of
      the first kind, \eqn{T_k}, as the default basis instead of the monomials,
      \eqn{x^k}. Monomials can still be used by passing the appropriate
      \code{basis} parameter.
      \item The internal functions underlying the Chebyshev polynomial
      calculations were ported to C for speed.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item \strong{Breaking:} Changed names in return object. What was
      \code{Basis} is now \code{Extrema}. The \code{basis} parameter is now used
      to capture whether monomials or Chebyshev polynomials should be/were used.
      \item \strong{Breaking:} User-facing functions now have \code{basis}
      parameter to determine which basis to use; see appropriate help files.
      Parameter positions have thus changed.
      \item The \code{print}, \code{plot}, and \code{coef} methods were updated
      to recognize or allow basis selection.
      \item The \code{print} method now shows 14 digits for the convergence
      ratio by default.
      \item Internal code was streamlined and tolerances were tightened.
      \item Updated documentation for new functionality.
      \item The \acronym{C} files were split and header files created to keep
      Chebyshev and Horner code separate from each other.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Fixed potential division-by-zero problem in the \dQuote{unchanging}
      test.
      \item Fixed outdated parameter name in documentation.
      \item Fixed incorrect \file{NEWS} entries.
    }
  }
}

\section{Version 0.3.0 (2024-01-25)}{
  \subsection{Changed}{
    \itemize{
      \item \strong{Breaking:} More descriptive names in return object. See the
      \strong{Value} section of \code{help("minimaxApprox")} for details.
      \item \strong{Breaking:} When using relative error for convergence and one
      of the basis points has a functional value of 0, instead of stopping with
      an error, the algorithm will perturb that basis point by \code{1e-12} in
      the appropriate direction and proceed with a warning.
      \item \strong{Breaking:} Forced convergence to look \emph{only} at the
      ratio if the difference between tested values is \eqn{\le}{<=}
      \code{.Machine$double.eps}.
      \item \strong{Breaking:} Increased default for \code{conviter} from 10 to
      30 as tests on some values now need more iterations to achieve ratio
      convergence given a distance \eqn{\le}{<=} \code{.Machine$double.eps}.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added unit tests for additional functionality.
    }
  }
  \subsection{Fixed}{
    \itemize{
      \item Fixed a misnamed option value in some unit tests.
    }
  }
}

\section{Version 0.2.2 (2023-10-12)}{
  \subsection{Changed}{
    \itemize{
      \item Refactored \file{README} to split out security and contributions
      wordings to their own respective files and link the citation to the file
      on \acronym{CRAN}.
      \item Changed some tests to be more acceptable for \acronym{CRAN}.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added \code{LTO_OPT} to \file{Makevars} and \code{UseLTO} to
      \file{DESCRIPTION}.
    }
  }
}

\section{Version 0.2.1 (2023-09-04)}{
  \subsection{Changed}{
    \itemize{
      \item Test for calculating polynomial of degree \code{n + 1} when one of
      degree \code{n} fails, while passing both R \acronym{Mac} Builder and
      Github's \acronym{Mac} platform, failed on \acronym{CRAN}'s testbed.
      So test is no longer run for \acronym{Mac} on \acronym{CRAN}.
    }
  }
}

\section{Version 0.2.0 (2023-09-04)}{
  \subsection{Added}{
    \itemize{
      \item When the polynomial algorithm fails with a singular matrix error
      for degree \code{n}, it will now look for a polynomial of degree
      \code{n + 1}. If found, \strong{and} the contribution of that coefficient
      is less than the new option \code{tailtol} (default \code{1e-10}), it will
      return the found polynomial without the uppermost coefficient---which is
      effectively 0. It will message the user to that effect. Passing
      \code{NULL} skips the test. (see
      \href{https://github.com/aadler/minimaxApprox/issues/2}{Issue 2} with
      thanks to \href{https://github.com/hwborchers}{Hans W. Borchers}).
      \item New option \code{ztol} which acts similarly to
      \code{\link{zapsmall}} in that it will zero out coefficient values---both
      polynomial and rational---whose contributions are \eqn{\le} \code{ztol}
      at each step of the iteration. The default, \code{NULL}, skips the test.
    }
  }
  \subsection{Changed}{
    \itemize{
      \item Ported Compensated Horner Scheme algorithms to C for speed. The
      algorithms are now roughly 10x faster.
      \item Changed linear solution algorithm to use \code{qr.solve} when the
      default \code{solve} fails with a singular matrix issue. This \emph{may}
      allow for higher degree polynomials in certain cases.
      \item \strong{Breaking:} Changed name of \sQuote{\code{round}} parameter
      in print method to \sQuote{\code{digits}} to prevent any confusion with
      base function \code{round}.
      \item Updated unit tests.
      \item Cleaned up and updated some documentation and various messages.
      \item Tightened tolerances and cleaned up some internal functions.
      \item Added \file{LICENSE} to \file{.Rbuildignore} and removed from
      \file{DESCRIPTION} as suggested on
      \href{https://stat.ethz.ch/pipermail/r-package-devel/2023q3/009417.html}{R-pkg-devel}.
    }
  }
}

\section{Version 0.1.0 (2023-08-07)}{
  \subsection{Changed}{
    \itemize{
      \item Internal polynomial evaluations now uses Compensated Horner Scheme
      of Langlois et al. (2006) for increased precision and stability.
      \item Clarified some error messages.
    }
  }
  \subsection{Removed}{
    \itemize{
      \item Remove one test which fails on M1 Mac (thanks to Professor Brian
      Ripley) and PowerPC
      (\href{https://github.com/aadler/minimaxApprox/issues/1}{Issue 1};
      thanks to \href{https://github.com/barracuda156}{Sergey Fedorov}) due to
      it returning a singular matrix where on Linux and Windows it does not.
    }
  }
  \subsection{Added}{
    \itemize{
      \item Added file LICENSE with text of MPL-2.0.
      \item Added copyright and license notices to source code where it was
      missing.
      \item Added unit tests.
    }
  }
}

\section{Version 0.0.2 (2023-07-20)}{
  \subsection{Changed}{
    \itemize{
      \item Correct documentation and make more consistent.
      \item Provide better examples for \code{minimaxApprox}.
      \item Updated README with live CRAN links.
    }
  }
}

\section{Version 0.0.1 (2023-07-18)}{
  \subsection{Added}{
    \itemize{
      \item Initial release.
    }
  }
}
