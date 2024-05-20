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

\section{Version 0.4.1.9000 (2024-03-28)}{
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
      Chebyshev and Horner code seperate from each other.
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
      30 as tests on some values now need more iterations to acheive ratio
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
