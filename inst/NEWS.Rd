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

\section{Version 0.2.0 (2023-09-01)}{
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
