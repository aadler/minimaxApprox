% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{NEWS}
\alias{NEWS}

\title{NEWS file for the minimaxApprox package}

\section{Version 0.0.3 (2023-07-21)}{
  \subsection{Changed}{
    \itemize{
      \item Trap one test which fails on M1 Mac due to it returning a singular
      matrix where Linux and Windows do not. Per \R
      \href{https://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-R-under-macOS}{Installation and Administration}, the native
      Apple build \dQuote{...may give different numerical results from the
      more common ‘x86_64’ platforms...as ARM hardware lacks extended-precision
      floating-point operations.}
      \item Fails on PowerPC as well per
      \href{https://github.com/aadler/minimaxApprox/issues/1}{Issue 1}; thanks
      to \href{https://github.com/barracuda156}{Sergey Fedorov}.
      \item Add a startup message to that effect for the M1 Mac and PowerPC.
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
