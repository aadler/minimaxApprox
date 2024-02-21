% Copyright Avraham Adler (c) 2023
% SPDX-License-Identifier: MPL-2.0+

\name{minimaxApprox-internal}
\encoding{UTF-8}

%% Internal-only functions need to be listed here wrapped in \alias{}

\title{Internal minimaxApprox Functions}
\description{Internal minimaxApprox functions}
\details{These are not to be called directly by the user.}

% In Chebyshev.R
\alias{chebMat}
\alias{chebCalc}
\alias{cheb2mon}

% In monomial.R
\alias{vanderMat}
\alias{polyCalc}

% In RemezPolynomial.R
\alias{polyMat}
\alias{polyCoeffs}
\alias{remPoly}

% In RemezRational.R
\alias{ratMat}
\alias{ratCoeffs}
\alias{remRat}

%In shared.R
\alias{fC}
\alias{chebNodes}
\alias{callFun}
\alias{isOscil}
\alias{evalFunc}
\alias{remErr}
\alias{findRoots}
\alias{switchX}
\alias{isConverged}
\alias{checkDenom}
\alias{checkIrrelevant}

\keyword{internal}
