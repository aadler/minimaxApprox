# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Printing convenience function
fC <- function(x, d = 6L, f = "g", w = -1) {
  formatC(x, digits = d, format = f, width = w)
}

# Default Chebyeshev nodes
chebNodes <- function(n, a, b) {
  sort(0.5 * (a + b + (b - a) * cos((2 * seq_len(n) - 1) * pi / (2 * n))))
}

# Create Vandermonde matrix to a polynomial degree n, or n + 1 terms
vanderMat <- function(x, n) {
  n1 <- n + 1L
  matrix(rep(x, each = n1) ^ (seq_len(n1) - 1L), ncol = n1, byrow = TRUE)
}

# Call the function being approximated on the points x
callFun <- function(fn, x) {
  if (!is.function(fn)) stop("Unable to parse function.")
  do.call(match.fun(fn), args = list(x = x))
}

# Check that the values passed are oscillating in sign
isOscil <- function(x) all(abs(diff(sign(x))) == 2)

# Calculate the polynomial approximation. Use in numer & denom for rationals
polyCalc <- function(x, a)  drop(vanderMat(x, length(a) - 1) %*% a)

# Check Remez iterations for convergence
isConverged <- function(errs, expE, cnvgRatio, tol) {
  aerrs <- abs(errs)
  mxae <- max(aerrs)
  mnae <- min(aerrs)

  # Check observed errors are close enough to expected by ratio or tolerance
  errDistance <- mxae / expE <= cnvgRatio || abs(mxae - expE) <= tol

  # Check observed errors are close enough to each other by ratio or tolerance
  errMag <- mxae / mnae <= cnvgRatio || mxae - mnae <= tol


  # Converged if magnitude and distance are close and error oscillates in sign
  isOscil(errs) && errDistance && errMag
}

# Check denominator polynomial for zero in the requested range
checkDenom <- function(a, l, u) {
  dngrRt <- tryCatch(uniroot(polyCalc, c(l, u), a = a),
                     error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(dngrRt, "simpleError")) {
    return(NULL)
  } else {
    return(dngrRt$root)
  }
}
