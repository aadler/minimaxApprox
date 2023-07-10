# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Print method
fN <- function(x, d = 6) formatC(x, digits = d, format = "e")

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
isConverged <- function(errs, E, cnvgRatio, tol) {
  abserrs <- abs(errs)
  mxae <- max(abserrs)
  # Check errors are of same magnitude by ratio or tolerance
  errMagnitude <- all(mxae / abserrs < cnvgRatio) || all(diff(abserrs) <= tol)
  # Check observed errors are close enough to expected by ratio or tolerance
  errDistance <- mxae / E < cnvgRatio || abs(mxae - E) <= tol

  # Converged if magnitude and distance are close and error oscillate
  isOscil(errs) && errMagnitude && errDistance
}

# Check denominator polynomial for zero in the requested range
checkDenom <- function(b, l, u) {
  dngrRt <- tryCatch(uniroot(polyCalc, c(l, u), a = b),
                     error = function(cond) simpleError(trimws(cond$message)))
  if (inherits(dngrRt, "simpleError")) {
    return(NULL)
  } else {
    return(dngrRt$root)
  }
}

# Print method (hide i and basis but leave in list and not attribute)
print.MiniMaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    ret <- list(b = x$b)
  } else {
    ret <- list(a = x$a, b = x$b)
  }

  ret <- c(ret, list(ExpectedError = x$EE,
                     ObservedError = x$OE,
                     Ratio = round(x$OE / x$EE, 6L),
                     Difference = abs(x$OE - x$EE),
                     Warnings = x$Warning))
  print(ret)
}

coef.MiniMaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    coef <- list(b = x$b)
  } else {
    coef <- list(a = x$a, b = x$b)
  }

  coef
}
# Plot method for errors and basis points
plot.MiniMaxApprox <- function(x, ...) {
  rng <- attr(x, "range")
  fn <- attr(x, "func")
  z <- seq(rng[1], rng[2], length.out = 1001L)

  if (attr(x, "type") == "Polynomial") {
    zz <- remPolyErr(z, x$b, fn)
    y <- remPolyErr(x$x, x$b, fn)
  } else {
    zz <- remRatErr(z, x$a, x$b, fn)
    y <- remRatErr(x$x, x$a, x$b, fn)
  }

  ybnd <- max(abs(x$EE), abs(x$OE))
  plot(z, zz, type = 'l',  xlab = "x", ylab = "Error", ylim = c(-ybnd, ybnd))
  abline(h = 0)
  points(x$x, y, col = "red", pch = 16)
  abline(h = c(-x$EE, x$EE), lty = 2, col = 'red')
  abline(h = c(-x$OE, x$OE), lty = 3, col = 'blue')
  legend(x = "bottomleft", inset = c(0.35, 1), col = c("red", "red", "blue"),
         lty = c(NA, 2, 3), legend = c("Basis", "Exp Err", "Obs Err"),
         pch = c(16, NA, NA), bg = "transparent", xpd = TRUE)
}
