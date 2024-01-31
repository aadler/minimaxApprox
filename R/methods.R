# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Print method (hide i and basis/x but leave in list and not in attribute)
print.minimaxApprox <- function(x, digits = 6L, ...) {
  basis <- list(`Polynomial Basis` = attr(x, "basis"))
  if (attr(x, "type") == "Polynomial") {
    coefficients <- list(a = x$a)
  } else {
    coefficients <- list(a = x$a, b = x$b)
  }

  if (attr(x, "basis") == "Chebyshev") {
    monomialEq <- list(aMono = x$aMono)
    if ((attr(x, "type") == "Rational")) {
      monomialEq <- c(monomialEq, list(bMono = x$bMono))
    }
  } else {
    monomialEq <- NULL
  }

  digits <- as.integer(digits)
  diagnostics <- list(x$ExpErr,
                      x$ObsErr,
                      Ratio = round(x$ObsErr / x$ExpErr, digits),
                      Difference = abs(x$ObsErr - x$ExpErr),
                      Warnings = x$Warning)

  names(diagnostics)[1:2] <- if (attr(x, "relErr")) {
    c("ExpectedRelError", "ObservedRelError")
  } else {
    c("ExpectedAbsError", "ObservedAbsError")
  }

  print(c(basis, coefficients, monomialEq, diagnostics))
}

coef.minimaxApprox <- function(object, ...) {
  if (attr(object, "type") == "Polynomial") {
    coef <- list(a = object$a)
  } else {
    coef <- list(a = object$a, b = object$b)
  }

  coef
}

# Plot method for errors and basis points
plot.minimaxApprox <- function(x, y = NULL, ...) {
  if (!is.null(y)) {
    message("The y values are taken from the minimaxApprox object. ",
            "Passed y values are ignored.")
  }
  args <- list(...)
  rng <- attr(x, "range")
  fn <- attr(x, "func")
  monoB <- attr(x, "basis") == "Monomial"
  relErr <- attr(x, "relErr")
  z <- seq(rng[1], rng[2], length.out = 1001L)
  zz <- remErr(z, x, fn, relErr, monoB)
  y <- remErr(x$Extrema, x, fn, relErr, monoB)

  # Default y-axis label
  ylab <- if (relErr) "Relative Error" else "Absolute Error"

  # Default y-axis limits
  ylim <- if ("ylim" %in% names(args)) { # nolint
    args$ylim
  } else {
    ybnd <- max(x$ExpErr, x$ObsErr)
    c(-ybnd, ybnd)
  }

  plot(z, zz, type = "l", xlab = "x", ylab = ylab, ...)
  abline(h = 0)
  points(x$Extrema, y, col = "red", pch = 16L)
  abline(h = c(-x$ObsErr, x$ObsErr), lty = 2L, col = "blue")
  abline(h = c(-x$ExpErr, x$ExpErr), lty = 3L, col = "red")
  legend(x = "bottomleft", inset = c(0.35, 1), col = c("red", "red", "blue"),
         lty = c(NA, 3L, 2L), legend = c("Extrema", "Exp Err", "Obs Err"),
         pch = c(16L, NA, NA), bg = "transparent", xpd = TRUE)
  title(sub = paste("Polynomial Basis:", attr(x, "basis")))
}
