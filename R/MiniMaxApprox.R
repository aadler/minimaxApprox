# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, errType = "abs", xi = NULL,
                          opts = list()) {

  # Handle configuration options
  if (!("maxiter" %in% names(opts))) {
    opts$maxiter <- 100L
  }

  if (!("miniter" %in% names(opts))) {
    opts$miniter <- 10L
  }

  if (!("conviter" %in% names(opts))) {
    opts$conviter <- 10L
  } else {
    # If actually passed then overwrite
    opts$maxiter <- max(opts$maxiter, opts$conviter)
  }

  if (!("showProgress" %in% names(opts))) {
    opts$showProgress <- FALSE
  }

  if (!("convRatio" %in% names(opts))) {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    opts$convRatio <- 1.000000001
  }

  if (!("tol" %in% names(opts))) {
    opts$tol <- 1e-14
  }

  absErr <- switch(tolower(substring(errType, 1, 3)),
                   abs = TRUE,
                   rel = FALSE,
                   stop("Error type must be either 'abs'olute or 'rel'ative."))

  if (length(degree) == 2L) {        # Rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE               # Polynomial approximation requested
    if (!is.null(xi)) {
      warning("Polynomial approximation uses Chebyeshev nodes for initial ",
              "guess. Any passed xi is ignored.")
    }
  } else {
    # All else is an error
    stop("Polynomial approximation takes one value for degree and rational ",
         "approximation takes a vector of two values for numerator and ",
         "denominator. Any other inputs are invalid.")
  }

  # CaLL Calculation Functions
  mmA <- if (ratApprox) {
    remRat(fn, lower, upper, numerd, denomd, absErr, xi, opts)
  } else {
    remPoly(fn, lower, upper, as.integer(degree), absErr, opts)
  }

  # Handle all warnings centrally
  gotWarning <- FALSE

  if (mmA$i >= opts$maxiter && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not acheived in ",
            mmA$i, " iterations.\n", "The ratio is ", fC(mmA$mxae / mmA$expe),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$unchanged && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not acheived in ",
            mmA$i, " iterations.\n", mmA$unchanging_i, " succesive calculated ",
            "errors were too close to each other to warrant further ",
            "iterations.\nThe ratio is ", fC(mmA$mxae / mmA$expe), " times ",
            "expected and the difference is ", fC(abs(mmA$mxae - mmA$expe)),
            " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$mxae < 10 * .Machine$double.eps) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal but should be best given the desired precision ",
            "and floating point limitations. Try a lower degree if needed.")
    gotWarning <- TRUE
  }

  coefficients <- if (ratApprox) list(a = mmA$a, b = mmA$b) else list(a = mmA$a)
  diagnostics <- list(EE = mmA$expe, OE = mmA$mxae,  iterations = mmA$i,
                      x = mmA$x, Warning = gotWarning)
  ret <- c(coefficients, diagnostics)
  attr(ret, "type") <- if (ratApprox) "Rational" else "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "absErr") <- absErr
  attr(ret, "tol") <- opts$tol
  attr(ret, "convRatio") <- opts$convRatio
  class(ret) <- c("minimaxApprox", class(ret))

  ret
}

# Evaluation convenience function
minimaxEval <- function(x, mmA) {
  if (attr(mmA, "type") == "Polynomial") {
    polyCalc(x, mmA$a)
  } else {
    remRatFunc(x, mmA$a, mmA$b)
  }
}

# Print method (hide i and basis/x but leave in list and not in attribute)
print.minimaxApprox <- function(x, ...) {
  if (attr(x, "type") == "Polynomial") {
    coefficients <- list(a = x$a)
  } else {
    coefficients <- list(a = x$a, b = x$b)
  }

  diagnostics <- list(x$EE,
                      x$OE,
                      Ratio = round(x$OE / x$EE, 6L),
                      Difference = abs(x$OE - x$EE),
                      Warnings = x$Warning)

  names(diagnostics)[1:2] <- if (attr(x, "absErr")) {
    c("ExpectedAbsError", "ObservedAbsError")
  } else {
    c("ExpectedRelError", "ObservedRelError")
  }

  print(c(coefficients, diagnostics))
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
plot.minimaxApprox <- function(x, y, ...) {
  args <- list(...)
  rng <- attr(x, "range")
  fn <- attr(x, "func")
  absErr <- attr(x, "absErr")
  z <- seq(rng[1], rng[2], length.out = 1001L)

  if (attr(x, "type") == "Polynomial") {
    zz <- remPolyErr(z, x$a, fn, absErr)
    y <- remPolyErr(x$x, x$a, fn, absErr)
  } else {
    zz <- remRatErr(z, x$a, x$b, fn, absErr)
    y <- remRatErr(x$x, x$a, x$b, fn, absErr)
  }

  ylab <- if (absErr) "Absolute Error" else "Relative Error"
  ybnd <- max(x$EE, x$OE)
  if ("ylim" %in% names(args)) {
    ylim <- args$ylim
  } else {
    ylim <- c(-ybnd, ybnd)
  }

  plot(z, zz, type = "l", xlab = "x", ylab = ylab, ...)
  abline(h = 0)
  points(x$x, y, col = "red", pch = 16L)
  abline(h = c(-x$OE, x$OE), lty = 2L, col = "blue")
  abline(h = c(-x$EE, x$EE), lty = 3L, col = "red")
  legend(x = "bottomleft", inset = c(0.35, 1), col = c("red", "red", "blue"),
         lty = c(NA, 3L, 2L), legend = c("Basis", "Exp Err", "Obs Err"),
         pch = c(16L, NA, NA), bg = "transparent", xpd = TRUE)
}
