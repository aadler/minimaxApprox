# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, relErr = FALSE, xi = NULL,
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
    # If actually passed then overwrite maxiter if conviter > maxiter
    opts$maxiter <- max(opts$maxiter, opts$conviter)
  }

  if (!("showProgress" %in% names(opts))) {
    opts$showProgress <- FALSE
  }

  if (!("convrat" %in% names(opts))) {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    opts$convrat <- 1.000000001
  }

  if (!("tol" %in% names(opts))) {
    opts$tol <- 1e-14
  }

  # Used for cases where we check polynomial degree n + 1.
  # See issue 2 https://github.com/aadler/minimaxApprox/issues/2
  if (!("tailtol" %in% names(opts))) {
    opts$tailtol <- min(1e-10, (upper - lower) / 1e6)
  }

  if (!("ztol" %in% names(opts))) {
    opts$ztol <- .Machine$double.eps
  }

  if (!is.logical(relErr)) {
    stop("Relative Error must be a logical value. ",
         "Default FALSE returns absolute error.")
  }

  if (length(degree) == 2L) {        # Rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE               # Polynomial approximation requested
    if (!is.null(xi)) {
      message("Polynomial approximation uses Chebyshev nodes for initial ",
              "guess. Any passed xi is ignored.")
    }
  } else {
    # All else is an error
    stop("Polynomial approximation takes one value for degree and rational ",
         "approximation takes a vector of two values for numerator and ",
         "denominator degrees. Any other inputs are invalid.")
  }

  # Call Calculation Functions
  mmA <- if (ratApprox) {
    remRat(fn, lower, upper, numerd, denomd, relErr, xi, opts)
  } else {
    tryCatch(remPoly(fn, lower, upper, as.integer(degree), relErr, opts),
             error = function(cond) simpleError(trimws(cond$message)))
  }

  # In response to issue 2, https://github.com/aadler/minimaxApprox/issues/2,
  # allow the polynomial algorithm to try degree n + 1 if it fails due to
  # singular error in degree n. IF the resulting highest coefficient contributes
  # less than tailtol (which defaults to 1e-10) to the result then consider it 0
  # and return the resulting degree n and message appropriately.

  if (inherits(mmA, "simpleError") &&
        grepl("singular", mmA$message, fixed = TRUE)) {
    if (!is.null(opts$tailtol)) {
      mmA <- tryCatch(remPoly(fn, lower, upper, as.integer(degree + 1L), relErr,
                              opts),
                      error = function(cond) simpleError(trimws(cond$message)))
      if (!inherits(mmA, "simpleError")) {
        xmax <- max(abs(lower), abs(upper))
        n <- length(mmA$a)
        if ((mmA$a[n] * xmax ^ (n - 1L)) <= opts$tailtol) {
          mess <- paste("The algorithm failed while looking for a polynomial",
                        "of degree", degree, "but successfully completed when",
                        "looking for a polynomial of degree", degree + 1L,
                        "with the largest coefficient's contribution to the",
                        "approximation <=", paste0(opts$tailtol, ":"), "the",
                        "tailtol option. The result is a polynomial of length",
                        degree, "as the uppermost coefficient is effectively",
                        "0.")
          mmA$a <- mmA$a[-length(mmA$a)]
          message(mess)
        } else {
          stop("The algorithm did not converge when looking for a polynomial ",
               "of length ", degree, " and when looking for a polynomial of ",
               "degree ", degree + 1L, " the uppermost coefficient is not ",
               "effectively zero.")
        }
      } else {
        stop("The algorithm neither converged when looking for a polynomial of",
             " length ", degree, " nor when looking for a polynomial of ",
             "degree ", degree + 1L, ".")
      }
    } else {
      stop("The algorithm did not converge when looking for a polynomial of ",
           "degree ", degree, " and NULL was passed to the tailtol option.")
    }
  }

  # Handle all warnings centrally
  gotWarning <- FALSE

  if (mmA$i >= opts$maxiter && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not achieved in ",
            mmA$i, " iterations.\n", "The ratio is ", fC(mmA$mxae / mmA$expe),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$unchanged && !mmA$converged) {
    warning("Convergence to requested ratio and tolerance not achieved in ",
            mmA$i, " iterations.\n", mmA$unchanging_i, " successive ",
            "calculated errors were too close to each other to warrant ",
            "further iterations.\nThe ratio is ", fC(mmA$mxae / mmA$expe),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$mxae < 10 * .Machine$double.eps) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal given floating point limitations.")
    gotWarning <- TRUE
  }

  coefficients <- if (ratApprox) list(a = mmA$a, b = mmA$b) else list(a = mmA$a)
  diagnostics <- list(EE = mmA$expe, OE = mmA$mxae,  iterations = mmA$i,
                      x = mmA$x, Warning = gotWarning)
  ret <- c(coefficients, diagnostics)
  attr(ret, "type") <- if (ratApprox) "Rational" else "Polynomial"
  attr(ret, "func") <- fn
  attr(ret, "range") <- c(lower, upper)
  attr(ret, "relErr") <- relErr
  attr(ret, "tol") <- opts$tol
  attr(ret, "convrat") <- opts$convrat
  class(ret) <- c("minimaxApprox", class(ret))

  ret
}

# Evaluation convenience function. Identical to evalFunc but tests for
# inheritance from minimaxApprox.
minimaxEval <- function(x, mmA) {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  evalFunc(x, mmA)
}
