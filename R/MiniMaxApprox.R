# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Master user-exposed function
minimaxApprox <- function(fn, lower, upper, degree, relErr = FALSE,
                          basis = "Chebyshev", xi = NULL, opts = list()) {

  basis <- tolower(substr(basis, 1L, 1L))
  if (!(basis %in% c("c", "m"))) {
    stop("Must select either 'C'hebyshev or 'm'onomial basis for analysis.")
  }

  # Handle configuration options
  nopts <- names(opts)
  if (!("maxiter" %in% nopts)) {
    opts$maxiter <- 100L
  }

  if (!("miniter" %in% nopts)) {
    opts$miniter <- 10L
  }

  if ("conviter" %in% nopts) {
    # If actually passed then overwrite both maxiter and miniter if conviter is
    # greather than either one.
    opts$maxiter <- max(opts$maxiter, opts$conviter)
    opts$miniter <- max(opts$miniter, opts$conviter)
  } else {
    opts$conviter <- 30L
  }

  if (!("showProgress" %in% nopts)) {
    opts$showProgress <- FALSE
  }

  if (!("convrat" %in% nopts)) {
    # Using 1 + 1e-9 - See Cody (1968) page 250. Can reasonably expect between
    # 9 & 12 significant figures.
    opts$convrat <- 1.000000001
  }

  if (!("tol" %in% nopts)) {
    opts$tol <- 1e-14
  }

  # Used for cases where we check polynomial degree n + 1.
  # See issue 2 https://github.com/aadler/minimaxApprox/issues/2
  if (!("tailtol" %in% nopts)) {
    opts$tailtol <- min(1e-10, (upper - lower) / 1e6)
  }

  if (!("ztol" %in% nopts)) {
    opts$ztol <- NULL
  }

  if (!is.logical(relErr)) {
    stop("Relative Error must be a logical value. ",
         "Default FALSE returns absolute error.")
  }

  if (any(degree < 0) || any(floor(degree) < degree)) {
    stop("Degrees must be integers of least 0 (constant).")
  }

  if (length(degree) == 2L) {         # Rational approximation requested
    numerd <- as.integer(degree[1L])
    denomd <- as.integer(degree[2L])
    ratApprox <- TRUE
  } else if (length(degree) == 1L) {
    ratApprox <- FALSE                # Polynomial approximation requested
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
    remRat(fn, lower, upper, numerd, denomd, relErr, basis, xi, opts)
  } else {
    tryCatch(remPoly(fn, lower, upper, as.integer(degree), relErr, basis, opts),
             error = function(e) simpleError(trimws(e$message)))
  }

  # In response to issue 2, https://github.com/aadler/minimaxApprox/issues/2,
  # allow the polynomial algorithm to try degree n + 1 if it fails due to
  # singular error in degree n. IF the resulting highest coefficient contributes
  # less than opts$tailtol to the result then consider it 0 and return the
  # resulting degree n and message appropriately.

  if (!ratApprox && inherits(mmA, "simpleError")) {
    if (grepl("singular", mmA$message, fixed = TRUE)) { # nolint unnecessary_nested_if_linter
      if (!is.null(opts$tailtol)) { # May be different error
        mmA <- tryCatch(remPoly(fn, lower, upper, as.integer(degree + 1L),
                                relErr, basis, opts),
                        error = function(e) simpleError(trimws(e$message)))
        if (inherits(mmA, "simpleError")) {
          stop("The algorithm neither converged when looking for a polynomial",
               " of length ", degree, " nor when looking for a polynomial of",
               " degree ", degree + 1L, ".")
        } else {
          xmax <- max(abs(lower), abs(upper))
          n <- length(mmA$a)
          if ((mmA$a[n] * xmax ^ (n - 1L)) <= opts$tailtol) {
            mess <- paste("The algorithm failed while looking for a polynomial",
                          "of degree", degree, "but successfully completed",
                          "when looking for a polynomial of degree",
                          degree + 1L, "with the largest coefficient's",
                          "contribution to the approximation <= the tailtol",
                          "option. The result is a polynomial of length",
                          degree, "as the uppermost coefficient is effectively",
                          "0.")
            mmA$a <- mmA$a[-n]
            message(mess)
          } else {
            stop("The algorithm did not converge when looking for a polynomial",
                 " of length ", degree, " and when looking for a polynomial of",
                 " degree ", degree + 1L, " the uppermost coefficient is not",
                 " effectively zero.")
          }
        }
      } else {
        stop("The algorithm did not converge when looking for a polynomial of ",
             "degree ", degree, " and NULL was passed to the tailtol option.")
      }
    }
  }

  # Handle all warnings centrally.
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
            "calculated solutions were too close to each other to warrant ",
            "further iterations.\nThe ratio is ",
            fC(mmA$mxae / mmA$expe, d = 14L),
            " times expected and the difference is ",
            fC(abs(mmA$mxae - mmA$expe)), " from the expected.")
    gotWarning <- TRUE
  }

  if (mmA$mxae < 10 * .Machine$double.eps) {
    warning("All errors very near machine double precision. The solution may ",
            "not be optimal given floating point limitations.")
    gotWarning <- TRUE
  }

  if (mmA$zeroBasisError) {
    warning("During convergence, the algorithm chose basis point(s) where the ",
            "functional value is 0. The basis point was perturbed by 1e-12, ",
            "but consider approximating using absolute---not relative---error.")
    gotWarning <- TRUE
  }

  coefficients <- if (ratApprox) {
    list(a = mmA$a, b = mmA$b)
  } else {
    list(a = mmA$a)
  }

  if (basis == "m") {
    monomialEq <- NULL
    polynomalBasis <- "Monomial"
  } else {
    monomialEq <- list(aMono = cheb2mon(mmA$a))
    polynomalBasis <- "Chebyshev"
    if (ratApprox) {
      monomialEq <- c(monomialEq, list(bMono = cheb2mon(mmA$b)))
      monomialEq <- mapply(`/`, monomialEq, monomialEq$bMono[1L],
                           SIMPLIFY = FALSE)
    }
  }

  diagnostics <- list(ExpErr = mmA$expe, ObsErr = mmA$mxae, iterations = mmA$i,
                      Extrema = mmA$x, Warning = gotWarning)
  ret <- c(coefficients, monomialEq, diagnostics)
  attr(ret, "type") <- if (ratApprox) "Rational" else "Polynomial"
  attr(ret, "basis") <- polynomalBasis
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
minimaxEval <- function(x, mmA, basis = "Chebyshev") {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  requestedbasis <- tolower(substr(basis, 1L, 1L))
  onlyMono <- attr(mmA, "basis") == "Monomial"
  if (!(requestedbasis %in% c("c", "m"))) {
    stop("Select either the 'M'onomial or 'C'hebyshev basis.")
  }
  if (requestedbasis == "c") {
    if (onlyMono) {
      message("Analysis was run using only the monomial basis. Calculating ",
              "errors using monomials.")
      evalFunc(x, mmA, "m")
    } else {
      evalFunc(x, mmA, "c")
    }
  } else {
    if (onlyMono) {
      evalFunc(x, mmA, "m")
    } else {
      RR <- list(a = mmA$aMono)
      if ("bMono" %in% names(mmA)) RR <- c(RR, list(b = mmA$bMono))
      evalFunc(x, RR, "m")
    }
  }
}

# Minimax approximation error convenience function. Based on remErr but takes a
# completed mmA object with relErr and basis as attributes.
minimaxErr <- function(x, mmA) {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  y <- callFun(attr(mmA, "func"), x)
  ret <- evalFunc(x, mmA, tolower(substr(attr(mmA, "basis"), 1L, 1L))) - y
  if (attr(mmA, "relErr")) ret <- ret / y

  ret
}
