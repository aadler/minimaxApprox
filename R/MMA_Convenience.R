# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

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
  } else if (onlyMono) {
    evalFunc(x, mmA, "m")
  } else {
    RR <- list(a = mmA$aMono)
    if ("bMono" %in% names(mmA)) RR <- c(RR, list(b = mmA$bMono))
    evalFunc(x, RR, "m")
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
