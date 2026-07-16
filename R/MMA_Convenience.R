# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Evaluation convenience function. Identical to evalFunc but tests for
# inheritance from minimaxApprox.
minimaxEval <- function(x, mmA, basis = "Chebyshev") {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  requestedbasis <- tolower(substr(basis, 1L, 1L))
  objBasis <- attr(mmA, "basis")
  onlyMono <- objBasis == "Monomial"
  isBary <- objBasis == "Barycentric"                              # M5
  if (!(requestedbasis %in% c("c", "m", "b"))) {                   # M5: & "b"
    stop("Select either the 'B'arycentric, 'M'onomial, or 'C'hebyshev basis.")
  }
  # M6 (F5 Option A): evalFunc's l/u are now required whenever the Chebyshev
  # branch is reached; the fitted range lives on the object as attr(range).
  rng <- attr(mmA, "range")

  # M5: for a barycentric object, the object's OWN basis wins (consistent with
  # current behavior where the object's basis is the default evaluator). The
  # stored barycentric representation is the most accurate path, so a defaulted
  # or explicit "b" request evaluates through it. An explicit "c"/"m" request
  # is honored via the converted coefficients, with a message (mirroring the
  # Monomial-only message pattern) since that conversion can be less accurate
  # (see the object's convResid).
  if (isBary) {
    if (missing(basis) || requestedbasis == "b") {
      return(evalFunc(x, mmA, "b", rng[1L], rng[2L]))
    }
    message("Analysis was run using the barycentric basis. Evaluating via the ",
            "converted ",
            if (requestedbasis == "c") "Chebyshev" else "monomial",
            " coefficients, which may be less accurate than the barycentric ",
            "representation (see the object's convResid).")
    if (requestedbasis == "c") {
      return(evalFunc(x, mmA, "c", rng[1L], rng[2L]))
    }
    # M5 Phase 2: a rational barycentric object also carries bMono; without
    # it this branch would evaluate the numerator polynomial alone.
    RR <- list(a = mmA$aMono)
    if ("bMono" %in% names(mmA)) RR <- c(RR, list(b = mmA$bMono))
    return(evalFunc(x, RR, "m", rng[1L], rng[2L]))
  }

  # Non-barycentric object: a "b" request has nothing to evaluate.  # M5
  if (requestedbasis == "b") {
    stop("Analysis was not run using the barycentric basis. Select the ",
         "'M'onomial or 'C'hebyshev basis.")
  }

  if (requestedbasis == "c") {
    if (onlyMono) {
      message("Analysis was run using only the monomial basis. Calculating ",
              "errors using monomials.")
      evalFunc(x, mmA, "m", rng[1L], rng[2L])
    } else {
      evalFunc(x, mmA, "c", rng[1L], rng[2L])
    }
  } else if (onlyMono) {
    evalFunc(x, mmA, "m", rng[1L], rng[2L])
  } else {
    # aMono/bMono are already in raw x (natural representation, unaffected by
    # M6's internal mapping), so this monomial-evaluation branch needs no map.
    RR <- list(a = mmA$aMono)
    if ("bMono" %in% names(mmA)) RR <- c(RR, list(b = mmA$bMono))
    evalFunc(x, RR, "m", rng[1L], rng[2L])
  }
}

# Minimax approximation error convenience function. Based on remErr but takes a
# completed mmA object with relErr and basis as attributes.
minimaxErr <- function(x, mmA) {
  if (!inherits(mmA, "minimaxApprox")) {
    stop("This function only works with 'minimaxApprox' objects.")
  }
  y <- callFun(attr(mmA, "func"), x)
  rng <- attr(mmA, "range")
  ret <- evalFunc(x, mmA, tolower(substr(attr(mmA, "basis"), 1L, 1L)),
                  rng[1L], rng[2L]) - y
  if (attr(mmA, "relErr")) ret <- ret / y

  ret
}
