# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

# Recursive function to calculate the kth Chebyshev polynomial T_k. Using the
# cos/cosh based version instead of the recursive definition for speed. Ported
# to C for speed since it and the outer call are the slowest elements.

chebPoly <- function(x, k) {
  stopifnot(exprs = {
    k >= 0
    is.numeric(x)
    # Need k to be a real for C purposes (pow) but still should LOOK like an
    # integer.
    all.equal(k - floor(k), 0)
  })
  .Call(chebPoly_c, as.double(x), floor(k))
}

# Create the equivalent of a Vandermonde matrix but using Chebyshev polynomials.

chebMat <- function(x, n) {
  outer(x, 0:n, chebPoly)
}

# Function to evaluate Chebyshev polynomials and their coefficient. May port
# this to C eventually. Going to use matrix multiplication for now. Checked with
# vanderMat, this invocation is the "fastest" equivalent to polyCalc although an
# order of magnitude slower, its an order of magnitude faster than rowSums on
# "sweep"ing the multiplication.

chebCalc <- function(x, a) {
  drop(chebMat(x, length(a) - 1L) %*% a)
}

evalFuncCheb <- function(x, R) {
  ret <- chebCalc(x, R$a)
  if ("b" %in% names(R)) {
    ret <- ret / chebCalc(x, R$b)
  }
  ret
}

# Based on open-source Netlib function "dconcm" in mathc90.

cheb2mon <- function(a) {
  n <- length(a)
  nm2 <- n - 2L
  tp <- 1
  for (j in 1:nm2) {
    # Cannot vectorize next step since already adjusted "a"s affect downstream
    # "a"s as part of the recursion.
    for (i in nm2:j) {
      a[i] <- a[i] - a[i + 2L]
    }
    a[j + 1L] <- a[j + 1L] / 2
    a[j] <- a[j] * tp
    tp <- tp * 2
  }
  a[c(n - 1L, n)] <- a[c(n - 1L, n)] * tp
  a
}
