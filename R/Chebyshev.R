# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

# Recursive function to calculate the kth Chebyshev polynomial T_k. Eventually,
# this should be split into pre-calculated known use cases and ported to C for
# speed, but for development purposes, keeping everything in R. The original is
# the recursive definition. Using the cos/cosh based version for speed.

# chebPoly <- function(x, k, maxExpo = 40L) {
#   k <- as.integer(k)
#   stopifnot(k >= 0,
#             k <= maxExpo,
#             is.numeric(x)
#   )
#   if (k == 0L) {
#     return(1)
#   } else if (k == 1L) {
#     return(x)
#   } else {
#     return(2 * x * chebPoly(x, k - 1L) - chebPoly(x, k - 2L))
#   }
# }

chebPoly <- function(x, k) {
  k <- as.integer(k)
  stopifnot(exprs = {
    k >= 0
    is.numeric(x)
  })
  suppressWarnings(ifelse(x <= -1,
                          (-1) ^ k * cosh(k * acosh(-x)),
                          ifelse(x <= 1L,
                                 cos(k * acos(x)),
                                 cosh(k * acosh(x)))))
}

# Create the equivalent of a Vandermonde matrix but using Chebyshev polynomials.
# Commented out version was original for development. Final version is an order
# of magnitude faster because it is more vectorized.
# chebMat <- function(x, n) {
#   k <- length(x)
#   np1 <- n + 1L
#   ret <- matrix(0, ncol = np1, nrow = k)
#   for (i in seq_len(np1)) {
#     ret[, i] <- chebPoly(x, i - 1L)
#   }
#   ret
# }

chebMat <- function(x, n) {
  outer(x, 0:n, chebPoly)
}

# Function to evaluate Chebyshev polynomials and their coefficient. Should port
# this to C eventually. Going to use matrix multiplication for now. Checked with
# VanderMat, this invocation is the "fastest" equivalent to polyCalc although an
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
