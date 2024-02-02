# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

# Create the equivalent of a Vandermonde matrix but using Chebyshev polynomials.
# The prior "outer" call on based on the old ChebPoly was ported to C to a
# nested loop. For some reason, calling outer on the C version of chebPoly did
# not work properly. But once the entirety of the matrix build was moved to C,
# and chebPoly removed completely, it became moot.

chebMat <- function(x, k) {
  .Call(chebMat_c, as.double(x), as.double(k))
}

# Function to evaluate Chebyshev polynomials and their coefficient. Originally
# was drop(chebMat(x, length(a) - 1L) %*% a). Ported to C and combines all three
# functions: chebPoly, chebMat, and chebCalc into one routine using DGEMV. This
# is 15%–25% faster than the ported chebMat and %*%. (AA: 2024-01-31)
#
# However, it is not good coding practice to have two seperate functions doing
# the same thing as that can introduce bugs if one is updated and one isn't, so
# while it bother me to calculate the matrix in C, pass it back to R, and pass
# it BACK to c for chebCalc, since chebMat has to live on its own for the matrix
# creation in polyMat and ratMat, it will be used for chebCalc as well. It
# remains 10%–15% faster than the %*% invocation, but not 15%–25% of the
# self-contained function. (AA: 2021-02-01)

chebCalc <- function(x, a) {
  .Call(chebCalc_c,
        .Call(chebMat_c, as.double(x), as.double(length(a) - 1)), as.double(a))
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
  # Chebyshev polynomials of order 0 and 1 ARE the monomials x^0 and x^1!
  if (n > 2L) {
    nm2 <- n - 2L
    tp <- 1
    for (j in seq_len(nm2)) {
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
  }
  a
}
