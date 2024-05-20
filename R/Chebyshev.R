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

# Function to evaluate Chebyshev polynomials and their coefficients. Originally
# was drop(chebMat(x, length(a) - 1L) %*% a). Ported to C and uses R's C
# interface to DGEMV. This is 15%–25% faster than the ported chebMat and %*%.

chebCalc <- function(x, a) {
  .Call(chebCalc_c, as.double(x), as.double(a))
}

# Below based on open-source Netlib function "dconcm" in mathc90.
cheb2mon <- function(a) {
  n <- length(a)
  # The Chebyshev polynomials of order 0 and 1 ARE the monomials x^0 and x^1!
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
