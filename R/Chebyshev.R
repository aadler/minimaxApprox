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
