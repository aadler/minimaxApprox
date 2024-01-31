# Copyright Avraham Adler (c) 2024
# SPDX-License-Identifier: MPL-2.0+

# Create Vandermonde matrix to a polynomial degree n, or n + 1 terms.
vanderMat <- function(x, n) {
  np1 <- n + 1L
  matrix(rep(x, each = np1) ^ (seq_len(np1) - 1L), ncol = np1, byrow = TRUE)
}

# polyCalc uses a Compensated Horner Method based on Langlois et al.(2006)
# https://drops.dagstuhl.de/opus/volltexte/2006/442/
# As primary bottleneck, it was ported to C for speed.
polyCalc <- function(x, a) {
  .Call(compHorner_c, as.double(x), as.double(a))
}

# Function to calculate value of minimax approximation at x given a & b.
evalFuncMono <- function(x, R) {
  ret <- polyCalc(x, R$a)
  if ("b" %in% names(R)) {
    ret <- ret / polyCalc(x, R$b)
  }
  ret
}
