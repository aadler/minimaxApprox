# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

tol <- sqrt(.Machine$double.eps)

# Test minimaxEval
x <- seq(0.1, 0.4, 0.025)

################################################################################
# Test minimaxEval
x <- seq(0.1, 0.4, 0.025)

## Check asking for Chebyshev when Chebyshev was run.
### Polynomial
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$ExpErr))

### Rational
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L))
expect_true(all(exp(x) - minimaxEval(x, mmA) <= mmA$ExpErr))

## Check asking for monomial when Chebyshev was run.
### Polynomial
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_true(all(exp(x) - minimaxEval(x, mmA, "m") <= mmA$ExpErr))

### Rational
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L))
expect_true(all(exp(x) - minimaxEval(x, mmA, "m") <= mmA$ExpErr))

## Check asking for monomial when only monomial was run.
### Polynomial
mmA <- minimaxApprox(exp, 0, 0.5, 5L, basis = "m")
expect_true(all(exp(x) - minimaxEval(x, mmA, "m") <= mmA$ExpErr))

### Rational
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L), basis = "m")
expect_true(all(exp(x) - minimaxEval(x, mmA, "m") <= mmA$ExpErr))

## Check asking for Chebyshev when only monomial was run.
msgMsg <- "Analysis was run using only the monomial basis."

### Polynomial
mmA <- minimaxApprox(exp, 0, 0.5, 5L, basis = "m")
expect_message(minimaxEval(x, mmA, basis = "Cheb"), msgMsg)

### Rational
mmA <- minimaxApprox(exp, 0, 0.5, c(2L, 3L), basis = "m")
expect_message(minimaxEval(x, mmA, basis = "Cheb"), msgMsg)

## Check error trap for mmA object
errMsg <- "This function only works with 'minimaxApprox' objects."
expect_error(minimaxEval(x, sin), errMsg)

## Check not selecting proper basis
errMsg <- "Select either the 'M'onomial or 'C'hebyshev basis."
expect_error(minimaxEval(x, mmA, basis = "A"), errMsg)
expect_error(minimaxEval(x, mmA, basis = 4), errMsg)

################################################################################
# Test minimaxErr
x <- seq(0.1, 0.4, 0.025)
## Absolute
mmA <- minimaxApprox(exp, 0, 0.5, 5L)
expect_equal(minimaxEval(x, mmA) - exp(x), minimaxErr(x, mmA), tolerance = tol)

## Relative
mmA <- minimaxApprox(exp, 0, 0.5, 5L, TRUE, basis = "c")
expect_equal((minimaxEval(x, mmA) - exp(x)) / exp(x), minimaxErr(x, mmA),
             tolerance = tol)

## Check error trap
errMsg <- "This function only works with 'minimaxApprox' objects."
expect_error(minimaxErr(x, sin), errMsg)
