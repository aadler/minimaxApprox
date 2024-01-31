# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

# Test print, plot, and coef methods
PP <- minimaxApprox(function(x) exp(x), 0, 1, 5L)
expect_identical(unlist(coef(PP), use.names = FALSE), PP$a)
expect_stdout(print(PP))
expect_stdout(plot(PP))

PP <- minimaxApprox(function(x) exp(x), 0, 1, 5L, TRUE)
expect_identical(unlist(coef(PP), use.names = FALSE), PP$a)
expect_stdout(print(PP))
expect_stdout(plot(PP))

PP <- minimaxApprox(function(x) exp(x), 0, 1, 5L, TRUE, "c")
expect_identical(unlist(coef(PP), use.names = FALSE), PP$a)
expect_stdout(print(PP))
expect_stdout(plot(PP))

RR <- minimaxApprox(function(x) exp(x), 0, 1, c(2L, 2L))
expect_identical(unlist(coef(RR)$a, use.names = FALSE), RR$a)
expect_identical(unlist(coef(RR)$b, use.names = FALSE), RR$b)
expect_stdout(print(RR))
expect_stdout(print(RR, digits = 2L))
expect_stdout(plot(RR))

RR <- suppressWarnings(minimaxApprox(exp, 0, 1, c(2L, 2L), TRUE))
expect_identical(unlist(coef(RR)$a, use.names = FALSE), RR$a)
expect_identical(unlist(coef(RR)$b, use.names = FALSE), RR$b)
expect_stdout(print(RR))
expect_stdout(plot(RR))
expect_stdout(plot(RR, ylim = c(-5e-6, 5e-6)))

RR <- suppressWarnings(minimaxApprox(exp, 0, 1, c(2L, 2L), TRUE, "c"))
expect_identical(unlist(coef(RR)$a, use.names = FALSE), RR$a)
expect_identical(unlist(coef(RR)$b, use.names = FALSE), RR$b)
expect_stdout(print(RR))
expect_stdout(plot(RR))

# Test ignore y message
mess <- paste("The y values are taken from the minimaxApprox object.",
              "Passed y values are ignored.")
expect_message(plot(RR, y = 1:4), mess)
