---
title: Package minimaxApprox
---

<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version-last-release/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![](http://cranlogs.r-pkg.org/badges/last-month/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![](https://cranlogs.r-pkg.org/badges/grand-total/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![R-CMD-check](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/aadler/minimaxApprox/graph/badge.svg?token=6T0B933SEL)](https://app.codecov.io/gh/aadler/minimaxApprox)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/7580/badge)](https://bestpractices.coreinfrastructure.org/projects/7580)
<!-- badges: end -->

## Description
**minimaxApprox** is an `R` package which implements minimax approximation of
functions via the Remez (1962) algorithm for polynomials and the
Cody-Fraser-Hart (1968) <doi:10.1007/BF02162506> algorithm for rational
functions, as well as their barycentric formulations: the Pachón-Trefethen
(2009) <doi:10.1007/s10543-009-0240-1> algorithm for polynomials and the
Filip-Nakatsukasa-Trefethen-Beckermann (2018) <doi:10.1137/17M1132409> algorithm
for rational functions, which provide improved numerical stability at higher
degrees and on wider intervals.

## Citation
If you use the package, please cite it as per
[CITATION](https://CRAN.R-project.org/package=minimaxApprox/citation.html).

## Acknowledgments
The author is grateful to [Martin Maechler](https://stat.ethz.ch/~maechler/) for
suggestions which helped the author's introduction to minimax approximation.

## Roadmap
### Major

  * Remove the `xi` argument (deprecated in 0.6.0): the redesigned exchange
 derives its reference from the error curve directly, making a
 user-supplied initial reference unnecessary.
 * Exploit even/odd symmetry via the t = x² substitution: half-degree
 solves with structurally exact zero coefficients for functions of known
 parity on symmetric intervals.
 * Report the dense-grid supremum on non-converged exits, so that
 warned results carry an honest error bound rather than the final
 reference's leveled error.

### Minor

  * Make the near-machine-precision warning threshold magnitude-relative
 rather than absolute.
 * Validate that the function evaluates finitely on the interval, with a
 clear message rather than a downstream error.
 * Assorted message formatting (e.g. the NaN convergence-ratio text when
 the observed error is exactly zero).

## Contributions
Please see
[CONTRIBUTING.md](https://github.com/aadler/minimaxApprox/blob/master/CONTRIBUTING.md).

## Security
Please see
[SECURITY.md](https://github.com/aadler/minimaxApprox/blob/master/SECURITY.md).
