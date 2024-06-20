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
**minimaxApprox** is an `R` package which implements the algorithm of Remez
(1962) for polynomial minimax approximation and of Cody et al. (1968)
<doi:10.1007/BF02162506> for rational minimax approximation.

## Citation
If you use the package, please cite it as per
[CITATION](https://CRAN.R-project.org/package=minimaxApprox/citation.html).

## Acknowledgments
The author is grateful to [Martin Maechler](https://stat.ethz.ch/~maechler/) for
suggestions which helped the author's introduction to minimax approximation.

## Roadmap
### Major

 * Research
 [barycentric representations](https://www.chebfun.org/publications/remez.pdf)
 to consider if possible to implement in R.
 * Alternatively, research if possible to use [Rmpfr](https://CRAN.R-project.org/package=Rmpfr)
 package for increased precision.

### Minor

 * Write a vignette.

## Contributions
Please see
[CONTRIBUTING.md](https://github.com/aadler/minimaxApprox/blob/master/CONTRIBUTING.md).

## Security
Please see
[SECURITY.md](https://github.com/aadler/minimaxApprox/blob/master/SECURITY.md).
