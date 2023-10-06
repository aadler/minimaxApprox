<!-- badges: start -->
[![](https://www.r-pkg.org/badges/version-last-release/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![](http://cranlogs.r-pkg.org/badges/last-month/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![](https://cranlogs.r-pkg.org/badges/grand-total/minimaxApprox)](https://cran.r-project.org/package=minimaxApprox)
[![R-CMD-check](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/aadler/minimaxApprox/branch/master/graph/badge.svg)](https://app.codecov.io/gh/aadler/minimaxApprox?branch=master)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/7580/badge)](https://bestpractices.coreinfrastructure.org/projects/7580)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8158855.svg)](https://doi.org/10.5281/zenodo.8158855)
<!-- badges: end -->

# minimaxApprox
**minimaxApprox** is an `R` package which implements the algorithm of Remez
(1962) for polynomial minimax approximation and of Cody et al. (1968)
<doi:10.1007/BF02162506> for rational minimax approximation.

## Citation
If you use the package, please cite it as:

  Adler A (2023). minimaxApprox: Implementation of Remez Algorithm for
  Polynomial and Rational Function Approximation.
  [doi: 10.5281/zenodo.8158855](https://doi.org/10.5281/zenodo.8158855),
  R package version 0.2.1.9000, https://CRAN.R-project.org/package=minimaxApprox

A BibTeX entry for LaTeX users is:

```
  @Manual{,
    author = {Avraham Adler},
    title = {minimaxApprox: Implementation of Remez Algorithm for Polynomial and Rational Function Approximation},
    year = {2023},
    url = {https://CRAN.R-project.org/package=minimaxApprox},
    doi = {10.5281/zenodo.8158855},
    note = {R package version 0.2.1.9000},
  }
```

## Acknowledgements
The author is grateful to [Martin Maechler](https://stat.ethz.ch/~maechler/) for
suggestions which helped the author with minimax approximation.

## Roadmap
### Major

 * Research if the "n + 1" check implemented for polynomial approximation can be
 applied to numerator and denominator of rational approximation.
 * Research
 [barycentric representations](https://www.chebfun.org/publications/remez.pdf)
 to consider if possible to implement in R.
 * Alternatively, research if possible to use [Rmpfr](https://CRAN.R-project.org/package=Rmpfr)
 package for increased precision.
 * Alternatively, consider using Chebyshev polynomials instead of monomials as
 basis for linear equation solutions.

### Minor

 * Write a vignette.

## Contributions
Please see
[CONTRIBUTING.md](https://github.com/aadler/minimaxApprox/blob/master/CONTRIBUTING.md).

## Security
Please see
[SECURITY.md](https://github.com/aadler/minimaxApprox/blob/master/SECURITY.md).
