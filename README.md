<!-- badges: start -->
<!-- Comment out for now
[![](https://www.r-pkg.org/badges/version-last-release/Delaporte)](https://cran.r-project.org/package=Delaporte)
[![](http://cranlogs.r-pkg.org/badges/last-month/Delaporte)](https://cran.r-project.org/package=Delaporte)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5880051.svg)](https://doi.org/10.5281/zenodo.5880051)
-->
 [![R-CMD-check](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aadler/minimaxApprox/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/aadler/minimaxApprox/branch/master/graph/badge.svg)](https://app.codecov.io/gh/aadler/minimaxApprox?branch=master)
[![OpenSSF Best Practices](https://bestpractices.coreinfrastructure.org/projects/7580/badge)](https://bestpractices.coreinfrastructure.org/projects/7580)
<!-- badges: end -->

# minimaxApprox
**minimaxApprox** is an `R` package which Implements the Remez algorithm for
polynomial function approximation and a form of the Fraser-Hart version
<doi:10.1007/BF02162506> for rational function approximation.

## Citation
If you use the package, please cite it as:

  Avraham Adler (2023). minimaxApprox: Implementation of Remez Algorithm for
                        Polynomial and Rational Function Approximation.
  R package version 0.0.1.

A BibTeX entry for LaTeX users is:

```
  @Manual{,
    author = {Avraham Adler},
    title = {minimaxApprox: Implementation of Remez Algorithm for Polynomial and Rational Function Approximation},
    year = {2023},
    note = {R package version 0.0.1},
  }
```

## Acknowledgements
The author is grateful to [Martin Maechler](https://stat.ethz.ch/~maechler/) for
suggestions which helped the author with MiniMax approximation.

## Contributions
Please ensure that all contributions comply with both
[R and CRAN standards for packages](https://cran.r-project.org/doc/manuals/r-release/R-exts.html).

### Versioning
This project attempts to follow [Semantic Versioning](https://semver.org/).

### Changelog
This project attempts to follow the changelog system at
[Keep a CHANGELOG](https://keepachangelog.com/).

### Dependencies
This project intends to have as few dependencies as possible. Please consider
that when writing code.

### Style
Please conform to this
[coding style guide](https://www.avrahamadler.com/coding-style-guide/) as best
possible.

### Documentation
Please provide valid .Rd files and **not** roxygen-style documentation.

### Tests
Please review the current test suite and supply similar `tinytest`-compatible
unit tests for all added functionality.

### Submission
If you would like to contribute to the project, it may be prudent to first
contact the maintainer via email. A request or suggestion may be raised as an
issue as well. To supply a pull request (PR), please:

 1. Fork the project and then clone into your own local repository
 2. Create a branch in your repository in which you will make your changes
 3. Ideally use -s to sign-off on commits under the
 [Developer Certificate of Origin](https://developercertificate.org/).
 4. If possible, sign commits using a GPG key.
 5. Push that branch and then create a pull request

At this point, the PR will be discussed and eventually accepted or rejected by
the lead maintainer.

## Roadmap
### Major

 * Consider moving some calculations into a compiled language (C, Fortran, C++)
 for efficiency
 * Research
 [barycentric representations](https://www.chebfun.org/publications/remez.pdf)
 to consider if possible to implement in R
 * Alternatively, research if possible to use [Rmpfr](https://CRAN.R-project.org/package=Rmpfr ) package for increased
 precision

### Minor

 * Write a vignette

## Security
### Expectations
This package is a calculation engine and requires no secrets or private
information. It is checked for memory leaks prior to releases to CRAN using
ASAN/UBSBAN. Dissemination is handled by CRAN. Bugs are reported via the tracker
and handled as soon as possible.

### Assurance
The threat model is that a malicious actor would "poison" the package code by
adding in elements having nothing to do with the package's purpose but which
would be used for malicious purposes. This is protected against by having the
email account of the maintainer—used for verification by CRAN—protected by a
physical 2FA device (Yubikey) which is carried by the lead maintainer.
