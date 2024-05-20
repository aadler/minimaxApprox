# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

if (Sys.info()["nodename"] == "HOME") {
  pV <- packageVersion("minimaxApprox")
  cit <- toBibtex(citation("minimaxApprox"))
  nws <- news(package = "minimaxApprox")
  myOtherPkgs <- c("Delaporte", "lamW", "Pade", "revss", "MBBEFDLite")

  # Test CITATION has most recent package version
  expect_true(any(grepl(pV, cit), fixed = TRUE))

  # Test NEWS has most recent package version
  expect_true(any(grepl(pV, nws, fixed = TRUE)))

  # Test that NEWS has an entry with DESCRIPTION's Date
  expect_true(any(grepl(packageDate("minimaxApprox"), nws, fixed = TRUE)))

  # Test that CITATION doesn't contain the name of any other of my packages
  expect_false(any(sapply(myOtherPkgs, grepl, x = cit, fixed = TRUE)))

  # Test that NEWS doesn't contain the name of any other of my packages
  expect_false(any(sapply(myOtherPkgs, grepl, x = nws, fixed = TRUE)))
}
