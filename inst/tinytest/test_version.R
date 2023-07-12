# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

pV <- packageVersion("minimaxApprox")

# Test CITATION has most recent package version
expect_true(any(grepl(pV, toBibtex(citation("minimaxApprox")), fixed = TRUE)))

# Test NEWS has most recent pacakge version
expect_true(any(grepl(pV, news(package = "minimaxApprox"), fixed = TRUE)))
