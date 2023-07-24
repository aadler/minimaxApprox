# Copyright Avraham Adler (c) 2023
# SPDX-License-Identifier: MPL-2.0+

pV <- packageVersion("minimaxApprox")

# Test CITATION has most recent package version
expect_true(any(grepl(pV, toBibtex(citation("minimaxApprox")), fixed = TRUE)))

# Test NEWS has most recent package version
expect_true(any(grepl(pV, news(package = "minimaxApprox"), fixed = TRUE)))

# Test that NEWS has an entry with DESCRIPTION's Date
expect_true(any(grepl(packageDate("minimaxApprox"),
                      news(package = "minimaxApprox"), fixed = TRUE)))
