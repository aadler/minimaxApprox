# Test Citation
expect_true(any(grepl(packageVersion("minimaxApprox"),
                      toBibtex(citation("minimaxApprox")), fixed = TRUE)))
