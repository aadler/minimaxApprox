// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#ifndef MMA_MMA_H
#define MMA_MMA_H

#include <R.h>
#include <Rinternals.h>

extern SEXP chebMat_c(SEXP x, SEXP k);
extern SEXP chebCalc_c(SEXP x, SEXP a);

#endif
