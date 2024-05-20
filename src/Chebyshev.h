// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#ifndef MMA_MMA_H
#define MMA_MMA_H

#include <R.h>
#include <Rinternals.h>

extern SEXP chebMat_c(SEXP x, SEXP k);
extern SEXP chebCalc_c(SEXP x, SEXP a);
void F77_NAME(chebM_f)(double *x, int m, int n, double *ret);
extern SEXP chebMat_fc(SEXP x, SEXP k);
void F77_NAME(chebC_f)(double *x, int m, double *a, int n, double *ret);
extern SEXP chebCalc_fc(SEXP x, SEXP a);

#endif
