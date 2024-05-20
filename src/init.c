// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "compHorner.h"
#include "Chebyshev.h"

void F77_NAME(chebM_f)(double *x, int m, int n, double *ret);

extern SEXP chebMat_fc(SEXP x, SEXP k) {
  const int m = LENGTH(x);
  const int n = asInteger(k) + 1;
  SEXP ret = PROTECT(allocMatrix(REALSXP, m, n));
  F77_CALL(chebM_f)(REAL(x), m, n, REAL(ret));
  UNPROTECT(1);
  return(ret);
}

void F77_NAME(chebC_f)(double *x, int m, double *a, int n, double *ret);

extern SEXP chebCalc_fc(SEXP x, SEXP a) {
  const int m = LENGTH(x);
  const int n = LENGTH(a);
  SEXP ret = PROTECT(allocVector(REALSXP, m));
  F77_CALL(chebC_f)(REAL(x), m, REAL(a), n, REAL(ret));
  UNPROTECT(1);
  return(ret);
}

static const R_CallMethodDef CallEntries[] = {
  {"compHorner_c",      (DL_FUNC) &compHorner_c,    2},
  {"chebMat_c",         (DL_FUNC) &chebMat_c,       2},
  {"chebCalc_c",        (DL_FUNC) &chebCalc_c,      2},
  {"chebMat_fc",        (DL_FUNC) &chebMat_fc,      2},
  {"chebCalc_fc",       (DL_FUNC) &chebCalc_fc,     2},
  {NULL,                NULL,                       0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
