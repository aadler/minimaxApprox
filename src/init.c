// Copyright Avraham Adler (c) 2024
// SPDX-License-Identifier: MPL-2.0+

#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "compHorner.h"
#include "Chebyshev.h"

static const R_CallMethodDef CallEntries[] = {
  {"compHorner_c",    (DL_FUNC) &compHorner_c,  2},
  {"chebMat_c",       (DL_FUNC) &chebMat_c,     2},
  {"chebCalc_c",      (DL_FUNC) &chebCalc_c,    2},
  {NULL,              NULL,                     0}
};

void R_init_minimaxApprox(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
