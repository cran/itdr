#include "itdr.h"
#include <R.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>
#include <stdlib.h>
#include "R_ext/BLAS.h"
#include <math.h>
#include <stddef.h> //for NULL

static const R_CallMethodDef callMethods[]={
  {"dlogden1",(DL_FUNC) &dlogden1,7},
  {"dlogden3",(DL_FUNC) &dlogden3,8},
  {"ITM_mean_norm",(DL_FUNC) &ITM_mean_norm,7},
  {"ITM_mean",(DL_FUNC) &ITM_mean,8},
  {"ITM_pdf_norm",(DL_FUNC) &ITM_pdf_norm,8},
  {"ITM_pdf",(DL_FUNC) &ITM_pdf,9},
  {"Fdlogden1",(DL_FUNC) &Fdlogden1,7},
  {"Fdlogden3",(DL_FUNC) &Fdlogden3,8},
  {"FM_mean_norm",(DL_FUNC) &FM_mean_norm,7},
  {"FM_mean",(DL_FUNC) &FM_mean,8},
  {"FM_pdf_norm",(DL_FUNC) &FM_pdf_norm,8},
  {"FM_pdf",(DL_FUNC) &FM_pdf,9},
  {NULL,NULL,0}
};

static R_NativePrimitiveArgType dlogden1_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType dlogden3_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType ITM_mean_norm_t[] = {
  REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType ITM_mean_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType ITM_pdf_norm_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType ITM_pdf_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};

static R_NativePrimitiveArgType Fdlogden1_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType Fdlogden3_t[] = {
  REALSXP,INTSXP,INTSXP,REALSXP,INTSXP,REALSXP,REALSXP,REALSXP
};
static R_NativePrimitiveArgType FM_mean_norm_t[] = {
  REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType FM_mean_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType FM_pdf_norm_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static R_NativePrimitiveArgType FM_pdf_t[] = {
  REALSXP,REALSXP,REALSXP,REALSXP,REALSXP,INTSXP,INTSXP,INTSXP,REALSXP
};
static const R_CMethodDef cMethods[] = {
  {"dlogden1", (DL_FUNC) &dlogden1, 7,dlogden1_t},
  {"dlogden3", (DL_FUNC) &dlogden3, 8,dlogden3_t},
  {"ITM_mean_norm", (DL_FUNC) &ITM_mean_norm, 7,ITM_mean_norm_t},
  {"ITM_mean", (DL_FUNC) &ITM_mean, 8,ITM_mean_t},
  {"ITM_pdf_norm", (DL_FUNC) &ITM_pdf_norm, 8,ITM_pdf_norm_t},
  {"ITM_pdf", (DL_FUNC) &ITM_pdf, 9,ITM_pdf_t},
  {"Fdlogden1", (DL_FUNC) &Fdlogden1, 7,Fdlogden1_t},
  {"Fdlogden3", (DL_FUNC) &Fdlogden3, 8,Fdlogden3_t},
  {"FM_mean_norm", (DL_FUNC) &FM_mean_norm, 7,FM_mean_norm_t},
  {"FM_mean", (DL_FUNC) &FM_mean, 8,FM_mean_t},
  {"FM_pdf_norm", (DL_FUNC) &FM_pdf_norm, 8,FM_pdf_norm_t},
  {"FM_pdf", (DL_FUNC) &FM_pdf, 9,FM_pdf_t},
  {NULL, NULL, 0, NULL}
};


//static R_NativePrimitiveArgType vecMfmn_t[] = {
//  REALSXP, REALSXP,INTSXP,INTSXP,REALSXP,REALSXP
//};

//static const R_CMethodDef cMethods[] = {
//  {"vecMfmn", (DL_FUNC) &dr2_interface, 6,vecMfmn_t},
//  {NULL, NULL, 0, NULL}
//};

//void R_init_dr2(DllInfo *info)
//{
//  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
//}


void R_init_itdr(DllInfo *dll)
{
  R_registerRoutines(dll, cMethods, NULL
                       , NULL,NULL);
  R_useDynamicSymbols(dll, FALSE);
  //R_forceSymbols(dll, TRUE);
}
