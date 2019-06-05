//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fluxes_emxutil.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 10:52:30
//
#ifndef FLUXES_EMXUTIL_H
#define FLUXES_EMXUTIL_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "fluxes_types.h"


// Function Declarations
extern void emxEnsureCapacity_int32_T(emxArray_int32_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int b_numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int b_numDimensions);

#endif

//
// File trailer for fluxes_emxutil.h
//
// [EOF]
//
