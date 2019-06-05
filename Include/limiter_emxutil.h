//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: limiter_emxutil.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 09:53:45
//
#ifndef LIMITER_EMXUTIL_H
#define LIMITER_EMXUTIL_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "limiter_types.h"


// Function Declarations
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);

#endif

//
// File trailer for limiter_emxutil.h
//
// [EOF]
//
