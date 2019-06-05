//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: limiter_types.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 09:53:45
//
#ifndef LIMITER_TYPES_H
#define LIMITER_TYPES_H

// Include Files
#include "rtwtypes.h"

// Type Definitions
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T
#endif

//
// File trailer for limiter_types.h
//
// [EOF]
//
