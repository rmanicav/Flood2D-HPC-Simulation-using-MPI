//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: fluxes.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 10:52:30
//

// Include Files
//#include "//solver.h"
//#include "rt_nonfinite.h"
#include "fluxes_emxutil.h"
#include<math.h>

class fluxes
{

	void ffluxes(const double UP[5292], double n, const double dwsex[1764], const
		double dwsey[1764], const double dux[1764], const double duy[1764],
		const double dvx[1764], const double dvy[1764], double hextra, const
		double zc[1764], emxArray_real_T* F, emxArray_real_T* G, double
		* amax);
};