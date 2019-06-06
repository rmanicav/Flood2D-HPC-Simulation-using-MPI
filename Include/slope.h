//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: slope.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 10:15:32
//

// Include Files
#include <cmath>
#include <math.h>
#include <string.h>
#include "fluxes_emxutil.h"

class slope
{
public:
	// Function Declarations
	double rt_hypotd_snf(double u0, double u1);
	double rt_powd_snf(double u0, double u1);
	void fslope(const double h[1764], const double u[1764], const double v[1764],
		double ManN, double hextra, const double dzcx[1764], const double
		dzcy[1764], double cellsize, double n, double sox[1764], double soy
		[1764], double sfx[100], double sfy[100]);
};
