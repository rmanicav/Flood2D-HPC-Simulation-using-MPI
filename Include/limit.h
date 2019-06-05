//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: limit.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 03-Jun-2019 14:16:12
//
#ifndef LIMIT_H
#define LIMIT_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
// Function Definitions

//
// Superbee limiterv  (df1x,df2x,df1m,df1y,df2y,df2m)
// Arguments    : double df1x
//                double df2x
//                double df1y
//                double df2y
//                double *f1
//                double *f2
// Return Type  : void
//
void limit(double df1x, double df2x, double df1y, double df2y, double* f1,
	double* f2)
{
	double s;
	double a;
	double b;
	double varargin_1[3];
	int idx;
	int k;
	bool exitg1;
	if (df1x * df2x < 0.0) {
		*f1 = 0.0;
	}
	else {
		s = df1x;
		if (df1x < 0.0) {
			s = -1.0;
		}
		else if (df1x > 0.0) {
			s = 1.0;
		}
		else {
			if (df1x == 0.0) {
				s = 0.0;
			}
		}

		a = std::abs(df1x);
		b = std::abs(df2x);
		varargin_1[0] = 2.0 * a;
		varargin_1[1] = 2.0 * b;
		varargin_1[2] = 0.5 * (a + b);
		if (varargin_1[0] != NAN) {
			idx = 1;
		}
		else {
			idx = 0;
			k = 2;
			exitg1 = false;
			while ((!exitg1) && (k < 4)) {
				if (varargin_1[k - 1] != NAN) {
					idx = k;
					exitg1 = true;
				}
				else {
					k++;
				}
			}
		}

		if (idx == 0) {
			a = varargin_1[0];
		}
		else {
			a = varargin_1[idx - 1];
			while (idx + 1 < 4) {
				if (a > varargin_1[idx]) {
					a = varargin_1[idx];
				}

				idx++;
			}
		}

		*f1 = s * a;
	}

	if (df1y * df2y < 0.0) {
		*f2 = 0.0;
	}
	else {
		s = df1y;
		if (df1y < 0.0) {
			s = -1.0;
		}
		else if (df1y > 0.0) {
			s = 1.0;
		}
		else {
			if (df1y == 0.0) {
				s = 0.0;
			}
		}

		a = std::abs(df1y);
		b = std::abs(df2y);
		varargin_1[0] = 2.0 * a;
		varargin_1[1] = 2.0 * b;
		varargin_1[2] = 0.5 * (a + b);
		if (varargin_1[0] !=NAN) {
			idx = 1;
		}
		else {
			idx = 0;
			k = 2;
			exitg1 = false;
			while ((!exitg1) && (k < 4)) {
				if (varargin_1[k - 1] != NAN) {
					idx = k;
					exitg1 = true;
				}
				else {
					k++;
				}
			}
		}

		if (idx == 0) {
			a = varargin_1[0];
		}
		else {
			a = varargin_1[idx - 1];
			while (idx + 1 < 4) {
				if (a > varargin_1[idx]) {
					a = varargin_1[idx];
				}

				idx++;
			}
		}

		*f2 = s * a;
	}

	//
	//  if(df1x*df2x < 0),
	//          f1=0;
	//    else
	//          s=sign(df1x);
	//          a=abs(df1x);
	//          b=abs(df2x);
	//          f1=s*min(max([a b]),2.0*min([a b]));
	//  end
	//
	//  if(df1y*df2y < 0),
	//          f2=0;
	//    else
	//          s=sign(df1y);
	//          a=abs(df1y);
	//          b=abs(df2y);
	//          f2=s*min(max([a b]),2.0*min([a b]));
	//  end
}

//
// File trailer for limit.cpp
//
// [EOF]
//

#endif

//
// File trailer for limit.h
//
// [EOF]
//
