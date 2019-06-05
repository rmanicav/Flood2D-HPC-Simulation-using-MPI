//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: corrector.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 05-Jun-2019 10:44:42
//

// Include Files


// Function Definitions

//
// Arguments    : double U[5292]
//                const double F[5292]
//                const double G[5292]
//                double n
//                double dt2
//                double dt
//                const double sox[1764]
//                const double sfx[1764]
//                const double soy[1764]
//                const double sfy[1764]
//                double grav
// Return Type  : void
//
void corrector(double U[5292], const double F[5292], const double G[5292],
	double n, double dt2, double dt, const double sox[1764], const
	double sfx[1764], const double soy[1764], const double sfy[1764],
	double grav)
{
	int j;
	int k;
	for (j = 0; j < (int)(n - 1.0); j++) {
		for (k = 0; k < (int)(n - 1.0); k++) {
			U[j + 42 * k] = (U[j + 42 * k] - dt2 * (F[((int)((1.0 + (double)j) + 1.0)
				+ 42 * k) - 1] - F[j + 42 * k])) - dt2 * (G[j + 42 * ((int)((1.0 +
				(double)k) + 1.0) - 1)] - G[j + 42 * k]);
			if (U[j + 42 * k] < 0.0) {
				U[j + 42 * k] = 0.0;
			}

			U[j + 42 * k] = ((U[j + 42 * k] - dt2 * (F[((int)((1.0 + (double)j) + 1.0)
				+ 42 * k) - 1] - F[j + 42 * k])) - dt2 * (G[j + 42 * ((int)((1.0 +
				(double)k) + 1.0) - 1)] - G[j + 42 * k])) - dt * grav * (sox[j + 42 * k]
					+ sfx[j + 42 * k]);
			U[j + 42 * k] = ((U[j + 42 * k] - dt2 * (F[((int)((1.0 + (double)j) + 1.0)
				+ 42 * k) - 1] - F[j + 42 * k])) - dt2 * (G[j + 42 * ((int)((1.0 +
				(double)k) + 1.0) - 1)] - G[j + 42 * k])) - dt * grav * (soy[j + 42 * k]
					+ sfy[j + 42 * k]);
		}
	}
}

//
// File trailer for corrector.cpp
//
// [EOF]
//
