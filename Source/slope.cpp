#include<cmath>
#include<iostream>
#include<string>

class slope
{
public:
	slope()
	{

	}

	// Function Definitions

	//
	// Arguments    : double u0
	//                double u1
	// Return Type  : double
	//
	double rt_hypotd_snf(double u0, double u1)
	{
		double y;
		double a;
		double b;
		a = std::abs(u0);
		b = std::abs(u1);
		if (a < b) {
			a /= b;
			y = b * std::sqrt(a * a + 1.0);
		}
		else if (a > b) {
			b /= a;
			y = a * std::sqrt(b * b + 1.0);
		}
		else if (b == NAN) {
			y = b;
		}
		else {
			y = a * 1.4142135623730951;
		}

		return y;
	}

	//
	// Arguments    : double u0
	//                double u1
	// Return Type  : double
	//
	double rt_powd_snf(double u0, double u1)
	{
		double y;
		double d0;
		double d1;
		if (u0 == NAN || u1 == NAN) {
			y = NAN;
		}
		else {
			d0 = std::abs(u0);
			d1 = std::abs(u1);
			if (isinf(u1)) {
				if (d0 == 1.0) {
					y = 1.0;
				}
				else if (d0 > 1.0) {
					if (u1 > 0.0) {
						y = 0;
					}
					else {
						y = 0.0;
					}
				}
				else if (u1 > 0.0) {
					y = 0.0;
				}
				else {
					y = 0;
				}
			}
			else if (d1 == 0.0) {
				y = 1.0;
			}
			else if (d1 == 1.0) {
				if (u1 > 0.0) {
					y = u0;
				}
				else {
					y = 1.0 / u0;
				}
			}
			else if (u1 == 2.0) {
				y = u0 * u0;
			}
			else if ((u1 == 0.5) && (u0 >= 0.0)) {
				y = std::sqrt(u0);
			}
			else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
				y = NAN;
			}
			else {
				y = pow(u0, u1);
			}
		}

		return y;
	}

	//
	// Intermediate Bed Slope values
	// Arguments    : const double h[1764]
	//                const double u[1764]
	//                const double v[1764]
	//                double ManN
	//                double hextra
	//                const double dzcx[1764]
	//                const double dzcy[1764]
	//                double cellsize
	//                double n
	//                double sox[1764]
	//                double soy[1764]
	//                double sfx[100]
	//                double sfy[100]
	// Return Type  : void
	//
	void fslope(const double h[1764], const double u[1764], const double v[1764],
		double ManN, double hextra, const double dzcx[1764], const double
		dzcy[1764], double cellsize, double n, double sox[1764], double soy
		[1764], double sfx[100], double sfy[100])
	{
		int j;
		int k;
		for (j = 0; j < 1764; j++) {
			sox[j] = dzcx[j] / cellsize;
			soy[j] = dzcy[j] / cellsize;
		}

		memset(&sfx[0], 0, 100U * sizeof(double));
		memset(&sfy[0], 0, 100U * sizeof(double));

		//  Intermediate frction slope values
		for (j = 0; j < (int)n; j++) {
			for (k = 0; k < (int)n; k++) {
				if (h[j + 42 * k] < 1.0) {
					sfx[j + 10 * k] = 0.0;
					sfy[j + 10 * k] = 0.0;
				}
				else {
					sfx[j + 10 * k] = u[j + 42 * k] * (ManN * ManN) * rt_hypotd_snf(u[j + 42
						* k], v[j + 42 * k]) / rt_powd_snf(h[j + 42 * k] + hextra,
							1.3333333333333333);
					sfy[j + 10 * k] = v[j + 42 * k] * (ManN * ManN) * rt_hypotd_snf(u[j + 42
						* k], v[j + 42 * k]) / rt_powd_snf(h[j + 42 * k] + hextra,
							1.3333333333333333);
				}
			}
		}
	}
};

//
// File trailer for slope.cpp
//
// [EOF]
//