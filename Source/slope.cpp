#include<cmath>
#include<iostream>
#include<string>
using namespace std;
class slope
{
public:
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
	void fslope(double** h, double** u, double** v,
		double ManN, double hextra, double** dzcx, double**
		dzcy, int cellsize, int n, double** sox, double** soy
		, double **sfx, double** sfy)
	{

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				sox[i][j] = dzcx[i][j] / cellsize;
				soy[i][j] = dzcy[i][j] / cellsize;
			}
		}	

		//  Intermediate frction slope values
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				if (h[j][k] < 1.0) {
					sfx[j][k] = 0.0;
					sfy[j][k] = 0.0;
				}
				else {
					sfx[j][k] = u[j][k] * (ManN * ManN) * hypot(u[j][k], v[j][k]) / pow(h[j][k] + hextra,
							1.3333333333333333);
					sfy[j][k] = v[j][k] * (ManN * ManN) * hypot(u[j][k], v[j][k]) / pow(h[j][k] + hextra,					1.3333333333333333);
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