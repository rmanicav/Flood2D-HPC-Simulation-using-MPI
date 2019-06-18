#include<cmath>
#include<iostream>
#include<string>
using namespace std;
/// <summary>
/// 
/// 
/// </summary>
class slope
{
public:
	//
	/// <summary>
	/// 
	/// </summary>
	/// <param name="h"></param>
	/// <param name="u"></param>
	/// <param name="v"></param>
	/// <param name="ManN"></param>
	/// <param name="hextra"></param>
	/// <param name="dzcx"></param>
	/// <param name="dzcy"></param>
	/// <param name="cellsize"></param>
	/// <param name="n"></param>
	/// <param name="sox"></param>
	/// <param name="soy"></param>
	/// <param name="sfx"></param>
	/// <param name="sfy"></param>
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