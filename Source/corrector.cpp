

class corrector
{
public:
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
	void fcorrector(double ***U, double ***F,  double ***G,
		double n, double dt2, double dt, double **sox, 
		double **sfx, double **soy, double **sfy,
		double grav)
	{
		int j;
		int k;
		for (j = 0; j < (int)(n - 1.0); j++) {
			for (k = 0; k < (int)(n - 1.0); k++) {
				U[j][k][0] = (U[j][k][0] - dt2 * (F[j][k - 1][0] - F[j][k][0])) - dt2 * (G[j][k - 1][0] - G[j][k][0]);
				if (U[j][k][0] < 0.0) {
					U[j][k][0] = 0.0;
				}

				U[j][k][1] = ((U[j][k][1] - dt2 * (F[j][k - 1][1] - F[j][k][1])) - dt2 * (G[j][k - 1][1] - G[j][k][1])) - dt * grav * (sox[j][k]
						+ sfx[j][k]);
				U[j][k][2] = ((U[j][k][2] - dt2 * (F[j][k - 1][2] - F[j][k][2]) - dt2 * G[j][k - 1][2] - G[j][k][2])) - dt * grav * (soy[j][k]
						+ sfy[j][k]);
			}
		}
	}
};
//
// File trailer for corrector.cpp
//
// [EOF]
//