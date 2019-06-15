

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
	double*** fcorrector(double ***U, double ***F,  double ***G,
		double n, double dt2, double dt, double **sox, 
		double **sfx, double **soy, double **sfy,
		double grav)
	{
		int j;
		int k;
		for (j = 1; j < n - 1 ; j++) {
			for (k = 1; k < n - 1 ; k++) {
				U[0][j][k] = (U[0][j][k] - dt2 * (F[0][j][k - 1] - F[0][j][k])) - dt2 * (G[0][j][k - 1] - G[0][j][k]);
				if (U[0][j][k] < 0.0) {
					U[0][j][k] = 0.0;
				}

				U[1][j][k] = ((U[1][j][k] - dt2 * (F[1][j][k - 1] - F[1][j][k])) - dt2 * (G[1][j][k - 1] - G[1][j][k])) - dt * grav * (sox[j][k]
						+ sfx[j][k]);
				U[2][j][k] = ((U[2][j][k] - dt2 * (F[2][j][k - 1] - F[2][j][k]) - dt2 * G[2][j][k - 1] - G[2][j][k])) - dt * grav * (soy[j][k]
						+ sfy[j][k]);
			}
		}
		return U;
	}
	
};
//
// File trailer for corrector.cpp
//
// [EOF]
//