#include <cmath>
class predictor
{
public:
	predictor()
	{

	}
	//
	// Arguments    : double nf
	//                 double wse[1764]
	//                 double h[1764]
	//                 double u[1764]
	//                 double v[1764]
	//                 double dwsex[1764]
	//                 double dwsey[1764]
	//                 double dux[1764]
	//                 double duy[1764]
	//                 double dvx[1764]
	//                 double dvy[1764]
	//                double dt2
	//                 double dzcx[1764]
	//                 double dzcy[1764]
	//                double epsilon
	//                 double zc[1764]
	//                 double sox[1764]
	//                 double sfx[1764]
	//                double dt
	//                 double soy[1764]
	//                 double sfy[1764]
	//                double wsep[1764]
	//                double up[1764]
	//                double vp[1764]
	// Return Type  : void
	//
	void fpredictor(int n, double grav, double nf, double** wse, double** hdouble,
		double **  u, double ** v, double ** dwsex,
		double ** dwsey, double ** dux, double **
		duy, double ** dvx, double ** dvy, double 
		dt2, double ** dzcx, double ** dzcy, double 
		epsilon, double ** zc, double ** sox,
		double ** sfx, double dt, double ** soy, double **
		sfy, double ** wsep, double ** up, double ** vp)
	{
		
		double** dhx = (double**)malloc(sizeof(double*) * n * n);
		double **dhy = (double**)malloc(sizeof(double*) * n * n);
		double **hp = (double**)malloc(sizeof(double*) * n * n);
		for (int i = 0; i < n; i++)
		{
			dhx[i] = (double*)malloc(sizeof(double) * n * n);
			dhy[i] = (double*)malloc(sizeof(double) * n * n);
			hp[i] = (double*)malloc(sizeof(double) * n * n);
		}
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				dhx[i][j] = dwsex[i][j] - dzcx[i][j];
				dhy[i][j] = dwsey[i][j] - dzcy[i][j];
				hp[i][j] = 0.0;
				wsep[i][j] = 0.0;
				up[i][j] = 0.0;
				vp[i][j] = 0.0;
			}
		}

		for (int j = 0; j < (int)nf; j++) {
			for (int k = 0; k < (int)nf; k++) {
				wsep[j][k] = wse[j][k] - 0.5 * dt2 * (((hp[j][k] * dux[j][k] + u[j][k] * dhx[j][k]) + hp[j][k] * dvy[j][k]) 
					+ v[j][k] * dhy[j][k]);
				hp[j][k] = wsep[j][k] - zc[j][k];

				//  continuity for wet/dry interface (roger et al page 60)
				if (hp[j][k] < epsilon) {
					up[j][k] = 0.0;
					vp[j][k] = 0.0;
				}
				else {
					up[j][k] = (u[j][k] - 0.5 * dt2 * (2.0 * u[j][k] *
						dux[j][k] + u[j][k] * dvy[j][k] + v[j][k] *
						duy[j][k] + grav * dhx[j][k]) +0.5 * grav * dt * (sox[j][k] + sfx[j][k]));
					vp[j][k] = (v[j][k] - 0.5 * dt2 * (((u[j][k] * dvx[j][k] + v[j][k] * dux[j][k]) + 2.0 * v[j][k] *
						dvy[j][k]) + grav * dhy[j][k])) + 0.5 * grav * dt * (soy[j][k] + sfy[j][k]);

					// up(j,k) = u(j,k)-0.5*dt2*(2*u(j,k)*dux(j,k)+u(j,k)*dvy(j,k)+v(j,k)*duy(j,k)+grav*dhx(j,k))-0.5*grav*dt*(sox(j,k)+sfx(j,k)); 
					// vp(j,k) = v(j,k)-0.5*dt2*(u(j,k)*dvx(j,k)+v(j,k)*dux(j,k)+2*v(j,k)*dvy(j,k)+grav*dhy(j,k))-0.5*grav*dt*(soy(j,k)+sfy(j,k)); 
				}
			}
		}
	}
};
//
// File trailer for predictor.cpp
//
// [EOF]
//
