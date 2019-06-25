#include <cmath>
#include<stdlib.h>
/// <summary>
/// 
/// </summary>
class predictor
{
public:
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="grav"></param>
	/// <param name="nf"></param>
	/// <param name="wse"></param>
	/// <param name="h"></param>
	/// <param name="u"></param>
	/// <param name="v"></param>
	/// <param name="dwsex"></param>
	/// <param name="dwsey"></param>
	/// <param name="dux"></param>
	/// <param name="duy"></param>
	/// <param name="dvx"></param>
	/// <param name="dvy"></param>
	/// <param name="dt2"></param>
	/// <param name="dzcx"></param>
	/// <param name="dzcy"></param>
	/// <param name="epsilon"></param>
	/// <param name="zc"></param>
	/// <param name="sox"></param>
	/// <param name="sfx"></param>
	/// <param name="dt"></param>
	/// <param name="soy"></param>
	/// <param name="sfy"></param>
	/// <param name="wsep"></param>
	/// <param name="up"></param>
	/// <param name="vp"></param>
	void fpredictor(int n, double grav, int nf, double** wse, double** h,
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

		for (int j = 0; j < nf; j++) {
			for (int k = 0; k < nf; k++) {
				//wse(j, k) - 0.5 * dt2 * (h(j, k) * dux(j, k) + u(j, k) * dhx(j, k) + h(j, k) * dvy(j, k) + v(j, k) * dhy(j, k));
				wsep[j][k] = wse[j][k] - 0.5 * dt2 * (h[j][k] * dux[j][k] + u[j][k] * dhx[j][k] + h[j][k] * dvy[j][k] + v[j][k] * dhy[j][k]);
				hp[j][k] = wsep[j][k] - zc[j][k];

				//  continuity for wet/dry interface (roger et al page 60)
				if (hp[j][k] < epsilon) {
					up[j][k] = 0.0;
					vp[j][k] = 0.0;
				}
				else {
					up[j][k] = u[j][k] - 0.5 * dt2 * (2.0 * u[j][k] *
						dux[j][k] + u[j][k] * dvy[j][k] + v[j][k] *
						duy[j][k] + grav * dhx[j][k]) - 0.5 * grav * dt * (sox[j][k] + sfx[j][k]);
					
					vp[j][k] = v[j][k] - 0.5 * dt2 * (u[j][k] *
						dvx[j][k] + v[j][k] * dux[j][k] + 2 * v[j][k] *
						dvy[j][k] + grav * dhy[j][k]) - 0.5 * grav * dt * (soy[j][k] + sfy[j][k]);
				}
			}
		}
	}
};

