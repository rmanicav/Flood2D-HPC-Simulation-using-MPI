#include "../Include/predictor.h"

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
void fpredictor(double nf,  double wse[1764],  double h[1764], 
	double u[1764],  double v[1764],  double dwsex[1764],
	 double dwsey[1764],  double dux[1764],  double
	duy[1764],  double dvx[1764],  double dvy[1764], double
	dt2,  double dzcx[1764],  double dzcy[1764], double
	epsilon,  double zc[1764],  double sox[1764], 
	double sfx[1764], double dt,  double soy[1764],  double
	sfy[1764], double wsep[1764], double up[1764], double vp[1764])
{
	double grav;
	int j;
	double dhx[1764];
	double dhy[1764];
	double hp[1764];
	int k;
	grav = 0.0;
	for (j = 0; j < 1764; j++) {
		dhx[j] = dwsex[j] - dzcx[j];
		dhy[j] = dwsey[j] - dzcy[j];
		hp[j] = 0.0;
		wsep[j] = 0.0;
		up[j] = 0.0;
		vp[j] = 0.0;
	}

	for (j = 0; j < (int)nf; j++) {
		for (k = 0; k < (int)nf; k++) {
			wsep[j + 42 * k] = wse[j + 42 * k] - 0.5 * dt2 * (((h[j + 42 * k] * dux[j
				+ 42 * k] + u[j + 42 * k] * dhx[j + 42 * k]) + h[j + 42 * k] * dvy[j +
				42 * k]) + v[j + 42 * k] * dhy[j + 42 * k]);
			hp[j + 42 * k] = wsep[j + 42 * k] - zc[j + 42 * k];

			//  continuity for wet/dry interface (roger et al page 60)
			if (hp[j + 42 * k] < epsilon) {
				up[j + 42 * k] = 0.0;
				vp[j + 42 * k] = 0.0;
			}
			else {
				up[j + 42 * k] = (u[j + 42 * k] - 0.5 * dt2 * (((2.0 * u[j + 42 * k] *
					dux[j + 42 * k] + u[j + 42 * k] * dvy[j + 42 * k]) + v[j + 42 * k] *
					duy[j + 42 * k]) + grav * dhx[j + 42 * k])) + 0.5 * grav * dt * (sox[j
						+ 42 * k] + sfx[j + 42 * k]);
				vp[j + 42 * k] = (v[j + 42 * k] - 0.5 * dt2 * (((u[j + 42 * k] * dvx[j +
					42 * k] + v[j + 42 * k] * dux[j + 42 * k]) + 2.0 * v[j + 42 * k] *
					dvy[j + 42 * k]) + grav * dhy[j + 42 * k])) + 0.5 * grav * dt * (soy[j
						+ 42 * k] + sfy[j + 42 * k]);

				// up(j,k) = u(j,k)-0.5*dt2*(2*u(j,k)*dux(j,k)+u(j,k)*dvy(j,k)+v(j,k)*duy(j,k)+grav*dhx(j,k))-0.5*grav*dt*(sox(j,k)+sfx(j,k)); 
				// vp(j,k) = v(j,k)-0.5*dt2*(u(j,k)*dvx(j,k)+v(j,k)*dux(j,k)+2*v(j,k)*dvy(j,k)+grav*dhy(j,k))-0.5*grav*dt*(soy(j,k)+sfy(j,k)); 
			}
		}
	}
}

//
// File trailer for predictor.cpp
//
// [EOF]
//
