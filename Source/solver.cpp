#include<cmath>
#include <string.h>
class solver
{
public:
	// Function Definitions

	//
	// Arguments    : double hl
	//                double hr
	//                double ul
	//                double ur
	//                double vl
	//                double vr
	//                double sn
	//                double cn
	//                double hextra
	//                double F[3]
	//                double *amax
	// Return Type  : void
	//
	void fsolver(double hl, double hr, double ul, double ur, double vl, double vr,
		double sn, double cn, double hextra, double F[3], double* amax)
	{
		double duml;
		double dumr;
		double cl;
		double cr;
		double hhat;
		double uhat;
		double vhat;
		double chat;
		double dh;
		double du;
		double dv;
		double duperp;
		double uperpl;
		double uperpr;
		double a1;
		double da1;
		double a2;
		double varargin_2;
		double b_da1;
		double d0;
		double c_da1;
		double d1;
		double v[3];
		double A[9];
		double b_uperpl[3];
		double b_uperpr[3];
		int j;
		double dv0[9];
		static const signed char iv0[3] = { 1, 0, 1 };

		int i0;
		double dv1[9];
		int i1;
		double grav = 0.0;

		// Compute Roe averages
		duml = std::sqrt(hl);
		dumr = std::sqrt(hr);
		cl = std::sqrt(grav * hl);
		cr = std::sqrt(grav * hr);
		hhat = duml * dumr;
		uhat = (duml * ul + dumr * ur) / ((duml + dumr) + hextra);
		vhat = (duml * vl + dumr * vr) / ((duml + dumr) + hextra);
		chat = std::sqrt(0.5 * grav * (hl + hr)) + hextra;
		dh = hr - hl;
		du = ur - ul;
		dv = vr - vl;
		duperp = du * cn + dv * sn;
		uperpl = ul * cn + vl * sn;
		uperpr = ur * cn + vr * sn;

		//  a1=abs(uperp-chat);
		//  a2=abs(uperp);
		//  a3=abs(uperp+chat);
		//  Dry bed condition (HLLC by Toro,.. kim 2007..)
		duml = (0.5 * (ul + ur) + cl) - cr;
		dumr = 0.5 * (cl + cr) + 0.25 * (ul - ur);
		if ((hl == 0.0) && (hr > 0.0)) {
			a1 = std::abs(ur - 2.0 * cr);
			a2 = std::abs(ur - 2.0 * cr);
			dumr = std::abs(ur + cr);
		}
		else if ((hr == 0.0) && (hl > 0.0)) {
			a1 = std::abs(ul - 2.0 * cl);
			dumr = std::abs(ul + cl);
			a2 = std::abs(ul + cl);
		}
		else {
			da1 = ul - cl;
			varargin_2 = duml - dumr;
			if ((da1 < varargin_2) || isnan(varargin_2)) {
				b_da1 = da1;
			}
			else {
				b_da1 = varargin_2;
			}

			a1 = std::abs(b_da1);
			a2 = std::abs((0.5 * (ul + ur) + cl) - cr);
			da1 = ur + cr;
			varargin_2 = duml + dumr;
			if ((da1 > varargin_2) || isnan(varargin_2)) {
				c_da1 = da1;
			}
			else {
				c_da1 = varargin_2;
			}

			dumr = std::abs(c_da1);
		}

		duml = 2.0 * ((uperpr - cr) - (uperpl - cl));
		if (0.0 < duml) {
			d0 = duml;
		}
		else {
			d0 = 0.0;
		}

		da1 = d0 + hextra;
		duml = 2.0 * ((uperpr + cr) - (uperpl + cl));
		if (0.0 < duml) {
			d1 = duml;
		}
		else {
			d1 = 0.0;
		}

		duml = d1 + hextra;

		// Critical flow fix
		if (a1 < da1) {
			a1 = 0.5 * (a1 * a1 / da1 + da1);
		}

		if (dumr < duml) {
			dumr = 0.5 * (dumr * dumr / duml + duml);
		}

		// Compute interface flux
		v[0] = a1;
		v[1] = a2;
		v[2] = dumr;
		memset(&A[0], 0, 9U * sizeof(double));
		b_uperpl[0] = uperpl * hl;
		b_uperpl[1] = ul * uperpl * hl + 0.5 * grav * hl * hl * cn;
		b_uperpl[2] = vl * uperpl * hl + 0.5 * grav * hl * hl * sn;
		b_uperpr[0] = uperpr * hr;
		b_uperpr[1] = ur * uperpr * hr + 0.5 * grav * hr * hr * cn;
		b_uperpr[2] = vr * uperpr * hr + 0.5 * grav * hr * hr * sn;
		for (j = 0; j < 3; j++) {
			A[j + 3 * j] = v[j];
			dv0[3 * j] = iv0[j];
		}

		dv0[1] = uhat - chat * cn;
		dv0[4] = -sn;
		dv0[7] = uhat + chat * cn;
		dv0[2] = vhat - chat * sn;
		dv0[5] = cn;
		dv0[8] = vhat + chat * sn;
		v[0] = 0.5 * (dh - hhat * duperp / chat);
		v[1] = hhat * (-du * sn + dv * cn);
		v[2] = 0.5 * (dh + hhat * duperp / chat);
		for (j = 0; j < 3; j++) {
			duml = 0.0;
			for (i0 = 0; i0 < 3; i0++) {
				dv1[j + 3 * i0] = 0.0;
				for (i1 = 0; i1 < 3; i1++) {
					dv1[j + 3 * i0] += dv0[j + 3 * i1] * A[i1 + 3 * i0];
				}

				duml += dv1[j + 3 * i0] * v[i0];
			}

			F[j] = 0.5 * ((b_uperpl[j] + b_uperpr[j]) - duml);
		}

		*amax = chat + std::abs(uhat * cn + vhat * sn);
	}
};
//
// File trailer for solver.cpp
//
// [EOF]
//