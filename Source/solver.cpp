#include<cmath>
#include <string.h>
#include<algorithm>
#include<stdlib.h>
#include<iostream>
using namespace std;
class solver
{
public:
	///

	void fsolver(double hl, double hr, double ul, double ur, double vl, double vr,
		double sn, double cn, double hextra, double F[3], double* amax)
	{
		double grav = 0.0;
		double duml = pow(hl, 0.5);;
		double dumr = pow(hr, 0.5);
		double cl = pow(grav * hl, 0.5);
		double cr = pow(grav * hl, 0.5);
		double hhat = duml * dumr;
		double uhat;
		double vhat;
		double chat;
		double dh;
		double du;
		double dv;
		double dupar;
		double dW[3][3];
		double duperp;
		double uperpl;
		double uperpr;
		double a1;
		double da1;
		double a2;



		uhat = (duml * ul + dumr * ur) / ((duml + dumr) + hextra);
		vhat = (duml * vl + dumr * vr) / ((duml + dumr) + hextra);
		chat = pow(0.5 * grav * (hl + hr), 0.5) + hextra;
		dh = hr - hl;
		du = ur - ul;
		dv = vr - vl;


		dupar = du * sn + dv * cn;
		duperp = du * cn + dv * sn;


		dW[0][0] = 0.5 * (dh - hhat * duperp / chat);
		dW[0][1] = hhat * dupar;
		dW[0][2] = 0.5 * (dh + hhat * duperp / chat);
		dW[1][0] = 0;
		dW[1][1] = 0;
		dW[1][2] = 0;
		dW[2][0] = 0;
		dW[2][1] = 0;
		dW[2][2] = 0;


		uperpl = ul * cn + vl * sn;
		uperpr = ur * cn + vr * sn;

		double R[3][3];
		R[0][0] = 1;
		R[0][1] = 0;
		R[0][2] = 1;
		R[1][0] = uhat - chat * cn;
		R[1][1] = -sn;
		R[1][2] = uhat + chat * sn;
		R[1][0] = vhat - chat * sn;
		R[1][1] = cn;
		R[1][2] = vhat + chat * sn;

		//  a1=abs(uperp-chat);
		//  a2=abs(uperp);
		//  a3=abs(uperp+chat);
		//  Dry bed condition (HLLC by Toro,.. kim 2007..)
		double ustar, cstar, a3;

		ustar = (0.5 * (ul + ur) + cl) - cr;
		cstar = 0.5 * (cl + cr) + 0.25 * (ul - ur);

		if ((hl == 0.0) && (hr > 0.0)) {
			a1 = abs(ur - 2.0 * cr);
			a2 = abs(ur - 2.0 * cr);
			a3 = abs(ur + cr);
		}
		else if ((hr == 0.0) && (hl > 0.0)) {
			a1 = abs(ul - 2.0 * cl);
			a3 = abs(ul + cl);
			a2 = abs(ul + cl);
		}
		else {
			a1 = abs(minmax(ul - cl, ustar - cstar).first);
			a2 = abs(0.5 * (ul + ur) + cl - cr);
			a3 = abs(minmax(ur + cr, ustar + cstar).second);
		}
		double al1 = uperpl - cl;
		double al3 = uperpl + cl;
		double ar1 = uperpr - cr;
		double ar3 = uperpr + cr;
		da1 = 0.2 * (ar1 - al1) + hextra;
		double da3 = 0.2 * (ar3 - al3) + hextra;

		if (a1 < da1)
		{
			a1 = 0.5 * ((a1 * a1) / (da1 + da1));
		}
		if (a3 < da3)
		{
			a3 = 0.5 * ((a3 * a3) / (da3 + da3));
		}
		// Compute interface flux
		double A[3][3], FL[3][3], FR[3][3];
		//assign A
		A[0][0] = a1;
		A[0][1] = 0;
		A[0][2] = 0;
		A[1][0] = 0;
		A[1][1] = a2;
		A[1][2] = 0;
		A[2][0] = 0;
		A[2][1] = 0;
		A[2][2] = a3;

		//F
		FL[0][0] = uperpl * hl;
		FL[0][1] = (ul * uperpl * hl) + (0.5 * grav * hl * hl * cn);
		FL[0][2] = (vl * uperpl * hl) + (0.5 * grav * hl * hl * sn);
		FL[1][0] = 0;
		FL[1][1] = 0;
		FL[1][2] = 0;
		FL[2][0] = 0;
		FL[2][1] = 0;
		FL[2][2] = 0;
		//FR
		FR[0][0] = uperpr * hr;
		FR[0][1] = (ul * uperpr * hr) + (0.5 * grav * hr * hr * cn);
		FR[0][2] = (vl * uperpr * hr) + (0.5 * grav * hr * hr * sn);
		FR[1][0] = 0;
		FR[1][1] = 0;
		FR[1][2] = 0;
		FR[2][0] = 0;
		FR[2][1] = 0;
		FR[2][2] = 0;

		double FSUM[3][3];
		FSUM[0][0] = FL[0][0] + FR[0][0];
		FSUM[0][1] = FL[0][1] + FR[0][1];
		FSUM[0][2] = FL[0][2] + FR[0][2];
		FSUM[1][0] = 0;
		FSUM[1][1] = 0;
		FSUM[1][2] = 0;
		FSUM[2][0] = 0;
		FSUM[2][1] = 0;
		FSUM[2][2] = 0;

		//F = 0.5 * ((FL + FR- R * A * dW);
		double temp[3][3], sub[3][3];
		double** c = new double* [3];
		int n = 3;
		for (int i = 0; i < n; i++)
		{
			c[i] = new double[n];
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				c[i][j] = 0;
				for (int k = 0; k < n; k++)
				{
					c[i][j] += R[i][k] * A[k][j];
				}
			}
		}
		
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				temp[i][j] = 0;
				for (int k = 0; k < n; k++)
				{
				temp[i][j] += c[i][k] * dW[k][j];
				}
			}
		}
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				sub[i][j] = FSUM[i][j] - temp[i][j];
			}
	    }

		
		*amax = chat + abs(uhat * cn + vhat * sn);
	}
	
	
};
//
// File trailer for solver.cpp
//
// [EOF]
//