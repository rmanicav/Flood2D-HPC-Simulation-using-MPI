#include<cmath>
#include <string.h>
#include<algorithm>
#include<stdlib.h>
#include<iostream>
using namespace std;
/// <summary>
/// 
/// </summary>
class solver
{
public:

	/// <summary>
	/// 
	/// </summary>
	/// <param name="hl"></param>
	/// <param name="hr"></param>
	/// <param name="ul"></param>
	/// <param name="ur"></param>
	/// <param name="vl"></param>
	/// <param name="vr"></param>
	/// <param name="sn"></param>
	/// <param name="cn"></param>
	/// <param name="hextra"></param>
	/// <param name="F"></param>
	/// <param name="amax"></param>
	void fsolver(double hl, double hr, double ul, double ur, double vl, double vr,
		double sn, double cn, double hextra, double(&a)[3], double* amax)
	{
		double grav = 9.806;
		double duml = pow(hl, 0.5);
		double dumr = pow(hr, 0.5);
		double cl = pow((grav * hl), 0.5);
		double cr = pow((grav * hr), 0.5);
		double hhat = duml * dumr;
		double uhat;
		double vhat;
		double chat;
		double dh;
		double du;
		double dv;
		double dupar;
		double dW[3];
		double duperp;
		double uperpl;
		double uperpr;
		double a1;
		double da1;
		double a2;



		uhat = ((duml * ul) + (dumr * ur)) / ((duml + dumr) + hextra);
		vhat = ((duml * vl) + (dumr * vr)) / ((duml + dumr) + hextra);
		chat = pow((0.5 * grav * (hl + hr)), 0.5) + (hextra);
		dh = hr - hl;
		du = ur - ul;
		dv = vr - vl;


		dupar = (dv * cn) - (du * sn);
		duperp = (du * cn) + (dv * sn);


		dW[0] = 0.5 * (dh - ((hhat * duperp) / chat));
		dW[1] = hhat * dupar;
		dW[2] = 0.5 * (dh + ((hhat * duperp) / chat));



		uperpl = (ul * cn) + (vl * sn);
		uperpr = (ur * cn) + (vr * sn);

		double R[3][3];
		R[0][0] = 1;
		R[0][1] = 0;
		R[0][2] = 1;
		R[1][0] = (uhat)-(chat * cn);
		R[1][1] = -sn;
		R[1][2] = (uhat)+(chat * cn);
		R[2][0] = (vhat)-(chat * sn);
		R[2][1] = cn;
		R[2][2] = (vhat)+(chat * sn);

		//  a1=abs(uperp-chat);
		//  a2=abs(uperp);
		//  a3=abs(uperp+chat);
		//  Dry bed condition (HLLC by Toro,.. kim 2007..)
		double ustar, cstar, a3;

		ustar = ((0.5 * (ul + ur)) + cl) - cr;
		cstar = (0.5 * (cl + cr)) + (0.25 * (ul - ur));

		if ((hl == 0.0) && (hr > 0.0)) {
			a1 = abs(ur - 2.0 * cr);
			a2 = abs(ur - 2.0 * cr);
			a3 = abs(ur + cr);
		}
		else if ((hr == 0.0) && (hl > 0.0)) {
			a1 = abs(ul - (2.0 * cl));
			a3 = abs(ul + cl);
			a2 = abs(ul + cl);
		}
		else {
			a1 = abs(minmax(ul - cl, ustar - cstar).first);
			a2 = abs(((0.5 * (ul + ur)) + cl) - cr);
			a3 = abs(minmax(ur + cr, ustar + cstar).second);
		}
		double al1 = uperpl - cl;
		double al3 = uperpl + cl;
		double ar1 = uperpr - cr;
		double ar3 = uperpr + cr;
		// da1 = 0.2 * (ar1 - al1) + hextra;
		//double da3 = 0.2 * (ar3 - al3) + hextra;
		da1 = minmax(0.0, (2 * (ar1 - al1))).second + (hextra);
		double da3 = minmax(0.0, (2 * (ar3 - al3))).second + (hextra);

		if (a1 < da1)
		{
			a1 = 0.5 * ((a1 * a1) / (da1 + da1));
		}
		if (a3 < da3)
		{
			a3 = 0.5 * ((a3 * a3) / (da3 + da3));
		}
		// Compute interface flux
		double A[3][3], FL[3], FR[3], FSUM[3];
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
		FL[0] = uperpl * hl;
		FL[1] = (ul * uperpl * hl) + (0.5 * grav * hl * hl * cn);
		FL[2] = (vl * uperpl * hl) + (0.5 * grav * hl * hl * sn);

		//FR
		FR[0] = uperpr * hr;
		FR[1] = (ur * uperpr * hr) + (0.5 * grav * hr * hr * cn);
		FR[2] = (vr * uperpr * hr) + (0.5 * grav * hr * hr * sn);


		FSUM[0] = FL[0] + FR[0];
		FSUM[1] = FL[1] + FR[1];
		FSUM[2] = FL[2] + FR[2];




		//F = 0.5 * ((FL + FR- R * A * dW);

		double temp[3];
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
			temp[i] = 0;
			for (int k = 0; k < n; k++)
			{
				temp[i] += c[i][k] * dW[k];
			}
		}

		a[0] = 0.5 * (FSUM[0] - temp[0]);
		a[1] = 0.5 * (FSUM[1] - temp[1]);
		a[2] = 0.5 * (FSUM[2] - temp[2]);

		*amax = chat + abs((uhat * cn) + (vhat * sn));
	}

};
