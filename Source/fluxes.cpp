#include "solver.cpp"
#include <stdlib.h>
using namespace std;
class fluxes
{

public:
	// Function Definitions

//
// Compute fluxes in X-direction
// Arguments    : const double UP[5292]
//                double n
//                const double dwsex[1764]
//                const double dwsey[1764]
//                const double dux[1764]
//                const double duy[1764]
//                const double dvx[1764]
//                const double dvy[1764]
//                double hextra
//                const double zc[1764]
//                emxArray_real_T *F
//                emxArray_real_T *G
//                double amax
// Return Type  : void
//
	void ffluxes(double*** UP, int n, double** dwsex,
		double** dwsey, double** dux, double** duy,
		double** dvx, double** dvy, double hextra,
		double** zc, double*** F, double*** G, double  &amax)
	{
		double hr=0.0;
		double hl=0.0;
		double ul=0.0;
		double vl=0.0;
		double ur=0.0;
		double vr=0.0;
		amax = 0.0;
		

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					F[i][j][k] = 0.0;
					G[i][j][k] = 0.0;
				}
			}
		}

		ul = 0.0;
		//  added
		vl = 0.0;

		//  added
		//  Enforce wall boundary on left side of box (west)
		amax =boundaryWest(n,UP,F,amax,zc,hextra);
		//  Enforce wall boundary on second rows along the X-axis
		amax=boundaryX(amax, UP, hextra, zc, F, n);
		// Compute the fluxe in the X-direction on the domain
		amax = xDirectionFlux(zc, UP, amax, F, dwsex,hextra,n,dux,dvx);
		// From 20th to 21st and 22nd rows
		amax = middleX(zc, UP, amax, F, dwsex, hextra, n, dux, dvx);
		//  % For the 23rd row
		amax = twentythirdX(zc, UP, amax, F, dwsex, hextra, n, dux, dvx);
		//boundary east
		amax = boundaryEast(n, UP, F, amax, zc,hextra);
	/**************************Y Direction***********************/
		//boundary south Y direction
		amax = boundarySouthY(n,UP,G,amax,zc,hextra);
		//second column y
		amax = secondColumnY(n, UP, G, amax, zc, hextra);
		//Y fluxes
		amax = yDirectionFlux(zc, UP, amax, G, dwsex, hextra, n, duy, dwsey, dvy);
		//23 
		amax =twentythreerowDownStream(zc, UP, amax, G, hextra, n);
		//24
		amax = twentyfourDam(zc,UP,amax,G,hextra,n);
		//35
		amax=thirtyFiveLeftDamn(zc,UP,amax,G,hextra,n);
		//left dam
		amax = leftDamn(zc, UP, amax, G, hextra, n);
	}
	double xDirectionFlux(double **zc,double ***UP,double amax,double ***F,double **dwsex,double hextra,int n,double **dux,
		double **dvx)
	{
		double zbc = 0.0;
		double hr = 0.0;
		double ur = 0.0;
		double vr = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
//					zbc = minmax(zc[j - 1][k]), zc[j],[k]).second();

					if ((UP[0][j][k] - zbc) > 0)
					{
					hl = UP[0][j][k] - zbc;
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hl = 0;
					}
					
					if ((UP[0][j][k] - zbc) > 0)
					{
						hr = UP[0][j][k] - zbc;
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hr = 0;
					}

					if (hr < 0.0) {
						hr = 0.0;
					}
					if (hl > 0)
					{
						hl = hl + 0.5 * dwsex[j][k];
						if (hl < 0)
						{
							hl = 0;
						}
					}
					if (hl == 0.0) {
						ul = 0.0;
						vl = 0.0;
					}
					else {
						ul = (UP[1][j][k] / (hr + hextra)) + 0.5 * dux[j][k];
						vl = (UP[2][j][k] / (hr + hextra)) + 0.5 * dvx[j][k];
					}
					if (hr < 0.0) {
						hr = 0.0;
					}
					if (hr > 0.0)
					{
						hr = hr - 0.5 * dwsex[j][k];
						if (hr < 0)
						{
							hr = 0;
						}
					}

					if (hr == 0.0) {
						ur = 0.0;
						vr = 0.0;
					}
					else {
						ur = (UP[1][j][k] / (hr + hextra)) - 0.5 * dux[j][k];
						vr = (UP[2][j][k] / (hr + hextra)) - 0.5 * dvx[j][k];
					}
					s.fsolver(hl, hr, ul, ur, vl, vr, 0.0, 1.0, hextra, dv0, &zbc);


					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							F[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}

			}

		}
		
		return amax;
	}
	// From 20th to 21st and 22nd rows
	double middleX(double** zc, double*** UP, double amax, double*** F, double** dwsex, double hextra, int n, double** dux,
		double** dvx)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 20; j < 21; j++)
			{
				for (int k = 1; k < 23 || k > 34 && k < n; k++)
				{
					//zbc = max(zc(j, k), zc(j, k));
					if ((UP[0][j][k] - zbc) > 0)
					{
						hl = UP[0][j][k] - zbc;

						if (hl < 0)
						{
							hl = 0;
						}
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[1][j][k]) / (hl + hextra);
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hl = 0; ul = 0; vl = 0;
					}
							
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
					
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							F[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}

			}

		}

		return amax;
	}
	//23 X direction
	double twentythirdX(double** zc, double*** UP, double amax, double*** F, double** dwsex, double hextra, int n, double** dux,
		double** dvx)
	{
		double zbc = 0.0;
		double hr = 0.0;
		double ur = 0.0;
		double vr = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 23; j < 24; j++)
			{
				for (int k = 1; (k < 24 || k >34) && k < n; k++)
				{

					//zbc = max(zc(j, k), zc(j, k));
					if ((UP[0][j][k] - zbc) > 0)
					{
						hl = UP[0][j][k] - zbc;
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hl = 0;
					}
					if ((UP[0][j][k] - zbc) > 0)
					{
						hr = UP[0][j][k] - zbc;
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hr = 0;
					}

					if (hl > 0)
					{
						hl = hl + 0.0;
						if (hl < 0)
							hl = 0;
					}

					if (hr > 0)
					{
						hr = hr - 0.5 * (UP[0][j][k] - UP[0][j][k]); //% wse difference bn 25 & 26
						if (hr < 0)
							hr = 0;
					}

					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					else
					{
						ul = ((UP[1][j][k]) / (hl + hextra));
						vl = ((UP[2][j][k]) / (hl + hextra));
					}
					if (hr == 0)
					{
						ur = 0; vr = 0;
					}
					else
					{
						ur = ((UP[1][j][k]) / (hr + hextra));
						vr = ((UP[2][j][k]) / (hr + hextra));
					}

					s.fsolver(hl, hr, ul, ur, vl, vr, 0.0, 1.0, hextra, dv0, &zbc);
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							F[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;
		}
		double boundaryX(double amax, double*** UP, double hextra, double** zc, double*** F,int n)
	{
		solver s;
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		
				for (int k = 0; k < n; k++)
				{
					hl = UP[0][2][k] - zc[2][k];
					if (hl < 0.0) {
						hl = 0.0;
					}

					if (hl == 0.0) {
						ul = 0.0;
						vl = 0.0;
					}
					else {
						ul = UP[1][1][k] / (hl + hextra);
						vl = UP[2][1][k] / (hl + hextra);
					}
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
				
						for (int j = 0; j < 3; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							F[2][0][j] = dv0[j];
						}
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
					
		return amax;
	}
	double boundaryWest(int n,double ***UP, double ***F,double &amax,double **zc,double hextra)
	{	
		double zbc = 0.0;
		double hr = 0.0;
		double ur = 0.0;
		double vr = 0.0;
		double dv0[3];
		solver s;
			for (int k = 1; k < n-1; k++)
				{
					hr = UP[0][1][k] - zc[1][k];
					if (hr < 0.0) {
						hr = 0.0;
					}

					if (hr == 0.0) {
						ur = 0.0;
						vr = 0.0;
					}
					else {
						ur = UP[1][1][k] / (hr + hextra);
						vr = UP[2][1][k] / (hr + hextra);
					}
					s.fsolver(hr, hr, ur, -ur, vr, vr, 0.0, 1.0, hextra, dv0, &zbc);
									
						for (int j = 0; j < 3; j++)
						{
							F[2][0][j] = dv0[j];
						}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
								
			return amax;
	}
	double boundaryNorth(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					// Enforce wall boundary on top of box (north)
				
						hl = UP[0][j][k]- zc[j][k];
						if (hl < 0.0) {
							hl = 0.0;
						}

						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
						s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
						for (int i = 0; i < n; i++) {
							for (int j = 0; j < n; j++)
							{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							G[0][i][j] = dv0[i];
							}
					    }
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;
	}
	double boundaryEast(int n, double*** UP, double*** F, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					hl = UP[0][j][k] - zc[j][k];
					if (hl < 0)
					{
						hl = 0;
					}
					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					else
					{
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
					}
					//s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 0.05, hextra);
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							F[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;
	}
	double boundarySouthY(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					//hl = UP(j, 2, 1) - zc(j, 2);
					hl = UP[0][j][k] - zc[j][k];
					if (hl < 0)
					{
						hl = 0;
					}
					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					else
					{
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
					}
					//s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 0.05, hextra);
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0, hextra, dv0, &zbc);
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;


	}
	double secondColumnY(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					//hl = UP(j, 2, 1) - zc(j, 2);
					hl = UP[0][j][k] - zc[j][k];
					if (hl < 0)
					{
						hl = 0;
					}
					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					else
					{
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
					}
					//s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 0.05, hextra);
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0, hextra, dv0, &zbc);
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;

	}
	double yDirectionFlux(double** zc, double*** UP, double amax, double*** G, double** dwsex, double hextra, int n, 
		double** duy,	double** dwsey,double ** dvy)
	{
		double zbc = 0.0;
		double hr = 0.0;
		double ur = 0.0;
		double vr = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
		
				if ((zc[j][k-1] > zc[j][k]) || !isnan(zc[j][k])) {
					zbc = zc[j][k];
				}
				else {
					zbc = zc[j][k];
				}

				if ((UP[0][j][k] - zbc) > 0.0) {
					hl = (UP[0][j][k] - zbc);
				}
				else {
					if ((UP[0][j][k] - zbc) <= 0.0)
					{
						hl = 0.0;
					}
				}

				if ((UP[0][j][k] - zbc) > 0.0) {
					hr = UP[0][j][k] - zbc;
				}
				else {
					if ((UP[0][j][k] - zbc) <= 0.0) {
						hr = 0.0;
					}
				}

				if (hl > 0.0) {
					hl += 0.5 * dwsey[j][k];
					if (hl < 0.0) {
						hl = 0.0;
					}
				}

				if (hr > 0.0) {
					hr -= 0.5 * dwsey[j][k];
					if (hr < 0.0) {
						hr = 0.0;
					}
				}

				if (hl == 0.0) {
					ul = 0.0;
					vl = 0.0;
				}
				else {
					ul = UP[1][j][k-1] + 0.5 * (duy[j][k-2] /
						(hl + hextra));
					vl = UP[2][j][k - 2] + 0.5 * (duy[j][k - 2] /
						(hl + hextra));
				}

				if (hr == 0.0) {
					ur = 0.0;
					vr = 0.0;
				}
				else {
					ur = UP[1][j][k] / (hr + hextra) - 0.5 * duy[j][k];
					vr = UP[2][j][k] / (hr + hextra) - 0.5 * dvy[j][k];
				}

				s.fsolver(hl, hr, ul, ur, vl, vr, 1.0, 0.0, hextra, dv0, &zbc);
				for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							//F->data[F->size[0] * k + F->size[0] * F->size[1] * r0[i0]] = dv0[i0];
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
			
		}
		}
		}
		return amax;
	}
	double twentythreerowDownStream(double** zc, double*** UP, double amax, double*** G,double hextra,int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 1; j < 23; j++)
			{
				for (int k = 21; k < 22; k++)
				{
					
								hl = UP[0][j][k] - zc[j][k];
								if (hl < 0.0) {
									hl = 0.0;
								}

								ul = UP[1][j][k] / (hl + hextra);
								vl = UP[2][j][k] / (hl + hextra);
								s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
								
								for (int i = 0; i < n; i++) {
									for (int j = 0; j < n; j++)
									{
										G[0][i][j] = dv0[i];
									}
								}
								if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
								}
								else {
									amax = zbc;
								}
							}
						
				}
			}
			return amax;
	}
	double twentyfourDam(double** zc, double*** UP, double amax, double*** G,  double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 24; j < 24; j++)
			{
				for (int k = 21; k < 22; k++)
				{

					//zbc = max(zc(j, k), zc(j, k));
					if ((UP[0][j][k] - zbc) > 0)
					{
						hl = UP[0][j][k] - zbc;
						if (hl < 0)
							hl = 0;
						
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hl = 0;
					}
					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					else
					{
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
					}
					
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);

					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;
	}
	double thirtyFiveLeftDamn(double** zc, double*** UP, double amax, double*** G, double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 35; j < 36; j++)
			{
				for (int k = 21; k < 22; k++)
				{

					//zbc = max(zc(j, k), zc(j, k));
					if ((UP[0][j][k] - zbc) > 0)
					{
						hl = UP[0][j][k] - zbc;
						if (hl < 0)
						{
							hl = 0;
						}
						ul = (UP[1][j][k]) / (hl + hextra);
						vl = (UP[2][j][k]) / (hl + hextra);
					}
					else if ((UP[0][j][k] - zbc) <= 0)
					{
						hl = 0;
					}
					if (hl == 0)
					{
						ul = 0; vl = 0;
					}
					
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);

					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			}
		}
		return amax;
	}
	double leftDamn(double** zc, double*** UP, double amax, double*** G, double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
		
					hl = UP[0][j][k] - zc[j][k];
					
				if (hl < 0.0) {
					hl = 0.0;
				}

				ul = UP[1][j][k] / (hl + hextra);
				vl = UP[2][j][k] / (hl + hextra);
				s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
				
				for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++)
						{
							G[0][i][j] = dv0[i];
						}
					}
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
	}
	}
	}
		return amax;
	}
	


	
		//
	// File trailer for fluxes.cpp
	//
	// [EOF]
	//
};