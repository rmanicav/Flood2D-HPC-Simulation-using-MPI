#include "solver.cpp"
#include<algorithm>
#include<fstream>

using namespace std;
/// <summary>
/// 
/// </summary>
class fluxes
{

public:
	/// <summary>
	/// 
	/// </summary>
	/// <param name="UP"></param>
	/// <param name="n"></param>
	/// <param name="dwsex"></param>
	/// <param name="dwsey"></param>
	/// <param name="dux"></param>
	/// <param name="duy"></param>
	/// <param name="dvx"></param>
	/// <param name="dvy"></param>
	/// <param name="hextra"></param>
	/// <param name="zc"></param>
	/// <param name="F"></param>
	/// <param name="G"></param>
	/// <param name="amax"></param>
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

		//write3dOutputFile(UP, n, "UP.txt");
		const int out = remove("F.txt");
		
		//  added
		//  Enforce wall boundary on left side of box (west)
		amax =boundaryWest(n,UP,F,amax,zc,hextra);
		cout << endl;
		
	//	print3dArray(F, n, "F");
		//  Enforce wall boundary on second rows along the X-axis
		amax=boundaryX(amax, UP, hextra, zc, F, n);
		//	print3dArray(F, n, "F");
		
		
		
		// Compute the fluxe in the X-direction on the domain
		amax = xDirectionFlux(zc, UP, amax, F, dwsex,hextra,n,dux,dvx);
		//print3dArray(F, n, "F");
		
		
		
		// From 20th to 21st and 22nd rows
		amax = middleX(zc, UP, amax, F, dwsex, hextra, n, dux, dvx);
			//print3dArray(F, n, "F");
		

		write3dOutputFile(F, n, "F.txt");
		//  % For the 23rd row
		amax = twentythirdX(zc, UP, amax, F, dwsex, hextra, n, dux, dvx);
		//print3dArray(F, n, "F");
		//boundary east
		
	  
		amax = boundaryEast(n, UP, F, amax, zc,hextra);
		
		
	//print3dArray(F, n, "F");
		
	/**************************Y Direction***********************/
		//boundary south Y direction
		amax = boundarySouthY(n,UP,G,amax,zc,hextra);
		
		//print3dArray(G, n, "G");
		//second column y
		amax = secondColumnY(n, UP, G, amax, zc, hextra);
		//write3dOutputFile(G, n, "G.txt");
		//Y fluxes
		amax = yDirectionFlux(zc, UP, amax, G, dwsex, hextra, n, duy, dwsey, dvy);
	//	print3dArray(G, n, "G");
		write3dOutputFile(G, n, "G.txt");
		//23 
	/*	amax =twentythreerowDownStream(zc, UP, amax, G, hextra, n);
		//print3dArray(G, n, "G");
		//24
		amax = twentyfourDam(zc,UP,amax,G,hextra,n);
		//print3dArray(G, n, "G");
		//35
		amax=thirtyFiveLeftDamn(zc,UP,amax,G,hextra,n);
		//print3dArray(G, n, "G");
		//left dam
		amax = leftDamn(zc, UP, amax, G, hextra, n);
		//print3dArray(G, n, "G");
		//boundary North direction
		amax = boundaryNorth(n, UP, G, amax, zc, hextra);
	//	print3dArray(G, n, "G");
		//write3dOutputFile(G, n, "G.txt");*/
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="F"></param>
	/// <param name="dwsex"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <param name="dux"></param>
	/// <param name="dvx"></param>
	/// <returns></returns>
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
		clearArray(dv0);

			for (int j = 2; j < n - 1; j++)
			{
				for (int k = 0; k < n ; k++)
				{
					zbc = minmax(zc[j - 1][k], zc[j][k]).second;

					if ((UP[0][j-1][k] - zbc) > 0)
					{
					hl = UP[0][j-1][k] - zbc;
					}
					else if ((UP[0][j-1][k] - zbc) <= 0)
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
						hl = hl + 0.5 * dwsex[j-1][k];
						if (hl < 0)
						{
							hl = 0;
						}
					}
					if (hr > 0.0)
					{
						hr = hr - 0.5 * dwsex[j][k];
						if (hr < 0)
						{
							hr = 0;
						}
					}

					if (hl == 0.0) {
						ul = 0.0;
						vl = 0.0;
					}
					else {
						ul = (UP[1][j-1][k] / (hl + hextra)) + 0.5 * dux[j-1][k];
						vl = (UP[2][j-1][k] / (hl + hextra)) + 0.5 * dvx[j-1][k];
					}
					if (hr < 0.0) {
						hr = 0.0;
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
					F[0][j][k] = dv0[0];
					F[1][j][k] = dv0[1];
					F[2][j][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}

			}

		
		
		return amax;
	}
	// From 20th to 21st and 22nd rows
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="F"></param>
	/// <param name="dwsex"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <param name="dux"></param>
	/// <param name="dvx"></param>
	/// <returns></returns>
	double middleX(double** zc, double*** UP, double amax, double*** F, double** dwsex, double hextra, int n, double** dux,
		double** dvx)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;
		clearArray(dv0);
		
			for (int j = 19; j < 21; j++)
			{
				for (int k = 0;k<23; k++)
				{
					zbc = minmax(zc[j][k], zc[j][k]).second;
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
						hl = 0; ul = 0; vl = 0;
					}
							
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
					
					
					F[0][j+1][k] = dv0[0];
					F[1][j+1][k] = dv0[1];
					F[2][j+1][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}

			}

			//columnx 35 - n-1
			for (int j = 19; j < 21; j++)
			{
				for (int k = 34; k < n - 1; k++)//////////////////n-1
				{
					zbc = minmax(zc[j][k], zc[j][k]).second;
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
						hl = 0; ul = 0; vl = 0;
					}

					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);


					F[0][j+1][k] = dv0[0];
					F[1][j+1][k] = dv0[1];
					F[2][j+1][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
				
			}

		

		return amax;
	}
	//23 X direction
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="F"></param>
	/// <param name="dwsex"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <param name="dux"></param>
	/// <param name="dvx"></param>
	/// <returns></returns>
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
		clearArray(dv0);
		
		int j = 22;
				for (int k = 0; k < 23 ; k++)
				{
					zbc = minmax(zc[j][k], zc[j][k]).second;
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
					
					F[0][j][k] = dv0[0];
					F[1][j][k] = dv0[1];
					F[2][j][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			
			//columns 35 - n
			
				for (int k = 34; k < n; k++)
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

					F[0][j][k] = dv0[0];
					F[1][j][k] = dv0[1];
					F[2][j][k] = dv0[2];

					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			
		
		return amax;
		}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="amax"></param>
	/// <param name="UP"></param>
	/// <param name="hextra"></param>
	/// <param name="zc"></param>
	/// <param name="F"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	double boundaryX(double amax, double*** UP, double hextra, double** zc, double*** F,int n)
	{
		solver s;
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		
				for (int k = 1; k < n - 1; k++)
				{
					hl = UP[0][1][k] - zc[1][k];
					if (hl < 0.0) {
						hl = 0.0;
					}

					if (hl == 0.0) {
						ul = 0.0;
						vl = 0.0;
					}
					else {
						ul = UP[1][0][k] / (hl + hextra);
						vl = UP[2][0][k] / (hl + hextra);
					}
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
				
					F[0][1][k] = dv0[0];
					F[1][1][k] = dv0[1];
					F[2][1][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
					
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="UP"></param>
	/// <param name="F"></param>
	/// <param name="amax"></param>
	/// <param name="zc"></param>
	/// <param name="hextra"></param>
	/// <returns></returns>
	double boundaryWest(int n,double ***UP, double ***F,double &amax,double **zc,double hextra)
	{	
		double zbc = 0.0;
		double hr = 0.0;
		double ur = 0.0;
		double vr = 0.0;
		double dv0[3];
		solver s;
			for (int k = 1; k < n - 1; k++)
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
					s.fsolver(hr, hr, -ur, ur, vr, vr, 0.0, 1.0, hextra, dv0, &zbc);
									
					F[0][0][k] = dv0[0];
					F[1][0][k] = dv0[1];
					F[2][0][k] = dv0[2];

					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
								
			return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="UP"></param>
	/// <param name="G"></param>
	/// <param name="amax"></param>
	/// <param name="zc"></param>
	/// <param name="hextra"></param>
	/// <returns></returns>
	double boundaryNorth(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		
			for (int j = 0; j < n; j++)
			{
				
					// Enforce wall boundary on top of box (north)
				
						hl = UP[0][j][n-1]- zc[j][n-1];
						if (hl < 0.0) {
							hl = 0.0;
						}

						ul = (UP[1][j][n-1]) / (hl + hextra);
						vl = (UP[2][j][n-1]) / (hl + hextra);
						s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
						G[0][j][n-1] = dv0[0];
						G[1][j][n-1] = dv0[1];
						G[2][j][n-1] = dv0[2];
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
			}
	
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="UP"></param>
	/// <param name="F"></param>
	/// <param name="amax"></param>
	/// <param name="zc"></param>
	/// <param name="hextra"></param>
	/// <returns></returns>
	double boundaryEast(int n, double*** UP, double*** F, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

	
				for (int k = 0; k < n; k++)
				{
					hl = UP[0][n-1][k] - zc[n-1][k];
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
						ul = (UP[1][n-1][k]) / (hl + hextra);
						vl = (UP[2][n-1][k]) / (hl + hextra);
					}
					
					s.fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
					
					F[0][n-1][k] = dv0[0];
					F[1][n-1][k] = dv0[1];
					F[2][n-1][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
			
		}
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="UP"></param>
	/// <param name="G"></param>
	/// <param name="amax"></param>
	/// <param name="zc"></param>
	/// <param name="hextra"></param>
	/// <returns></returns>
	double boundarySouthY(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		double hr = 0.0;
		solver s;

		clearArray(dv0);
			for (int j = 1; j < n -1; j++)
			{
						
					hl = UP[0][j][1] - zc[j][1];
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
						ul = (UP[1][j][0]) / (hr + hextra);
						vl = (UP[2][j][0]) / (hr + hextra);
					}
					
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0, hextra, dv0, &zbc);
					
					G[0][j][0] = dv0[0];
					G[1][j][0] = dv0[1];
					G[2][j][0] = dv0[2];
						
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
			
			}		
		return amax;
 	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="UP"></param>
	/// <param name="G"></param>
	/// <param name="amax"></param>
	/// <param name="zc"></param>
	/// <param name="hextra"></param>
	/// <returns></returns>
	double secondColumnY(int n, double*** UP, double*** G, double amax, double** zc, double hextra)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;
		clearArray(dv0);

		for (int j = 1; j < n -1; j++)
			{			
					hl = UP[0][j][1] - zc[j][1];
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
						ul = (UP[1][j][0]) / (hl + hextra);
						vl = (UP[2][j][0]) / (hl + hextra);
					}
					
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0, hextra, dv0, &zbc);
					G[0][j][1] = dv0[0];
					G[1][j][1] = dv0[1];
					G[2][j][1] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				
			}
		
		return amax;

	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="G"></param>
	/// <param name="dwsex"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <param name="duy"></param>
	/// <param name="dwsey"></param>
	/// <param name="dvy"></param>
	/// <returns></returns>
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

		clearArray(dv0);
			for (int k = 2; k < n-1; k++)
			{
				for (int j = 0; j < n ; j++)
				{
					zbc = minmax(zc[j][k - 1], zc[j][k]).second;
				
				if ((UP[0][j][k - 1] - zbc) > 0.0) {
					hl = (UP[0][j][k-1] - zbc);
				}
				else if ((UP[0][j][k-1] - zbc) <= 0.0)
					{
						hl = 0.0;
					}			

				if ((UP[0][j][k] - zbc) > 0.0) {
					hr = UP[0][j][k] - zbc;
				}
				else if ((UP[0][j][k] - zbc) <= 0.0) {
						hr = 0.0;
					}
				

				if (hl > 0.0) {
					hl += 0.5 * dwsey[j][k-1];
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
					ul = UP[1][j][k - 1] + 0.5 * (duy[j][k-1] /
						(hl + hextra));
					vl = UP[2][j][k -1] + 0.5 * (dvy[j][k - 1] /
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
				/*if (dv0[0] > 0 && k == 40)
				{
					dv0[0] = 0;
				}*/
				G[0][j][k] = dv0[0];
				G[1][j][k] = dv0[1];
				G[2][j][k] = dv0[2];
				
				
							
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
			
		}
		}
		
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="G"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	double twentythreerowDownStream(double** zc, double*** UP, double amax, double*** G,double hextra,int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;
				
			for (int k = 1; k < 23; k++)
			{
				for (int j = 20; j < 22; j++)
				{					
								hl = UP[0][j][k - 1] - zc[j][k-1];
								if (hl < 0.0) {
									hl = 0.0;
								}
								ul = UP[1][j][k-1] / (hl + hextra);
								vl = UP[2][j][k-1] / (hl + hextra);
								s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
								
								G[0][j][k] = dv0[0];
								G[1][j][k] = dv0[1];
								G[2][j][k] = dv0[2];

								if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
								}
								else {
									amax = zbc;
								}
							}
						
				}
			
			return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="G"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	double twentyfourDam(double** zc, double*** UP, double amax, double*** G,  double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		
		int k = 23;
				for (int j = 20; j < 22; j++)
				{
					zbc = minmax(zc[j][k], zc[j][k]).second;
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

					G[0][j][k] = dv0[0];
					G[1][j][k] = dv0[1];
					G[2][j][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}			
		
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="G"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	double thirtyFiveLeftDamn(double** zc, double*** UP, double amax, double*** G, double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;
		int k = 34;
			for (int j = 20; j < 22; j++)
				{
				zbc = minmax(zc[j][k - 1], zc[j][k - 1]).second;
					
					if ((UP[0][j][k-1] - zbc) > 0)
					{
						hl = UP[0][j][k-1] - zbc;
						if (hl < 0)
						{
							hl = 0;
						}
						ul = (UP[1][j][k-1]) / (hl + hextra);
						vl = (UP[2][j][k-1]) / (hl + hextra);
					}
					else if ((UP[0][j][k-1] - zbc) <= 0)
					{
						hl = 0;
						if (hl == 0)
						{
							ul = 0; vl = 0;
						}
					}					
					
					s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);

					G[0][j][k] = dv0[0];
					G[1][j][k] = dv0[1];
					G[2][j][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
				}
			
		
		return amax;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="zc"></param>
	/// <param name="UP"></param>
	/// <param name="amax"></param>
	/// <param name="G"></param>
	/// <param name="hextra"></param>
	/// <param name="n"></param>
	/// <returns></returns>
	double leftDamn(double** zc, double*** UP, double amax, double*** G, double hextra, int n)
	{
		double zbc = 0.0;
		double hl = 0.0;
		double ul = 0.0;
		double vl = 0.0;
		double dv0[3];
		solver s;

		
			for (int k = 35; k < n -1; k++)
			{
				for (int j = 20; j < 22; j++)
				{
		
					hl = UP[0][j][k-1] - zc[j][k-1];
					
				if (hl < 0.0) {
					hl = 0.0;
				}

				ul = UP[1][j-1][k] / (hl + hextra);
				vl = UP[2][j-1][k] / (hl + hextra);
				s.fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
				
				
				G[0][j][k] = dv0[0];
				G[1][j][k] = dv0[1];
				G[2][j][k] = dv0[2];
					
					if ((zbc < amax) || (isnan(zbc) && (!isnan(amax)))) {
					}
					else {
						amax = zbc;
					}
	}
	}
	
		return amax;
	}
	/// <summary>
/// 
/// </summary>
/// <param name="arr"></param>
/// <param name="n"></param>
/// <param name="name"></param>
	void print3dArray(double*** arr, int n, string name)
	{
		cout << "*********************************************************************************************" << endl;
		cout << " Arr 0 dim" << endl;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				cout << arr[0][j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << " Arr 1 dim" << endl;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				cout << arr[1][j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << " Arr 2 dim" << endl;
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				//cout << j << k << endl;
				cout << arr[2][j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << "*********************************************************************************************" << endl;
		cout << "*********************************************************************************************" << endl;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="dv0"></param>
	void clearArray(double dv0[3])
	{
		dv0[0] = 0.0;
		dv0[1] = 0.0;
		dv0[2] = 0.0;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
	/// <param name="fileName"></param>
	void write3dOutputFile(double*** arr, int n, string fileName)
	{
		ofstream outStream;
		outStream.open("Output/" + fileName, ios::out);
		if (outStream.is_open())
		{

		//	outStream << "*********************************************************************************************" << endl;
	//		outStream << " Arr 0 dim" << endl;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					outStream << arr[0][j][k] << "\t";
				}
				outStream << endl;
			}
		  outStream << endl;
			outStream << " Arr 1 dim" << endl;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					outStream << arr[1][j][k] << "\t";
				}
				outStream << endl;
			}
			outStream << endl;
			outStream << " Arr 2 dim" << endl;
			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					outStream << arr[2][j][k] << "\t";
				}
				outStream << endl;
			}
			outStream << endl;

			outStream << "*********************************************************************************************" << endl;
			outStream.close();
		}
		else
		{
			cout << "Unable to open file" << endl;
		}

	}
};