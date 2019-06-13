#include<iostream>
#include<fstream>
#include<ostream>
#include<algorithm>
#include "helper.h"
#include "corrector.cpp"
#include "limiter.cpp"
#include "predictor.cpp"
#include "slope.cpp"
#include "fluxes.cpp"
using namespace std;


int main(void)
{
	helper help;
	floodData fd;
	help.readFromFile(&fd);
	int iv0[31];
	int simtime;
	int iv1[2];
	static double unusedExpr[5292];
	int n;
	copy(begin(fd.iv1), end(fd.iv1), begin(iv1));
	copy(begin(fd.iv0), end(fd.iv0), begin(iv0));
	double grav = fd.gravity;
	double ManN = fd.manN;
	double hextra = fd.hextra;
	double epsilon = fd.epsilon;// minimum depth of water(threshold for dry bed)(1cm)
	double	cellsize = fd.cellSize;// Cell size(dx = dy)m
	int	L = fd.L;// box size(200 * 200m)
	n = L / cellsize; // (n = 40) number of cells(nx = ny)
	int	m = n;
	int nf = n;// number of edges(will be used to compute fluxes at the interface)
	double	cr = fd.cr;

	// x = 0:cellsize:L;// array for edge coordinates in x
	// y = 0:cellsize:L;// array for edge coordinates in y

	double xc = cellsize / 2; //cellsize:L;// x coordinate of cell center
	double	yc = cellsize / 2; //:cellsize:L;// y coordinate of cell center

	double	nt = fd.nt;// number of time steps
	double	ntplot = fd.ntPlot;// plotting interval
	double	dt = fd.dt;// time steps
	double	dt2 = dt / cellsize;// ratio of dt / dx or dt / dy
	double hp = 0.0;
	double amax = 0.0;
	double** zc = help.allocateMemory(n);
	help.clearArray(zc, n);

	for (int j = 0; j < n; j++) {
		for (int k = 0; k < n; k++) {
			zc[1][k] = fd.initV;
			zc[n][k] = fd.initV;
			zc[j][1] = fd.initV;
			zc[j][n] = fd.initV;
		}
	}

	double** wse = help.allocateMemory(n);
	// Set up initial conditions
	//  wse=1.01*ones(n,n);                           % Initial water surface elevation
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			wse[i][j] = fd.initWSE;
		}

	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			wse[i][j] = fd.hWL;
		}
	}
	/*

	//  0.5*n = 20 or higher water level
	for (i0 = 0; i0 < 31; i0++) {
		for (j = 0; j < 2; j++) {
			wse[(j + 42 * iv0[i0]) + 20] = 99999.0;
		}
	}*/
	double** h = help.allocateMemory(n);
	double** dzcx = help.allocateMemory(n);
	help.clearArray(dzcx, n);
	double** dzcy = help.allocateMemory(n);
	help.clearArray(dzcy, n);
	double** dwsex = help.allocateMemory(n);
	help.clearArray(dwsex, n);
	double** dwsey = help.allocateMemory(n);
	help.clearArray(dwsey, n);
	double** dux = help.allocateMemory(n);
	help.clearArray(dux, n);
	double** duy = help.allocateMemory(n);
	help.clearArray(duy, n);
	double** dvx = help.allocateMemory(n);
	help.clearArray(dvx, n);
	double** dvy = help.allocateMemory(n);
	help.clearArray(dvy, n);
	double** u = help.allocateMemory(n);
	help.clearArray(u, n);
	double** v = help.allocateMemory(n);
	help.clearArray(v, n);
	double** sox = help.allocateMemory(n);
	help.clearArray(sox, n);
	double** soy = help.allocateMemory(n);
	help.clearArray(soy, n);
	double** sfx = help.allocateMemory(n);
	help.clearArray(sfx, n);
	double** sfy = help.allocateMemory(n);
	help.clearArray(sfy, n);
	help.clearArray(h, n);
	double** wsep = help.allocateMemory(n);
	help.clearArray(wsep, n);
	double** up = help.allocateMemory(n);
	help.clearArray(up, n);
	double** vp = help.allocateMemory(n);
	help.clearArray(vp, n);
	double*** U = help.allocate3dMemory(n, fd.dim);
	double*** F = help.allocate3dMemory(n, fd.dim);
	double*** G = help.allocate3dMemory(n, fd.dim);
	//  0.5*n = 20 or higher water level
	//  water depth to water surface elevatin relation ship
	for (int j = 0; j < 42; j++) {
		for (int k = 0; k < 42; k++) {
			h[j][k] = wse[j][k] - zc[j][k];
			h[j][k] = 0.0;
			h[j][k] = 0.0;
			h[j][k] = 0.0;
			h[j][k] = 0.0;
			if (h[j][k] < 0.0) {
				h[j][k] = 0.0;
			}
		}
	}





	//  fluxes in the $x$ direction
	//  fluxes in the $y$ direction
	//  friction slope Sf
	//  Bed slope So
	//  variables on the new time level
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			U[i][i][0] = h[i][j];
			u[i][j] = U[i][j][1] / (U[i][j][1] + 0.1);
			v[i][j] = U[i][j][2] / (U[i][j][2] + 0.1);
		}
	}

	limiter l;
	slope s;
	predictor p;
	fluxes f;
	corrector c;


	for (int j = 0; j < 1; j++) {
		simtime = 1 + j;

		l.flimiter(n, zc, dzcx, dzcy);
		l.flimiter(n, wse, dwsex, dwsey);
		l.flimiter(n, u, dux, duy);
		l.flimiter(n, v, dvx, dvy);

		//  Slope calculation
		s.fslope(h, u, v, ManN, hextra, dzcx, dzcy, cellsize, n, sox, soy, sfx, sfy);



		//  predictor step (estimate the values at half timestep)
		p.fpredictor(n, fd.gravity, nf, wse, h, u, v, dwsex, dwsey, dux, duy, dvx, dvy, dt2, dzcx, dzcy, epsilon, zc, sox, sfx, dt, soy, sfy, wsep, up, vp);

		double*** UP = help.allocate3dMemory(n, fd.dim);

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				UP[i][0][0] = wsep[i][j];
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				UP[0][j][0] = u[i][j];
			}
		}
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				UP[0][0][j] = v[i][j];
			}
		}


		//    Compute fluxes at the interfaces
		f.ffluxes(UP, n, dwsex, dwsey, dux, duy, dvx, dvy, hextra, zc, F, G, amax);



		//  Estimate the flux vectors on the next time step
		c.fcorrector(U, F, G, n, dt2, dt, sox, sfx, soy, sfy, grav);


		//simulate(simtime, dt, ntplot);

	}
		return 0;
	}
