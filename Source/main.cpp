// Main MPI-based flood simulation driver.
// Each MPI rank operates on a subdomain of the 2D grid.

#include<iostream>
#include<fstream>
#include<ostream>
#include<algorithm>
#include<time.h>
#include<math.h>
#include "helper.h"
#include "corrector.cpp"
#include "limiter.cpp"
#include "predictor.cpp"
#include "slope.cpp"
#include "fluxes.cpp"
#include "hypograph.h"
#include "dem.h"
using namespace std;
/// <summary>
/// 
/// </summary>
/// <returns></returns>
int main(void)
{

	helper help;
	floodData fd;
	dem d;
	hydrograph hg;
	hyg shyg;
	double xloc, yloc;
	/***********************************Taum Sauke************************************************/
	/*d.readDemFile("Input/Dem.txt");
	d.readLocData("Input/inflowlocations.txt", xloc, yloc);*/
	//hg.readHygData(&shyg,"Input/hydrographs.txt");
	/*int nrows = d.nrows;
	int ncols = d.ncols;
	double cellSize = d.cellSize;
	double xllCorner = d.xllCorner;
	double yllCorner = d.yllCorner;
	double noData = d.nodata;*/

	
	/*****************************************zc****************************************************/
	help.readFromFile(&fd);
	int simtime;
	int n;

	double grav = fd.gravity;
	double ManN = fd.manN;
	double hextra = fd.hextra;
	double epsilon = fd.epsilon;// minimum depth of water(threshold for dry bed)(1cm)
	int	cellsize = fd.cellSize;// Cell size(dx = dy)m
	int	L = fd.L;// box size(200 * 200m)
	n = L / cellsize; // (n = 40) number of cells(nx = ny)
	int	m = n;
	int nrows = fd.x;
	int ncols = fd.y;
	int ndepth = fd.z;
	int nf = n;// number of edges(will be used to compute fluxes at the interface)
	double	cr = fd.cr;
	double xc = cellsize / 2; //cellsize:L;// x coordinate of cell center
	double	yc = cellsize / 2; //:cellsize:L;// y coordinate of cell center
	double	nt = fd.nt;// number of time steps
	double	ntplot = fd.ntPlot;// plotting interval
	double	dt = fd.dt;// time steps
	double	dt2 = dt / cellsize;// ratio of dt / dx or dt / dy
	double amax = 0.0;
	double** zc = help.allocateMemoryTest(nrows,ncols);
	help.clearArray(zc, nrows,ncols);
	zc = fd.zc;
	std::cout << "n:" << n << endl;
	std::cout << "Completed ZC initialization" << endl;

	double** wse = help.allocateMemory(nrows,ncols);
	// Set up initial conditions
	//  wse=1.01*ones(n,n);                           % Initial water surface elevation


	//0-20 - assign 11
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < n; j++) {
			wse[i][j] = fd.hWL;// 11 upstream
		}
	}


	//21,22 - assign 9999
	for (int i = 20; i < 22; i++) {
		for (int j = 0; j < n; j++) {
			wse[i][j] = fd.initV;// 9999
		}
	}

	//23 -42 rows - assign 6
	for (int i = 22; i < n; i++) {
		for (int j = 0; j < n; j++) {
			wse[i][j] = fd.initWSE;//6 downstream
		}

	}
	//row 21 - columns(24-34) - assign 11
	for (int j = 23; j < 34; j++)
	{
		wse[20][j] = fd.hWL;// 11 upstream
	}

	//row 22 columns(24-34)- assign 6
	for (int j = 23; j < 34; j++) {
		wse[21][j] = fd.initWSE;//6 downstream
	}

	//	help.printArray(wse, n, "wse");
	help.writeOutputFile(wse, nrows,ncols, "wse.txt");
	double** h = help.allocateMemory(nrows,ncols);
	help.clearArray(h, nrows,ncols);

	//  0.5*n = 20 or higher water level
	//  water depth to water surface elevatin relation ship
	for (int j = 0; j < n; j++) {
		for (int k = 0; k < n; k++) {
			h[j][k] = wse[j][k] - zc[j][k];
			h[0][k] = 0.0;
			h[n - 1][k] = 0.0;
			h[j][0] = 0.0;
			h[j][n - 1] = 0.0;
			if (h[j][k] < 0.0) {
				h[j][k] = 0.0;
			}
		}
	}
	//help.writeOutputFile(h, n, "h.txt");
	//help.printArray(h, n, "h");
	double*** U = help.allocate3dMemory(nrows,ncols,ndepth);
	double*** F = help.allocate3dMemory(nrows,ncols,ndepth);
	double*** G = help.allocate3dMemory(nrows,ncols,ndepth);
	double** u = help.allocateMemory(nrows,ncols);
	help.clearArray(u, nrows, ncols);
	double** v = help.allocateMemory(nrows,ncols);
	help.clearArray(v, nrows, ncols);

	//  fluxes in the $x$ direction
	//  fluxes in the $y$ direction
	//  friction slope Sf
	//  Bed slope So
	//  variables on the new time level
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			U[0][i][j] = h[i][j];
			U[1][i][j] = u[i][j] * h[i][j];
			U[2][i][j] = v[i][j] * h[i][j];
		}
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {

			u[i][j] = U[1][i][j] / (U[0][i][j] + fd.hextra);
			v[i][j] = U[2][i][j] / (U[0][i][j] + fd.hextra);
		}
	}

	limiter l;
	slope s;
	predictor p;
	fluxes f;
	corrector c;
	double** dzcx = help.allocateMemory(nrows,ncols);
	help.clearArray(dzcx, nrows, ncols);
	double** dzcy = help.allocateMemory(nrows,ncols);
	help.clearArray(dzcy, nrows, ncols);
	double** dwsex = help.allocateMemory(nrows,ncols);
	help.clearArray(dwsex, nrows, ncols);
	double** dwsey = help.allocateMemory(nrows,ncols);
	help.clearArray(dwsey, nrows, ncols);
	double** dux = help.allocateMemory(nrows,ncols);
	help.clearArray(dux, nrows, ncols);
	double** duy = help.allocateMemory(nrows,ncols);
	help.clearArray(duy, nrows, ncols);
	double** dvx = help.allocateMemory(nrows,ncols);
	help.clearArray(dvx, nrows, ncols);
	double** dvy = help.allocateMemory(nrows,ncols);
	help.clearArray(dvy, nrows, ncols);
	double** sox = help.allocateMemory(nrows,ncols);
	help.clearArray(sox, nrows, ncols);
	double** soy = help.allocateMemory(nrows,ncols);
	help.clearArray(soy, nrows, ncols);
	double** sfx = help.allocateMemory(nrows,ncols);
	help.clearArray(sfx, nrows, ncols);
	double** sfy = help.allocateMemory(nrows,ncols);
	help.clearArray(sfy, nrows, ncols);
	double** wsep = help.allocateMemory(nrows,ncols);
	help.clearArray(wsep, nrows, ncols);
	double** up = help.allocateMemory(nrows,ncols);
	help.clearArray(up, nrows, ncols);
	double** vp = help.allocateMemory(nrows,ncols);
	help.clearArray(vp, nrows, ncols);
	double** hp = help.allocateMemory(nrows,ncols);
	help.clearArray(hp, nrows, ncols);
	//% Bed slope along Xand Y
	for (int j = 0; j < n - 1; j++)
	{
		for (int k = 0; k < n - 1; k++)
		{
			dzcx[j][k] = zc[j + 1][k] - zc[j][k];
			dzcy[j][k] = zc[j][k + 1] - zc[j][k];
			dzcx[1][k] = 0;    dzcy[1][k] = 0;
			dzcx[n - 1][k] = 0;    dzcy[n - 1][k] = 0;
			dzcx[j][1] = 0;    dzcy[j][1] = 0;
			dzcx[j][n - 1] = 0;    dzcy[j][n - 1] = 0;
		}
	}
	clock_t start, end;
	try
	{
		int count = 0;
		double t = 0.0;
		double time_counter = 0;
		int startIter = 0;
		/************************hot start - get last successfull number of iterations***************/
		int numOfFiles = help.getFileCount("Output", zc, wse, h, u, v, nrows,ncols);
		cout << "Total Number of files in Output:" << numOfFiles << endl;
		if (numOfFiles <= nt)
		{
			startIter = numOfFiles;
		}
		for (int j = startIter; j < nt; j++) {
			//start clock
			start = clock();
			simtime = j;
			std::cout << endl << "Iteration number :" << j << endl;
			/****limiter******************************************/
			l.flimiter(n, zc, dzcx, dzcy);
			std::cout << endl << "Completed Limiter 1 Function" << endl;
			l.flimiter(n, wse, dwsex, dwsey);
			std::cout << endl << "Completed Limiter 2 Function" << endl;
			l.flimiter(n, u, dux, duy);
			std::cout << endl << "Completed Limiter 3 Function" << endl;
			l.flimiter(n, v, dvx, dvy);
			std::cout << endl << "Completed Limiter 4 Function" << endl;



			/************Slope calculation**********************************************/
			s.fslope(h, u, v, ManN, hextra, dzcx, dzcy, cellsize, n, sox, soy, sfx, sfy);
			std::cout << endl << "Completed Slope Function" << endl;

			/***********predictor step (estimate the values at half timestep)***********************************************/
			p.fpredictor(n, fd.gravity, nf, wse, h, u, v, dwsex, dwsey, dux, duy, dvx, dvy, dt2, dzcx, dzcy, epsilon, zc, sox, sfx, dt, soy, sfy, wsep, up, vp);
			std::cout << endl << "Completed Predictor Function" << endl;
			double*** UP = help.allocate3dMemory(nrows,ncols,ndepth);
			//assign 0 dim  with wsep value
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					UP[0][i][j] = wsep[i][j];
					help.checkForNan(UP[0][i][j]);
				}
			}
			//hp = wsep - zc
			//assign 1 dim with up value, 2 - vp
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					hp[i][j] = wsep[i][j] - zc[i][j];
					if (hp[i][j] < 0)
					{
						hp[i][j] = 0.0;
					}
					help.checkForNan(hp[i][j]);
				}
			}
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					help.fixvp(vp, i, j);
					help.fixup(up, i, j);

					UP[1][i][j] = up[i][j] * hp[i][j];
					UP[2][i][j] = vp[i][j] * hp[i][j];
					help.checkForNan(UP[1][i][j]);
					help.checkForNan(UP[2][i][j]);
				}
			}
			//    Compute fluxes at the interfaces
			f.ffluxes(UP, n, dwsex, dwsey, dux, duy, dvx, dvy, hextra, zc, F, G, amax);
			std::cout << endl << "Completed Fluxes Function" << endl;
			//  Estimate the flux vectors on the next time step
			double*** uNew = help.allocate3dMemory(nrows,ncols,ndepth);

			uNew = c.fcorrector(U, F, G, n, dt2, dt, sox, sfx, soy, sfy, grav);
			std::cout << endl << "Completed Corrector Function" << endl;

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					for (int k = 0; k < n; k++)
					{
						U[i][j][k] = uNew[i][j][k];
						help.checkForNan(U[i][j][k]);
					}
				}
			}

			//reassign values after correction
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{

					//check epsilon
					if (uNew[0][i][j] < epsilon)
					{
						u[i][j] = 0;
						v[i][j] = 0;
					}
					else
					{
						u[i][j] = uNew[1][i][j] / (uNew[0][i][j] + fd.hextra);
						v[i][j] = uNew[2][i][j] / (uNew[0][i][j] + fd.hextra);
						help.checkForNan(u[i][j]);
						help.checkForNan(v[i][j]);
					}
				}
			}

			for (int j = 0; j < n; j++)
			{
				for (int k = 0; k < n; k++)
				{
					//check epsilon
					if (uNew[0][j][k] < epsilon)
					{
						h[j][k] = 0;
					}
					else
					{
						// computed water depth(water level)
						h[j][k] = uNew[0][j][k];
						help.checkForNan(h[j][k]);
					}
				}
			}

			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					//compute the new free surface height
					wse[i][j] = h[i][j] + zc[i][j];
					if (wse[i][j] < 0)
					{
						wse[i][j] = 0;
					}
					help.checkForNan(wse[i][j]);
				}
			}
			std::cout << "Correct and re-assign values completed" << endl;

			//time steps
			t = t + dt;
			cr = (amax * dt) / cellsize;
			double ctrs;
			ctrs = simtime * dt;
			std::cout << "Simtime: " << simtime << endl;
			std::cout << "Ctrs: " << ctrs << endl;
			std::cout << "dts: " << dt << endl;
			std::cout << "epsilon: " << epsilon << endl;
			//help.freeMemory3d(uNew, n);
		  //if (fmod(ctrs,1) == 0)
		  //{

				
			help.writeOutputFile(h, nrows,ncols, "hOut_" + to_string(count) + ".txt");
			help.writeOutputFile(u, nrows,ncols, "uOut_" + to_string(count) + ".txt");
			help.writeOutputFile(v, nrows,ncols, "vOut_" + to_string(count) + ".txt");
			help.writeOutputFile(dzcx, nrows,ncols, "dzcx.txt");
			help.writeOutputFile(dzcy, nrows,ncols, "dzcy.txt");
			help.writeOutputFile(dwsex, nrows,ncols, "dwsex.txt");
			help.writeOutputFile(dwsey, nrows,ncols, "dwsey.txt");
			count++;
			//}
			  		
		}
		help.freeMemory(h, nrows);
		help.freeMemory(dzcx, nrows);
		help.freeMemory(dzcy, nrows);
		help.freeMemory(dwsex, nrows);
		help.freeMemory(dwsey, nrows);
		help.freeMemory(dux, nrows);
		help.freeMemory(duy, nrows);
		help.freeMemory(dvx, nrows);
		help.freeMemory(dvy, nrows);
		help.freeMemory(u, nrows);
		help.freeMemory(v, nrows);
		help.freeMemory(sox, nrows);
		help.freeMemory(soy, nrows);
		help.freeMemory(sfx, nrows);
		help.freeMemory(sfy, nrows);
		help.freeMemory(wsep, nrows);
		help.freeMemory(up, nrows);
		help.freeMemory(vp, nrows);
		help.freeMemory3d(U, nrows,ncols);
		help.freeMemory3d(F, nrows,ncols);
		help.freeMemory3d(G, nrows,ncols);
		help.freeMemory(hp, nrows);
	}
	catch (exception ex)
	{
		std::cout << "Exception occured -->" << ex.what() << endl;

	}
	std::cout << endl;
	std::cout << "Flood 2d Completed" << endl;
	end = clock();

	double cpu_time_used = ((double)end - start) / CLOCKS_PER_SEC;
	std::cout << "Time taken in seconds : " << cpu_time_used;

	return 0;
}

