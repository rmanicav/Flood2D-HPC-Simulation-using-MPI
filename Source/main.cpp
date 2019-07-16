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
#include "MPIXMatrix.h"
#include"omp.h"
#include<mpi.h>

using namespace std;
/// <summary>
/// 
/// </summary>
/// <returns></returns>
int main()
{

	helper help;
	floodData fd;
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
	int nf = n;// number of edges(will be used to compute fluxes at the interface)
	double	cr = fd.cr;
	double xc = cellsize / 2; //cellsize:L;// x coordinate of cell center
	double	yc = cellsize / 2; //:cellsize:L;// y coordinate of cell center
	double	nt = fd.nt;// number of time steps
	double	ntplot = fd.ntPlot;// plotting interval
	double	dt = fd.dt;// time steps
	double	dt2 = dt / cellsize;// ratio of dt / dx or dt / dy
	double amax = 0.0;
	double** zc = help.allocateMemory(n);
	help.clearArray(zc, n);
	double** h = help.allocateMemory(n);
	help.clearArray(h, n);

	double*** U = help.allocate3dMemory(n, n, n);
	double*** F = help.allocate3dMemory(n, n, n);
	double*** G = help.allocate3dMemory(n, n, n);
	double** u = help.allocateMemory(n);
	help.clearArray(u, n);
	double** v = help.allocateMemory(n);
       	double** wse = help.allocateMemory(n);
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
	double** sox = help.allocateMemory(n);
	help.clearArray(sox, n);
	double** soy = help.allocateMemory(n);
	help.clearArray(soy, n);
	double** sfx = help.allocateMemory(n);
	help.clearArray(sfx, n);
	double** sfy = help.allocateMemory(n);
	help.clearArray(sfy, n);
	double** wsep = help.allocateMemory(n);
	help.clearArray(wsep, n);
	double** up = help.allocateMemory(n);
	help.clearArray(up, n);
	double** vp = help.allocateMemory(n);
	help.clearArray(vp, n);
	double** hp = help.allocateMemory(n);
	help.clearArray(hp, n);

       	slope s;
	predictor p;
	fluxes f;
	corrector c;
      	limiter l;

   /*************MPI*********************/        
        int rank, size;
	MPI_Comm cartcomm = MPI_COMM_WORLD;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(cartcomm, &rank);

        MPIXMatrix mp;
        partition_data_t pd(size,n,n);




//if master initiate zc,wse,h,u v values
     if(rank ==0)
     {

	zc = fd.zc;
	cout << "n:" << n << endl;
	cout << "Completed ZC initialization" << endl;

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

                cout<<"Error occured "<<endl;

	//row 22 columns(24-34)- assign 6
	for (int j = 23; j < 34; j++) {
		wse[21][j] = fd.initWSE;//6 downstream
	}

	help.writeOutputFile(wse, n, "wse.txt");

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

      }



        /****************************ZC********************************/
        double *subzc= (double*)malloc(sizeof(double)*n*n);
        double *subwse= (double*)malloc(sizeof(double)*n*n);
        double *subh= (double*)malloc(sizeof(double)*n*n);
        double *subu= (double*)malloc(sizeof(double)*n*n);
        double *subv= (double*)malloc(sizeof(double)*n*n);


        for(int i =0;i<n;i++)
        {
          for(int j=0;j<n;j++)
          {
             subzc[i*n+j] = zc[i][j];
             subwse[i*n+j] = wse[i][j];
             subh[i*n+j] = h[i][j];
             subu[i*n+j] = u[i][j];
             subv[i*n+j] = v[i][j];             
          }
        }
        subzc = mp.scatter_exchange(subzc,pd,rank,cartcomm,n);
        subwse = mp.scatter_exchange(subwse,pd,rank,cartcomm,n);
        subh=mp.scatter_exchange(subh,pd,rank,cartcomm,n);
        subu=mp.scatter_exchange(subu,pd,rank,cartcomm,n);
        subv=mp.scatter_exchange(subv,pd,rank,cartcomm,n);
        
	clock_t start, end;
	try
	{
		int count = 0;
		double t = 0.0;
		double time_counter = 0;
                int startIter =0;

              /************************hot start - get last successfull number of iterations***************/
                int numOfFiles = help.getFileCount("Output",zc,wse,h,u,v,n);
                cout<<"Total Number of files in Output:"<<numOfFiles<<endl;
		if(numOfFiles < nt)
                {
                 startIter = numOfFiles;
                } 

             /************************start of simulation ************************************************/
                cout<<"Start iteration:"<<startIter<<endl;
                cout<<"Total iteration:"<<nt<<endl;
                //begin simulation
		for (int j = numOfFiles; j < nt; j++) {
			//start clock
			start = clock();
			simtime = j;
			cout << endl << "Iteration number :" << j << endl;
			/****limiter******************************************/
			l.flimiter(n, zc, dzcx, dzcy);
			cout << endl << "Completed Limiter 1 Function" << endl;
			l.flimiter(n, wse, dwsex, dwsey);
			cout << endl << "Completed Limiter 2 Function" << endl;
			l.flimiter(n, u, dux, duy);
			cout << endl << "Completed Limiter 3 Function" << endl;
			l.flimiter(n, v, dvx, dvy);
			cout << endl << "Completed Limiter 4 Function" << endl;
                         
                        //MPT_Barrier(MPI_COMM_WORLD);

			/************Slope calculation**********************************************/
			s.fslope(h, u, v, ManN, hextra, dzcx, dzcy, cellsize, n, sox, soy, sfx, sfy);
			cout << endl << "Completed Slope Function" << endl;

			/***********predictor step (estimate the values at half timestep)***********************************************/
			p.fpredictor(n, fd.gravity, nf, wse, h, u, v, dwsex, dwsey, dux, duy, dvx, dvy, dt2, dzcx, dzcy, epsilon, zc, sox, sfx, dt, soy, sfy, wsep, up, vp);
			cout << endl << "Completed Predictor Function" << endl;

			double*** UP = help.allocate3dMemory(n, n, n);
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

			cout << endl << "Completed Fluxes Function" << endl;
			//  Estimate the flux vectors on the next time step
			double*** uNew = help.allocate3dMemory(n, n, n);

			uNew = c.fcorrector(U, F, G, n, dt2, dt, sox, sfx, soy, sfy, grav);
			cout << endl << "Completed Corrector Function" << endl;
                        #pragma omp parallel num_threads(4)
                        {
                        #pragma omp sections
                        {
                        #pragma omp section
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
                        #pragma omp section
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
                        #pragma omp section
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
                        #pragma omp section
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
                       }
                       }

			cout << "Correct and re-assign values completed" << endl;
                  
                 if(rank ==0)
                 {
	           	//time steps
			t = t + dt;
			cr = (amax * dt) / cellsize;
			double ctrs;
			ctrs = simtime * dt;
			cout << "Simtime: " << simtime << endl;
			cout << "Ctrs: " << ctrs << endl;
			cout << "dts: " << dt << endl;
			cout << "epsilon: " << epsilon << endl;
			//help.freeMemory3d(uNew, n);
		  //if (fmod(ctrs,1) == 0)
		  //{

				
			help.writeOutputFile(h, n, "hOut_" + to_string(count) + ".txt");
			help.writeOutputFile(u, n, "uOut_" + to_string(count) + ".txt");
			help.writeOutputFile(v, n, "vOut_" + to_string(count) + ".txt");
			help.writeOutputFile(dzcx, n, "dzcx.txt");
			help.writeOutputFile(dzcy, n, "dzcy.txt");
			help.writeOutputFile(dwsex, n, "dwsex.txt");
			help.writeOutputFile(dwsey, n, "dwsey.txt");
			count++;
			//}
                }
		}
                if(rank ==0)
                {
		help.freeMemory(h, n);
		help.freeMemory(dzcx, n);
		help.freeMemory(dzcy, n);
		help.freeMemory(dwsex, n);
		help.freeMemory(dwsey, n);
		help.freeMemory(dux, n);
		help.freeMemory(duy, n);
		help.freeMemory(dvx, n);
		help.freeMemory(dvy, n);
		help.freeMemory(u, n);
		help.freeMemory(v, n);
		help.freeMemory(sox, n);
		help.freeMemory(soy, n);
		help.freeMemory(sfx, n);
		help.freeMemory(sfy, n);
		help.freeMemory(wsep, n);
		help.freeMemory(up, n);
		help.freeMemory(vp, n);
		help.freeMemory3d(U, n);
		help.freeMemory3d(F, n);
		help.freeMemory3d(G, n);
		help.freeMemory(hp, n);
                free(subzc);
                free(subh);
                free(subu);
                free(subv);

                }
	}
	catch (exception ex)
	{
		cout << "Exception occured -->" << ex.what() << endl;

	}
	cout << endl;
	cout << "Flood 2d Completed" << endl;
	end = clock();

	double cpu_time_used = ((double)end - start) / CLOCKS_PER_SEC;
	cout << "Time taken in seconds : " << cpu_time_used<<endl;

	return 0;
}

