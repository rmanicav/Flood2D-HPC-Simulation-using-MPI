#include<iostream>
#include<fstream>
#include<ostream>
#include "corrector.cpp"
#include "limiter.cpp"
#include "predictor.cpp"
#include "slope.cpp"
#include "fluxes.cpp"
using namespace std;

int main(void)
{
	int i0;
	int j;
	double zc[1764];
	int k;
	double wse[1764];
	const signed char iv0[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 34, 35, 36, 37, 38, 39, 40, 41 };

	double U[5292];
	double h[1764];
	double u[1764];
	double v[1764];
	int simtime;
	//emxArray_real_T *dzcx;
	double dzcx[1764];
	double dzcy[1764];
	double dwsex[1764];
	double dwsey[1764]; 
	double dux[1764];
	double duy[1764];
	double dvx[1764];
	double dvy[1764];
	double sox[1764];
	double soy[1764];
	double sfx[100];
	double sfy[100];
	static int iv1[2] = { 10, 10 };

	double wsep[1764];
	double up[1764];
	double vp[1764];
	//signed char fileid;
	static double unusedExpr[5292];
	int n = 42;
	//
	double grav = 9.806;
	double ManN = 0.001;
	
	double hextra = 0.1;
	double epsilon = 0.1;// minimum depth of water(threshold for dry bed)(1cm)
	int	cellsize = 5;// Cell size(dx = dy)m
	int	L = 210;// box size(200 * 200m)
	n = L / cellsize; // (n = 40) number of cells(nx = ny)
	int	m = n;
	int nf = n;// number of edges(will be used to compute fluxes at the interface)
	double	cr = 0.05;

	// x = 0:cellsize:L;// array for edge coordinates in x
		// y = 0:cellsize:L;// array for edge coordinates in y

	double xc = cellsize / 2; //cellsize:L;// x coordinate of cell center
	double	yc = cellsize / 2; //:cellsize:L;// y coordinate of cell center

	double	nt = 20000;// number of time steps
	double	ntplot = 2;// plotting interval
	double	dt = 0.05;// time steps
	double	dt2 = dt / cellsize;// ratio of dt / dx or dt / dy
	double hp = 0.0;
	double amax = 0.0;
	for (i0 = 0; i0 < 1764; i0++) {
		zc[i0] = 1.0;
	}

	for (j = 0; j < 42; j++) {
		for (k = 0; k < 42; k++) {
			zc[42 * k] = 99999.0;
			zc[41 + 42 * k] = 99999.0;
			zc[j] = 99999.0;
			zc[1722 + j] = 99999.0;
			for (i0 = 0; i0 < 23; i0++) {
				zc[20 + 42 * i0] = 99999.0;
			}

			for (i0 = 0; i0 < 8; i0++) {
				zc[20 + 42 * (34 + i0)] = 99999.0;
			}

			for (i0 = 0; i0 < 23; i0++) {
				zc[21 + 42 * i0] = 99999.0;
			}

			for (i0 = 0; i0 < 8; i0++) {
				zc[21 + 42 * (34 + i0)] = 99999.0;
			}
		}
	}

	//  Set up initial conditions
	//  wse=1.01*ones(n,n);                           % Initial water surface elevation  
	for (i0 = 0; i0 < 1764; i0++) {
		wse[i0] = 6.0;
	}

	for (i0 = 0; i0 < 42; i0++) {
		for (j = 0; j < 21; j++) {
			wse[j + 42 * i0] = 11.0;
		}
	}

	//  0.5*n = 20 or higher water level
	for (i0 = 0; i0 < 31; i0++) {
		for (j = 0; j < 2; j++) {
			wse[(j + 42 * iv0[i0]) + 20] = 99999.0;
		}
	}

	//  0.5*n = 20 or higher water level
	//  water depth to water surface elevatin relation ship
	for (j = 0; j < 42; j++) {
		for (k = 0; k < 42; k++) {
			h[j + 42 * k] = wse[j + 42 * k] - zc[j + 42 * k];
			h[42 * k] = 0.0;
			h[41 + 42 * k] = 0.0;
			h[j] = 0.0;
			h[1722 + j] = 0.0;
			if (h[j + 42 * k] < 0.0) {
				h[j + 42 * k] = 0.0;
			}
		}
	}

	// Intialize the arrays
	std::memset(&U[0], 0, 5292U * sizeof(double));

	//  fluxes in the $x$ direction
	//  fluxes in the $y$ direction
	//  friction slope Sf
	//  Bed slope So
	//  variables on the new time level
	for (i0 = 0; i0 < 42; i0++) {
		for (j = 0; j < 42; j++) {
			U[j + 42 * i0] = h[j + 42 * i0];
			U[1764 + (j + 42 * i0)] = 0.0 * h[j + 42 * i0];
			U[3528 + (j + 42 * i0)] = 0.0 * h[j + 42 * i0];
			u[j + 42 * i0] = U[1764 + (j + 42 * i0)] / (U[j + 42 * i0] + 0.1);
			v[j + 42 * i0] = U[3528 + (j + 42 * i0)] / (U[j + 42 * i0] + 0.1);
		}
	}

	fstream fileid,fileid1,fileid2;
	const char* path = "C:/Users/raj/source/repos/Flood2dOutput/hsensor1.txt";
	const char* path1 = "C:/Users/raj/source/repos/Flood2dOutput/hsensor2.txt";
	const char* path2 = "C:/Users/raj/source/repos/Flood2dOutput/hsensor3.txt";
	//  Bed slope along X and Y
	//  writes the outputs for the sensors
	//  hnorm=zeros(size(num));     hsens_1=zeros(m,3); hsens_2=zeros(m,3); hsens_3=zeros(m,3);; 
	fileid.open(path,ios::out);
	fileid1.open(path, ios::out);
	fileid2.open(path, ios::out);

        
	//  tcntr=0;
	for (j = 0; j < 20000; j++) {
		simtime = 1 + j;

		flimiter(n,zc, dzcx, dzcy);
		flimiter(n,wse, dwsex, dwsey);
		flimiter(n,u, dux, duy);
		flimiter(n,v, dvx, dvy);

		//  Slope calculation
		fslope(h, u, v, ManN, hextra, dzcx, dzcy, cellsize, n, sox, soy, sfx, sfy);
		

		//  predictor step (estimate the values at half timestep)
		fpredictor(nf,wse,h,u,v,dwsex,dwsey,dux,duy,dvx,dvy,dt2,dzcx,dzcy,epsilon,zc,sox,sfx,dt,soy,sfy,wsep,up,vp);
		

		//    Compute fluxes at the interfaces
		// ffluxes(UP,n,dwsex,dwsey,dux,duy,dvx,dvy,hextra,zc);
		//  Estimate the flux vectors on the next time step
		//fcorrector(U, sox, sfx, iv1, soy, sfy, iv1, grav, unusedExpr);

		//  U(:,:,:)=Unew(:,:,:);
		//  h(:,:)= Unew(:,:,1);                    % computed water depth (water level)      
		//  if Unew(:,:,1)< epsilon
		std::memset(&u[0], 0, 1764U * sizeof(double));
		std::memset(&v[0], 0, 1764U * sizeof(double));

		// else
		// u(:,:)= Unew(:,:,2)./(Unew(:,:,1)+hextra);
		// v(:,:)= Unew(:,:,3)./(Unew(:,:,1)+hextra);
	}

	// wse = Unew(:,:,1)+zc;                  % compute the new free surface height 
	// if wse < 0
	//    wse = 0;
	// end
	   //t=t+dt;
	   dt= cr*cellsize/amax;
	   cr=amax*dt/cellsize;
	//  Writes the output at each counters seconds
	   double ctrs = simtime * dt;
	   ofstream hfile;
	if (fmod(ctrs,ntplot) == 0.0) {
		hfile.open("C:/Users/raj/source/repos/Flood2dOutput/hOut_",ios::out);
		//    fprintf(fid1,'%.2f\r\n',h(:,:));
		//   dlmwrite(outfile,x,'delimiter','\t','precision',12)
		hfile.close();

		//  % %   outputs the 2D-water depth plot
		//  ctrs=sprintf('L%s',ctrs);
		// saveas(gca,output);
	}
	double ctrs2 = simtime * dt;
	if (fmod(ctrs2,0.5) == 0) {
		//fileid.open("C:/Users/raj/source/repos/Flood2d/Output/hsensor1.txt","w");
		fileid.close();

		// dlmwrite(outfile2,hsens_1,'-append','newline','pc','delimiter','\t','precision',12) 
		// ctrs2=num2str(ctrs2);
	}
	double ctrs3= simtime * dt;
	if (fmod(ctrs3, 0.5) == 0) {
		//fileid2.open("C:\Users\raj\source\repos\Flood2dOutput\hsensor2.txt","w");

		// dlmwrite(outfile3,hsens_2,'-append','newline','pc','delimiter','\t','precision',12) 
		fileid1.close();

		//  ctrs3=num2str(ctrs3);
	}

	double ctrs4 = simtime * dt;
	if (fmod(ctrs4, 0.5)) {
		//fileid2.open("C:\Users\raj\source\repos\Flood2dOutput\hsensor2.txt","w");

		// dlmwrite(outfile3,hsens_2,'-append','newline','pc','delimiter','\t','precision',12) 
		fileid2.close();

		//  ctrs3=num2str(ctrs3);
	}

	return 0;
}