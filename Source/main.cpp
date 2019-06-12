#include<iostream>
#include<fstream>
#include<ostream>
#include<algorithm>

#include "corrector.cpp"
#include "limiter.cpp"
#include "predictor.cpp"
#include "slope.cpp"
#include "fluxes.cpp"
using namespace std;
struct floodData
{
	int iv0[31],iv1[2];
	double gravity,manN;
	int cellSize,initV;
	double hextra, epsilon, L, nt, ntPlot, dt, initWSE, hWL;
	double cr;

};
void clearArray(double** arr, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			arr[i][j] = 0.0;
		}
	}
}
double ** allocateMemory(int n)
{
	double** arr;
	arr = (double**)malloc(sizeof(double*) * (int)n);
	for (int i = 0; i < n; i++)
	{

		arr[i] = (double*)malloc(sizeof(double) * (int)n);
	}
	return arr;
}
void readFromFile(floodData &fd)
{
		/*const signed char iv0[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 34, 35, 36, 37, 38, 39, 40, 41 };*
	  static int iv1[2] = { 10, 10 };*/
	ifstream infile;
	int count;
	infile.open("iv.txt");
	infile >> fd.gravity;
	cout << "----------------------------------\n";
	cout << "Gravity is " << fd.gravity << endl;
	infile >> fd.manN;
	cout << "Manning is " << fd.manN << endl;
	infile >> fd.hextra;
	cout << "hextra is " << fd.hextra << endl;
	infile >> fd.epsilon;
	cout << "epsilon is " << fd.epsilon << endl;
	infile >> fd.cellSize;
	cout << "Cell size is " << fd.cellSize << endl;
	infile >> fd.L;
	cout << "Box size is " << fd.L << endl;
	infile >> fd.nt;
	cout << "Number of time steps is " << fd.nt << endl;
	infile >> fd.ntPlot;
	cout << "Plotting interval is " << fd.ntPlot<< endl;
	infile >> fd.dt;
	cout << "Time steps is " << fd.dt << endl;
	cout << "----------------------------------\n\n";
	
	
	infile >> count;
	cout << "----iv0------\n";
	for (int i = 0; i < count; i++)
	{
		infile >> fd.iv0[i];
		cout << fd.iv0[i] << "\t";
	}
	cout << "\n-------------------------------\n\n";
	cout << endl;
	infile >> count;
	cout << "----iv1------\n";
	for (int i = 0; i < count; i++)
	{
		infile >> fd.iv1[i];
		cout << fd.iv1[i] << "\t";
	}
	cout << "\n-------------------------------\n\n";
	cout << "\n-------------------------------\n\n";
	cout << endl;
	infile >> fd.initV;
	cout << "Init value is " << fd.initV << endl;
	infile >> fd.initWSE;
	cout << "Initial value of WSE is " << fd.initWSE << endl;
	infile >> fd.hWL;
	cout << "Higher water level is  " << fd.hWL << endl;
	infile >> fd.cr;
	cout << "cr value is  " << fd.cr << endl;
	cout << "\n-------------------------------\n\n";

}
void simulate(int simtime,double dt,int ntplot)
{
	fstream fileid, fileid1, fileid2;
	const char* path1 = "C:/Users/raj/source/repos/Flood2d/Output/hsensor1.txt";
	const char* path2 = "C:/Users/raj/source/repos/Flood2d/Output/hsensor2.txt";
	const char* path3 = "C:/Users/raj/source/repos/Flood2d/Output/hsensor3.txt";
	//  Bed slope along X and Y
	//  writes the outputs for the sensors
	//  hnorm=zeros(size(num));     hsens_1=zeros(m,3); hsens_2=zeros(m,3); hsens_3=zeros(m,3);;

	//t=t+dt;
	 // dt= cr*cellsize/amax;
	 // cr=amax*dt/cellsize;
   //  Writes the output at each counters seconds
   
	  double ctrs = simtime * dt;
	  ofstream hfile;
   if (fmod(ctrs,ntplot) == 0.0) {
	   hfile.open("C:/Users/raj/source/repos/Flood2d/Output/hOut_",ios::out);
	   //    fprintf(fid1,'%.2f\r\n',h(:,:));
	   //   dlmwrite(outfile,x,'delimiter','\t','precision',12)
	   hfile.close();

	   //  % %   outputs the 2D-water depth plot
	   //  ctrs=sprintf('L%s',ctrs);
	   // saveas(gca,output);
   }
   double ctrs2 = simtime * dt;
   if (fmod(ctrs2,0.5) == 0) {
	   fileid.open(path1, ios::out);
	   fileid.close();

	   // dlmwrite(outfile2,hsens_1,'-append','newline','pc','delimiter','\t','precision',12)
	   // ctrs2=num2str(ctrs2);
   }
   double ctrs3= simtime * dt;
   if (fmod(ctrs3, 0.5) == 0) {
	   //fileid2.open("C:\Users\raj\source\repos\Flood2dOutput\hsensor2.txt","w");
	   fileid1.open(path2, ios::out);
	   // dlmwrite(outfile3,hsens_2,'-append','newline','pc','delimiter','\t','precision',12)
	   fileid1.close();

	   //  ctrs3=num2str(ctrs3);
   }

   double ctrs4 = simtime * dt;
   if (fmod(ctrs4, 0.5)) {
	   fileid2.open(path3, ios::out);
	   // dlmwrite(outfile3,hsens_2,'-append','newline','pc','delimiter','\t','precision',12)
	   fileid2.close();

	   //  ctrs3=num2str(ctrs3);
   }

}
int main(void)
{
	
	floodData fd;
	readFromFile(fd);
	int iv0[31];
	;
	
	
	int simtime;
	
	
	int iv1[2];
	
	//int iv1[2] = { 10, 10 };
	
	double UP[1764];
	
	
	static double unusedExpr[5292];
	int n;
	
	copy(begin(fd.iv1), end(fd.iv1), begin(iv1));
	copy(begin(fd.iv0), end(fd.iv0), begin(iv0));
	//double grav = 9.806;
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
	double** zc =	allocateMemory(n);
	clearArray(zc,n);
	
	for (int j = 0; j < n; j++) {
		for (int k = 0; k < n; k++) {
			zc[1][k] = fd.initV;
			zc[n][k] = fd.initV;
			zc[j][1] = fd.initV;
			zc[j][n] = fd.initV;
		}
	}

	double** wse = allocateMemory(n);
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
	*/
	limiter l;
	slope s;
	predictor p;
	fluxes f;
	corrector c;
	double** dzcx =allocateMemory(n);
	clearArray(dzcx, n);
	double** dzcy = allocateMemory(n);
	clearArray(dzcy, n);
	double** dwsex = allocateMemory(n);
	clearArray(dwsex, n);
	double** dwsey = allocateMemory(n);
	clearArray(dwsey, n);
	double** dux = allocateMemory(n);
	clearArray(dux, n);
	double** duy = allocateMemory(n);
	clearArray(duy, n);
	double** dvx = allocateMemory(n);
	clearArray(dvx, n);
	double** dvy = allocateMemory(n);
	clearArray(dvy, n);
	double** u  =allocateMemory(n);
	clearArray(u, n);
	double** v =allocateMemory(n);
	clearArray(v, n);
	double **sox = allocateMemory(n);
	clearArray(sox, n);
	double **soy = allocateMemory(n);
	clearArray(soy, n);
	double **sfx = allocateMemory(n);
	clearArray(sfx, n);
	double **sfy = allocateMemory(n);
	clearArray(sfy, n);
	double** h = allocateMemory(n);
	clearArray(h, n);
	for (int j = 0; j < 1; j++) {
		simtime = 1 + j;

		l.flimiter(n,zc, dzcx, dzcy);
	    l.flimiter(n,wse, dwsex, dwsey);
	    l.flimiter(n,u, dux, duy);
	    l.flimiter(n,v, dvx, dvy);

		//  Slope calculation
	  s.fslope(h, u, v, ManN, hextra, dzcx, dzcy, cellsize, n, sox, soy, sfx, sfy);

	  double **wsep =allocateMemory(n);
	  clearArray(wsep,n);
	  double **up = allocateMemory(n);
	  clearArray(up, n);
	  double **vp = allocateMemory(n);
	  clearArray(vp, n);

		//  predictor step (estimate the values at half timestep)
	   p.fpredictor(n,fd.gravity,nf,wse,h,u,v,dwsex,dwsey,dux,duy,dvx,dvy,dt2,dzcx,dzcy,epsilon,zc,sox,sfx,dt,soy,sfy,wsep,up,vp);


		//    Compute fluxes at the interfaces
		//f.ffluxes(UP, n, dwsex, dwsey, dux, duy, dvx, dvy, hextra, zc, F, G, amax);

	   double **U =allocateMemory(n); 
	   double **F = allocateMemory(n); 
	   double **G =allocateMemory(n);
		//  Estimate the flux vectors on the next time step
		c.fcorrector(U, F, G, n, dt2, dt, sox, sfx, soy, sfy, grav);


		//  U(:,:,:)=Unew(:,:,:);
		//  h(:,:)= Unew(:,:,1);                    % computed water depth (water level)
		//  if Unew(:,:,1)< epsilon
		//std::memset(&u[0], 0, 1764U * sizeof(double));
	//	std::memset(&v[0], 0, 1764U * sizeof(double));

		// else
		// u(:,:)= Unew(:,:,2)./(Unew(:,:,1)+hextra);
		// v(:,:)= Unew(:,:,3)./(Unew(:,:,1)+hextra);
	}

	// wse = Unew(:,:,1)+zc;                  % compute the new free surface height
	// if wse < 0
	//    wse = 0;
	// end*/

	//simulate(simtime, dt, ntplot);

	
	return 0;
}