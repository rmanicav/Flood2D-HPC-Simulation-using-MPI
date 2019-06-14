#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<ostream>
#include<iomanip>
using namespace std;
struct floodData
{
	int iv0[31], iv1[2];
	double gravity, manN, cellSize, initV;
	int  L, dim;
	double hextra, epsilon, nt, ntPlot, dt, initWSE, hWL;
	double cr;
	double **zc;
};
class helper
{
public:
	
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
	double** allocateMemory(int n)
	{
		double** arr;
		arr = (double**)malloc(sizeof(double*) * (int)n);
		for (int i = 0; i < n; i++)
		{

			arr[i] = (double*)malloc(sizeof(double) * (int)n);
		}
		return arr;
	}
	double*** allocate3dMemory(int n, int dim)
	{
		double*** arr = new double** [n];
		
		for (int i = 0; i < n; i++)
		{
			arr[i] = new double* [n];
			for (int j = 0; j < n; j++)
			{
				arr[i][j] = new double[n];
			}			
		}
		return arr;
	}
	void readFromFile(floodData* fd)
	{
		/*const signed char iv0[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 34, 35, 36, 37, 38, 39, 40, 41 };*
	  static int iv1[2] = { 10, 10 };*/
		ifstream infile;
		int count;
		infile.open("iv.txt");
		cout << "----------------------------------\n";
		infile >> fd->gravity;
		cout << "Gravity:" << fd->gravity << endl;
		infile >> fd->manN;
		cout << "Manning:" << fd->manN << endl;
		infile >> fd->hextra;
		cout << "hextra:" << fd->hextra << endl;
		infile >> fd->epsilon;
		cout << "epsilon:" << fd->epsilon << endl;
		infile >> fd->cellSize;
		cout << "Cell size:" << fd->cellSize << endl;
		infile >> fd->L;
		cout << "Box size:" << fd->L << endl;
		infile >> fd->nt;
		cout << "Number of time steps:" << fd->nt << endl;
		infile >> fd->ntPlot;
		cout << "Plotting interval:" << fd->ntPlot << endl;
		infile >> fd->dt;
		cout << "Time steps:" << fd->dt << endl;
		//cout << "----------------------------------\n\n";
		infile >> count;
		cout << "----iv0------\n";
		for (int i = 0; i < count; i++)
		{
			infile >> fd->iv0[i];
			//		cout << fd->iv0[i] << "\t";
		}
		//cout << "\n-------------------------------\n\n";
		cout << endl;
		infile >> count;
		//cout << "----iv1------\n";
		for (int i = 0; i < count; i++)
		{
			infile >> fd->iv1[i];
			//		cout << fd->iv1[i] << "\t";
		}
		//cout << "\n-------------------------------\n\n";
		//cout << "\n-------------------------------\n\n";
		cout << endl;
		infile >> fd->initV;
		cout << "Init value:" << fd->initV << endl;
		infile >> fd->initWSE;
		cout << "WSE downstream value:" << fd->initWSE << endl;
		infile >> fd->hWL;
		cout << "WSE upstream value:" << fd->hWL << endl;
		infile >> fd->cr;
		cout << "cr value:" << fd->cr << endl;
		infile >> fd->dim;
		cout << "Dimension:" << fd->dim << endl;
		cout << "\n-------------------------------\n\n";
		fd->zc = allocateMemory(42);
		clearArray(fd->zc,42);
		for (int i = 0; i < 42; i++)
		{
			for (int j = 0; j < 42; j++)
			{
				infile >> fd->zc[i][j];
			}
		}
		printArray(fd->zc, 42, "zc");
	}
	void simulate(int simtime, double dt, int ntplot)
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
		if (fmod(ctrs, ntplot) == 0.0) {
			hfile.open("C:/Users/raj/source/repos/Flood2d/Output/hOut_", ios::out);
			//    fprintf(fid1,'%.2f\r\n',h(:,:));
			//   dlmwrite(outfile,x,'delimiter','\t','precision',12)
			hfile.close();

			//  % %   outputs the 2D-water depth plot
			//  ctrs=sprintf('L%s',ctrs);
			// saveas(gca,output);
		}
		double ctrs2 = simtime * dt;
		if (fmod(ctrs2, 0.5) == 0) {
			fileid.open(path1, ios::out);
			fileid.close();

			// dlmwrite(outfile2,hsens_1,'-append','newline','pc','delimiter','\t','precision',12)
			// ctrs2=num2str(ctrs2);
		}
		double ctrs3 = simtime * dt;
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

	void printArray(double** arr, int n, string name)
	{
		cout << "*********************************************************************************************" << endl;
		cout << "\t\t\t" << name << " Array"<<endl;
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				cout << arr[i][j] << "\t";
			}
			cout << endl;
		}
		cout << "*********************************************************************************************" << endl;
	}
};