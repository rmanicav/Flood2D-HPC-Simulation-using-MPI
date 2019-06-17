#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<ostream>
#include<iomanip>
using namespace std;


struct floodData
{
	int iv0[31], iv1[2];
	double gravity, manN,  initV;
	int  L, dim,cellSize;
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
	///
	void freeMemory(double** arr,int n)
	{
		for (int i = 0; i < n; i++)
		{
			free(arr[i]);
		}
		free(arr);
	}
	void freeMemory3d(double*** arr, int n)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				free(arr[i][j]);
			}
			free(arr[i]);
		}
		free(arr);
	}
	///
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
	void writeHout(double **h,int n)
	{
		//hsens_1 = h(25, 29); hsens_2 = h(30, 29); hsens_3 = h(37, 29);
		ofstream outStream;;
		outStream.open("C:/Users/raj/source/repos/FVM/output/hOutput1.txt",ios::out);
		if (outStream.is_open())
		{
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					
					outStream << h[i][j] << "\t";
				}
				outStream << endl;
			}
			outStream.close();
		}
		else
		{
			cout << "Unable to open file" << endl;
		}
		
	}
	///

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