#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<ostream>
#include<string>
#include<iomanip>
#include<vector>
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
		string name;
		infile.open("Input/iv.txt");
		cout << "----------------------------------\n";
		getline(infile, name);
		size_t pos = name.find(":");
		fd->gravity = atof(name.substr(pos+1).c_str());
		cout << "Gravity:" << fd->gravity << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->manN = atof(name.substr(pos + 1).c_str());
		cout << "Manning:" << fd->manN << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->hextra = atof(name.substr(pos + 1).c_str());
		cout << "hextra:" << fd->hextra << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->epsilon = atof(name.substr(pos + 1).c_str());
		cout << "epsilon:" << fd->epsilon << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->cellSize = atoi(name.substr(pos + 1).c_str());
		cout << "Cell size:" << fd->cellSize << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->L = atoi(name.substr(pos + 1).c_str());
		cout << "Box size:" << fd->L << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->nt = atof(name.substr(pos + 1).c_str());
		cout << "Number of time steps:" << fd->nt << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->ntPlot = atof(name.substr(pos + 1).c_str());
		cout << "Plotting interval:" << fd->ntPlot << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->dt = atof(name.substr(pos + 1).c_str());
		cout << "Time steps:" << fd->dt << endl;		
		getline(infile, name);
		pos = name.find(":");
		fd->initV = atof(name.substr(pos + 1).c_str());
		cout << "Init value:" << fd->initV << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->initWSE= atof(name.substr(pos + 1).c_str());
		cout << "WSE downstream value:" << fd->initWSE << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->hWL= atof(name.substr(pos + 1).c_str());
		cout << "WSE upstream value:" << fd->hWL << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->cr = atof(name.substr(pos + 1).c_str());
		cout << "cr value:" << fd->cr << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->dim= atoi(name.substr(pos + 1).c_str());
		cout << "Dimension:" << fd->dim << endl;
		cout << "\n-------------------------------\n\n";
		infile.close();
		infile.open("Input/zc_dem.txt");
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
		infile.close();
	}
	void writeHout(double **h,int n, string fileName)
	{
		//hsens_1 = h(25, 29); hsens_2 = h(30, 29); hsens_3 = h(37, 29);
		ofstream outStream;;
		outStream.open(fileName,ios::out);
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