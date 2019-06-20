#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<ostream>
#include<string>
#include<iomanip>
#include<ctime>
using namespace std;

/// <summary>
/// 
/// </summary>
struct floodData
{
	int  L, dim,cellSize,x,y,z;
	double hextra, epsilon, nt, ntPlot, dt, initWSE, hWL,cr,gravity, manN, initV;
	double **zc;
};

/// <summary>
/// 
/// </summary>
class helper
{
public:

	/// <summary>
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
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
	/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <returns></returns>
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
	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
	void freeMemory(double** arr,int n)
	{
		for (int i = 0; i < n; i++)
		{
			free(arr[i]);
		}
		free(arr);
	}
	/// <summary>
	/// 
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
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
	/// <summary>
	/// 
	/// </summary>
	/// <param name="a"></param>
	/// <param name="b"></param>
	/// <param name="c"></param>
	/// <returns></returns>
	double*** allocate3dMemory(int a,int b, int c)
	{
		double*** arr = new double** [a];
		
		for (int i = 0; i < b; i++)
		{
			arr[i] = new double* [b];
			for (int j = 0; j < b; j++)
			{
				arr[i][j] = new double[b];
			}			
		}
		for (int i = 0; i < b; i++)
		{
			for (int j = 0; j < b; j++)
			{
				for (int k = 0; k < b; k++)
				{
					arr[i][j][k] = 0.0;
				}
				
			}
		}

		return arr;
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="fd"></param>
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
		getline(infile, name);
		pos = name.find(":");
		fd->x = atoi(name.substr(pos + 1).c_str());
		cout << "X:" << fd->x<< endl;
		getline(infile, name);
		pos = name.find(":");
		fd->y = atoi(name.substr(pos + 1).c_str());
		cout << "Y:" << fd->y << endl;
		getline(infile, name);
		pos = name.find(":");
		fd->z = atoi(name.substr(pos + 1).c_str());
		cout << "Z:" << fd->z << endl;
		cout << "\n-------------------------------\n\n";
		infile.close();
		infile.open("Input/zc_dem.txt");
		fd->zc = allocateMemory(fd->x);
		clearArray(fd->zc,fd->x);
		for (int i = 0; i < fd->x; i++)
		{			
	      for (int j = 0; j < fd->y; j++)
			{
				infile >> fd->zc[i][j];
			}
		}
	//	printArray(fd->zc, 42, "zc");
		infile.close();
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="h"></param>
	/// <param name="n"></param>
	/// <param name="fileName"></param>
	void writeOutputFile(double **h,int n, string fileName)
	{
		ofstream outStream;
		outStream.open("Output/" + fileName,ios::out);
		if (outStream.is_open())
		{
			
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					
					outStream <<setw(5)<< h[i][j] << "\t";
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

	void write3dOutputFile(double*** arr, int n, string fileName)
	{
		ofstream outStream;
		outStream.open("Output/" + fileName, ios::out);
		if (outStream.is_open())
		{

			outStream << "*********************************************************************************************" << endl;
			outStream << " Arr 0 dim" << endl;
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
	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
	/// <param name="name"></param>
	void printArray(double** arr, int n, string name)
	{
		cout << "*********************************************************************************************" << endl;
		cout << "\t\t\t" << name << " Array"<<endl;
		for (int i = 0; i < n; i++)
		{
			cout <<i + 1 << "\t";
		}
		cout << endl;
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
				cout << arr[2][j][k] << "\t";
			}
			cout << endl;
		}
		cout << endl;
		cout << "*********************************************************************************************" << endl;
		cout << "*********************************************************************************************" << endl;
	}
	/// <summary>
/// <summary>
/// 
/// 
/// </summary>
/// <param name="h"></param>
/// <param name="u"></param>
/// <param name="v"></param>
/// <param name="n"></param>
/// <param name="interval"></param>
	void printArrayInterval(double** h,double **u,double **v, int n,int iterNum)
	{
		if (iterNum % 2 == 10)
		{
			cout << "*********************************************************************************************" << endl;
			cout<<"**************************Iteration "<<iterNum<<" **********************************************"<<endl;
			cout << "\t\t\t h Array" << endl;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					cout << h[i][j] << "\t";
				}
				cout << endl;
			}
			cout << "*********************************************************************************************" << endl;
			cout << "*********************************************************************************************" << endl;
			cout << "\t\t\t u Array" << endl;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					cout << h[i][j] << "\t";
				}
				cout << endl;
			}
			cout << "*********************************************************************************************" << endl; cout << "*********************************************************************************************" << endl;
			cout << "\t\t\t v Array" << endl;
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					cout << h[i][j] << "\t";
				}
				cout << endl;
			}
			cout << "*********************************************************************************************" << endl;
		}
	}

	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="error"></param>
	void writeErrorLog(string error)
	{
		ofstream outStream;
		time_t now = time(0);
		//char* dt = ctime(&now);

		outStream.open("Output/errorLog.txt", ios::out);
		if (outStream.is_open())
		{
			outStream << "*********************************************************" << endl;
			//outStream << "Date is: " << dt << endl;
			outStream << endl;
			outStream << "Error Message:" << error << endl;
			outStream << "*********************************************************" << endl;
			outStream.close();
		}
		else
		{
			cout << "Unable to open file" << endl;
		}
	}
	/// <summary>
	/// 
	/// </summary>
	/// <param name="h"></param>
	/// <param name="ntPlot"></param>
	/// <param name="simTime"></param>
	/// <param name="dt"></param>
	void writeSensor(double **h,double ntPlot,int simTime, double dt)
	{
		ofstream outStream,outStream2,outStream3;
		double ctrs, ctrs2, ctrs3;

		
		ctrs = simTime * dt;
		double** hsens_1 = allocateMemory(1);
		double** hsens_2 = allocateMemory(1);
		double** hsens_3 = allocateMemory(1);
		hsens_1[0][0] = h[24][28];
		hsens_2[0][0] = h[29][28];
		hsens_3[0][0] = h[36][28];

			ctrs = simTime * dt;
			if (fmod(ctrs, 0.5) == 0)
			{
				outStream.open("Output/hsensor1.txt", ios::out);
				outStream << hsens_1[0][0] << endl;
				outStream.close();
			}
			ctrs2 = simTime * dt;
			if (fmod(ctrs2, 0.5) == 0)
			{
				outStream2.open("Output/hsensor2.txt", ios::out);
				outStream2 << hsens_2[0][0] << endl;
				outStream2.close();
			}
			ctrs3 = simTime * dt;
			if (fmod(ctrs3, 0.5) == 0)
			{
				outStream3.open("Output/hsensor3.txt", ios::out);
				outStream3 << hsens_2[0][0] << endl;
				outStream3.close();
			}
			freeMemory(hsens_1,1);
			freeMemory(hsens_2, 1);
			freeMemory(hsens_3, 1);
	}
};