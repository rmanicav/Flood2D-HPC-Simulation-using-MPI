using namespace std;
class dem
{
public:
	int ncols, nrows;
	double xllCorner, yllCorner, cellSize, xloc, yloc;
	int nodata;
	double** dm;
	void readDemFile(string fileName);
	void readLocData(string fileName, double& xloc, double& yloc);
private:
	double** allocateMemory(int n, int m);
};
/// <summary>
	/// 
	/// </summary>
	/// <param name="fd"></param>
	/// <param name="fileName"></param>
void dem::readDemFile(string fileName)
{
	ifstream infile(fileName);
	string line;
	if (!infile.good())
	{
		std::cout << "Error Reading DEM file" << endl;

	}
	else
	{
		string name;
		infile >> name;
		infile >> ncols;
		infile >> name;
		infile >> nrows;
		infile >> name;
		infile >> xllCorner;
		infile >> name;
		infile >> yllCorner;
		infile >> name;
		infile >> cellSize;
		infile >> name;
		infile >> nodata;

			dm = allocateMemory(nrows,ncols);
			for (int i = 0; i < nrows ; i++)
			{
				for (int j = 0; j < ncols ; j++)
				{
					infile >> dm[i][j];
					cout << dm[i][j] << "\t";
				}
				cout << endl;
			}
			infile.close();
	}
}

/// <summary>
/// read location data
/// </summary>
/// <param name="fd"></param>
/// <param name="fileName"></param>
void dem::readLocData(string fileName, double& xloc, double& yloc)
{
	ifstream infile(fileName);
	if (!infile.good())
	{
		std::cout << "Error Reading Location file" << endl;
	}
	else
	{
		string line;
		getline(infile, line);
		std::cout << "************************************************" << endl;
		std::cout << "Inflow Location" << endl;
		getline(infile, line,',');
		xloc = stod(line);
		getline(infile, line, ',');
		yloc = stod(line);
		std::cout << "xLocation:" << xloc << "\t" << "xLocation:" << yloc << endl;
		std::cout << "************************************************" << endl;
		infile.close();
	}

}

/// <summary>
	/// 
	/// </summary>
	/// <param name="n"></param>
	/// <param name="m"></param>
	/// <returns></returns>
double** dem::allocateMemory(int n, int m)
{
	double** arr;
	arr = (double**)malloc(sizeof(double*) * (int)n);
	for (int i = 0; i < n; i++)
	{
		arr[i] = (double*)malloc(sizeof(double) * (int)m);
	}
	return arr;
}