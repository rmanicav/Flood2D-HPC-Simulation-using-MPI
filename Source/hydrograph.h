#include<string>
using namespace std;
struct hyg
{
	double* dtime;
	double* discharge;
};
///
///
/// 
 	
class hydrograph
	{
	public:
		void readHygData(hyg *h,string fileName);		
		void printHyg(hyg * h,size_t lCount);
		void clearArray(double* arr, size_t n);
	};
	/// <summary>
	/// Read hydrograph data
				/// </summary>
				/// <param name="fd"></param>
				/// <param name="fileName"></param>
	void hydrograph :: readHygData(hyg* hg,string fileName)
	{
		ifstream infile(fileName.c_str());
		int i = 0;
		size_t lCount;
		string line;

		if (!infile.good())
		{
			std::cout << "Error Reading Hydrograph file" << endl;
		}
		else
		{
			
			lCount = std::count(std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>(), '\n');
			infile.clear();
			infile.seekg(0, ios::beg);
			lCount = lCount - 2;
			hg->dtime = (double *)malloc(sizeof(double) * lCount);
			clearArray(hg->dtime,lCount);
			hg->discharge = (double *)malloc(sizeof(double) * lCount);
			clearArray(hg->discharge,lCount);
			while (getline(infile,line))
			{
				if (line.size() > 0 && line[0] != '%')
				{
					size_t pos = 0;
					pos = line.find(',');
					hg->dtime[i] = stod(line.substr(0,pos));
					hg->discharge[i] = stod(line.substr(pos + 1, line.length() - 1));
					line.clear();
					i++;
				}
			}	
		}
		infile.close();
		printHyg(hg, lCount);
		free(hg->dtime);
		free(hg->discharge);
	}
	
	/// <summary>
	/// 
	/// </summary>
	void hydrograph::printHyg(hyg *h,size_t lCount)
	{
		cout << "********************************************************" << endl;
		cout << "\tHydrograph" << endl;
		cout << "Time\tDischarge" << endl;
		for (int i = 0; i < lCount; i++)
		{
			cout << h->dtime[i] << "," << h->discharge[i]<< endl;
		}
		cout << "********************************************************" << endl;
	}


	/// <summary>
	/// 
	/// </summary>
	/// <param name="arr"></param>
	/// <param name="n"></param>
	void hydrograph:: clearArray(double * arr, size_t n)
	{
			for (int i = 0; i < n; i++)
			{
				arr[i]= 0.0;
			}
	}