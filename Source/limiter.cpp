#include<iostream>
#include<cmath>
#include<algorithm>

using namespace std;


class limiter
{
public:
	// Function Definitions

	//
	// Arguments    : double n
	//                double f[1764]
	//                double f[1764]
	//                double df2[1764]
	// Return Type  : void
	//
	void flimiter(int n, double** f, double** df1,
		double** df2)
	{

		double df1x;
		double df2x;
		double df1y;
		double df2y;
		double s;
		double a, b,c;
		double f1,f2;
		
		//zero array
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++)
			{
				df1[i][j] = 0.0;
				df2[i][j] = 0.0;
			}
		}

		for (int j = 0; j < n; j++)
		{
			df1[1][j] = 0.0;
			df2[j][1] = 0.0;
			df1[1][j] = 0.0;
			df2[j][1] = 0.0;
		}

		for (int j = 0; j < n; j++)
		{
			df1[n - 1][j] = 0.0;
			df2[j][n - 1] = 0.0;
			df1[n - 1][j] = 0.0;
			df2[j][n - 1] = 0.0;
		}

		for (int j = 1; j < n - 1; j++) {
			for (int k = 1; k < n - 1; k++)
			{
				df1x = f[j + 1][k] - f[j][k];
				df2x = f[j][k] - f[j - 1][k];
				df1y = f[j][k + 1] - f[j][k];
				df2y = f[j][k] - f[j][k - 1];

				if (df1x * df2x < 0.0) {
					f1 = 0.0;
				}
				else {
					s = sign(df1x);
					a = abs(df1x);
					b = abs(df2x);
					c = 0.5 * (a + b);
					auto min = minmax({ 2 * a, 2 * b, c });
					f1 = s * min.first;
				}
				if (df1y * df2y < 0.0) {
					f2 = 0.0;
				}
				else {
					s = sign(df1y);
					a = abs(df1y);
					b = abs(df2y);
					c = 0.5 * (a + b);
					auto min = minmax({ 2 * a, 2 * b, c });
					f2 = s * min.first;
				}
				df1[j][k] = f1;
				df2[j][k] = f2;
			}
		}
	}
	///
	int sign(double v)
	{
		return (v < 0) ? -1 : (v > 0) ? 1 : 0;
	}
};