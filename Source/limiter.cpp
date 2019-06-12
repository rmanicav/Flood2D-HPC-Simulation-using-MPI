#include<iostream>
#include<cmath>


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
		double a;
		double varargin_1[3];
		int idx;
		int b_k;
		bool exitg1;
		
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
			df1[n][j] = 0.0;
			df2[j][n] = 0.0;
			df1[n][j] = 0.0;
			df2[j][n] = 0.0;
		}
		
		for (int j = 1; j < n; j++) {
			for (int k = 1; k < n ; k++)
			{
				df1x = f[j + 1][k] - f[j][k];
				df2x = f[j][k] - f[j -1][k];
				df1y = f[j][k + 1] - f[j][k];
				df2y = f[j][k] - f[j][k - 1];
	
				
				//  Superbee limiterv  (df1x,df2x,df1m,df1y,df2y,df2m)
				if (df1x * df2x < 0.0) {
					df2x = 0.0;
				}
				else {
					s = df1x;
					if (df1x < 0.0) {
						s = -1.0;
					}
					else if (df1x > 0.0) {
						s = 1.0;
					}
					else {
						if (df1x == 0.0) {
							s = 0.0;
						}
					}

					a = std::abs(df1x);
					df1x = std::abs(df2x);
					varargin_1[0] = 2.0 * a;
					varargin_1[1] = 2.0 * df1x;
					varargin_1[2] = 0.5 * (a + df1x);
					if (varargin_1[0] != NAN) {
						idx = 1;
					}
					else {
						idx = 0;
						b_k = 2;
						exitg1 = false;
						while ((!exitg1) && (b_k < 4)) {
							if (varargin_1[b_k - 1] != NAN) {
								idx = b_k;
								exitg1 = true;
							}
							else {
								b_k++;
							}
						}
					}

					if (idx == 0) {
						df1x = varargin_1[0];
					}
					else {
						df1x = varargin_1[idx - 1];
						while (idx + 1 < 4) {
							if (df1x > varargin_1[idx]) {
								df1x = varargin_1[idx];
							}

							idx++;
						}
					}

					df2x = s * df1x;
				}

				if (df1y * df2y < 0.0) {
					df1x = 0.0;
				}
				else {
					s = df1y;
					if (df1y < 0.0) {
						s = -1.0;
					}
					else if (df1y > 0.0) {
						s = 1.0;
					}
					else {
						if (df1y == 0.0) {
							s = 0.0;
						}
					}

					a = std::abs(df1y);
					df1x = std::abs(df2y);
					varargin_1[0] = 2.0 * a;
					varargin_1[1] = 2.0 * df1x;
					varargin_1[2] = 0.5 * (a + df1x);
					if (varargin_1[0] != NAN) {
						idx = 1;
					}
					else {
						idx = 0;
						b_k = 2;
						exitg1 = false;
						while ((!exitg1) && (b_k < 4)) {
							if (varargin_1[b_k - 1] != NAN) {
								idx = b_k;
								exitg1 = true;
							}
							else {
								b_k++;
							}
						}
					}

					if (idx == 0) {
						df1x = varargin_1[0];
					}
					else {
						df1x = varargin_1[idx - 1];
						while (idx + 1 < 4) {
							if (df1x > varargin_1[idx]) {
								df1x = varargin_1[idx];
							}

							idx++;
						}
					}

					df1x *= s;
				}
				df1[j][k] = df2x;
				df2[j][k] = df1x;
			}
		}
	}
};