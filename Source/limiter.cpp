#include<iostream>
#include<cmath>


class limiter
{
public:
	limiter()
	{

	}
	// Function Definitions

	// Function Definitions

	//
	// Arguments    : double n
	//                const double f[1764]
	//                emxArray_real_T *df1
	//                emxArray_real_T *df2
	// Return Type  : void
	//
	void flimiter(double n, double f[1764], double df1[1764],
		double df2[1764])
	{
		int k;
		int loop_ub;
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

		loop_ub = (int)n * (int)n;
		for (k = 0; k < loop_ub; k++) {
			df1[k] = 0.0;
			df2[k] = 0.0;
		}



		/*for (k = 0; k < (int)n; k++) {
			df1->data[((int)n + df1->size[0] * k) - 1] = 0.0;
			df1->data[k + df1->size[0] * ((int)n - 1)] = 0.0;
			df2->data[((int)n + df2->size[0] * k) - 1] = 0.0;
			df2->data[k + df2->size[0] * ((int)n - 1)] = 0.0;
		}*/

		for (loop_ub = 1; loop_ub - 1 < (int)((n - 1.0) + -1.0); loop_ub++) {
			for (k = 1; k - 1 < (int)((n - 1.0) + -1.0); k++) {
				df1x = f[((int)((2.0 + (double)(loop_ub - 1)) + 1.0) + 42 * k) - 1] -
					f[loop_ub + 42 * k];
				df2x = f[loop_ub + 42 * k] - f[((int)((2.0 + (double)(loop_ub - 1)) - 1.0)
					+ 42 * k) - 1];
				df1y = f[loop_ub + 42 * ((int)((2.0 + (double)(k - 1)) + 1.0) - 1)] -
					f[loop_ub + 42 * k];
				df2y = f[loop_ub + 42 * k] - f[loop_ub + 42 * ((int)((2.0 + (double)(k - 1))
					- 1.0) - 1)];

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

				//
				//  if(df1x*df2x < 0),
				//          f1=0;
				//    else
				//          s=sign(df1x);
				//          a=abs(df1x);
				//          b=abs(df2x);
				//          f1=s*min(max([a b]),2.0*min([a b]));
				//  end
				//
				//  if(df1y*df2y < 0),
				//          f2=0;
				//    else
				//          s=sign(df1y);
				//          a=abs(df1y);
				//          b=abs(df2y);
				//          f2=s*min(max([a b]),2.0*min([a b]));
				//  end
				df1[loop_ub] = df2x;
				df2[loop_ub] = df1x;
			}
		}
	}
};