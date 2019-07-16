#include<cmath>
#include<iostream>
#include "omp.h"
using namespace std;

/// <summary>
/// 
/// </summary>
class corrector
{
public:
	
	/// <summary>
	/// 
	/// </summary>
	/// <param name="U"></param>
	/// <param name="F"></param>
	/// <param name="G"></param>
	/// <param name="n"></param>
	/// <param name="dt2"></param>
	/// <param name="dt"></param>
	/// <param name="sox"></param>
	/// <param name="sfx"></param>
	/// <param name="soy"></param>
	/// <param name="sfy"></param>
	/// <param name="grav"></param>
	/// <returns></returns>
	double*** fcorrector(double ***U, double ***F,  double ***G,double n, double dt2, double dt, double **sox, 
		double **sfx, double **soy, double **sfy,double grav)
	{
         int nThreads =3;
         double wtime = omp_get_wtime();
         #pragma omp parallel num_threads(nThreads)
          {
      		int j;
		int k;
                int tid = omp_get_thread_num();
                printf("Thread number is %d\n",tid);
                #pragma omp sections 
				{
                 #pragma omp section
					for (j = 0; j < n - 1; j++) {
						for (k = 0; k < n - 1; k++) {

							U[0][j][k] = U[0][j][k] - dt2 * (F[0][j + 1][k] - F[0][j][k]) - dt2 * (G[0][j][k + 1] - G[0][j][k]);
							if (U[0][j][k] < 0.0) {
								U[0][j][k] = 0.0;
							}
						}
					}

                  #pragma omp section
					for (j = 0; j < n - 1; j++) {
						for (k = 0; k < n - 1; k++) {
							U[1][j][k] = ((U[1][j][k] - dt2 * (F[1][j + 1][k] - F[1][j][k])) - dt2 * (G[1][j][k + 1] - G[1][j][k])) - (dt * grav * (sox[j][k] + sfx[j][k]));
						}
					}
					
                   #pragma omp section
					for (j = 0; j < n - 1; j++) {
						for (k = 0; k < n - 1; k++) {
							U[2][j][k] = ((U[2][j][k] - dt2 * (F[2][j + 1][k] - F[2][j][k]) - dt2 * (G[2][j][k + 1] - G[2][j][k]))) - (dt * grav * (soy[j][k] + sfy[j][k]));
						}
					}
				}

	}
        wtime = omp_get_wtime() - wtime;
        printf( "Time taken for corrector is  %f\n",wtime);
       	return U;
         
	}
	
};
