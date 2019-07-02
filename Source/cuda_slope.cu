#include<cmath>
#include<iostream>
#include<string>
#include<stdio.h>
#include<stdlib.h>
#include<cuda.h>
using namespace std;
/// <summary>
/// 
/// 
/// </summary>

__global__ void slopeKernel(double* h, double* u, double* v,
		double ManN, double hextra, int n, double *sfx, double* sfy, double *sox,double* soy,int cellsize,double *dzcx, double* dzcy)
{
 int id = blockIdx.x * blockDim.x + threadIdx.x;

          if(id < n)
          {
                      sox[id] = dzcx[id] / cellsize;
                      soy[id] = dzcy[id] / cellsize;

		//  Intermediate frction slope values
				if (h[id] < 1.0) {
					sfx[id] = 0.0;
					sfy[id] = 0.0;
				}
				else {
					sfx[id] = u[id] * (ManN * ManN) * hypot(u[id], v[id]) / pow(h[id] + hextra,
							1.3333333333333333);
					sfy[id] = v[id] * (ManN * ManN) * hypot(u[id], v[id]) / pow(h[id] + hextra,					1.3333333333333333);
				}
          
         }
}
	//
	/// <summary>
	/// 
	/// </summary>
	/// <param name="h"></param>
	/// <param name="u"></param>
	/// <param name="v"></param>
	/// <param name="ManN"></param>
	/// <param name="hextra"></param>
	/// <param name="dzcx"></param>
	/// <param name="dzcy"></param>
	/// <param name="cellsize"></param>
	/// <param name="n"></param>
	/// <param name="sox"></param>
	/// <param name="soy"></param>
	/// <param name="sfx"></param>
	/// <param name="sfy"></param>
extern "C" void fslope(double** h, double** u, double** v,
                double ManN, double hextra, double** dzcx, double**
                dzcy, int cellsize, int n, double** sox, double** soy
                , double **sfx, double** sfy)
        {

              double* d_v;
              double* d_h;
              double* d_u;
              double* d_sox;
              double* d_soy;
              double* d_sfx;
              double* d_sfy;
              double* d_zcx;
              double* d_zcy;
         
              int size = n * n * sizeof(double);
             
              //allocate memory
              cudaMalloc(&d_h,size);
              cudaMalloc(&d_u,size);
              cudaMalloc(&d_v,size);
              cudaMalloc(&d_sox,size);
              cudaMalloc(&d_soy,size);
              cudaMalloc(&d_sfx,size);
              cudaMalloc(&d_sfy,size);
              cudaMalloc(&d_zcx,size);
              cudaMalloc(&d_zcy,size);
              

              //copy host to device
                cudaMemcpy(d_h,h,size,cudaMemcpyHostToDevice);
                cudaMemcpy(d_u,u,size,cudaMemcpyHostToDevice);
                cudaMemcpy(d_v,v,size,cudaMemcpyHostToDevice);
                cudaMemcpy(d_zcx,dzcx,size,cudaMemcpyHostToDevice);
                cudaMemcpy(d_zcy,dzcy,size,cudaMemcpyHostToDevice);

                dim3 grid(1*1);//number of blocks
                dim3 block(n*n);//number of threads

                slopeKernel<<<grid,block>>>(d_h,d_u, d_v,ManN,hextra,n,d_sfx,d_sfy,d_sox,d_soy,cellsize,d_zcx,d_zcy);

                //copy device to host
                cudaMemcpy(sfx,d_sfx,size,cudaMemcpyDeviceToHost);
                cudaMemcpy(sfy,d_sfy,size,cudaMemcpyDeviceToHost);
                cudaMemcpy(sox,d_sox,size,cudaMemcpyDeviceToHost);
                cudaMemcpy(soy,d_soy,size,cudaMemcpyDeviceToHost);
                              
                //print
                cout<<"SFY"<<endl;
                for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                       cout<<sfx[i][j]<<endl;
                        }
                cout<<endl;
                }
                cout<<"SFY"<<endl;
                for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                       cout<<sfy[i][j]<<endl;
                        }
                cout<<endl;
                }
                
                cout<<"SOX"<<endl;
                for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                       cout<<sox[i][j]<<endl;
                        }
                cout<<endl;
                }
                cout<<"SOY"<<endl;
                for (int i = 0; i < n; i++) {
                        for (int j = 0; j < n; j++) {
                       cout<<soy[i][j]<<endl;
                        }
                cout<<endl;
                }
                
                //free memory
                cudaFree(d_h);
                cudaFree(d_u);
                cudaFree(d_v);
                cudaFree(d_sox);
                cudaFree(d_soy);
                cudaFree(d_sfx);
                cudaFree(d_sfy);
                cudaFree(d_zcx);
                cudaFree(d_zcy);

                 

}
//
// File trailer for slope.cpp
//
// [EOF]
//
