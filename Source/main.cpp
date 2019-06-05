#include<iostream>
#include "Include/corrector.h"
#include "Include/limit.h"
#include "Include/limiter.h"
#include "Include/predictor.h"
#include "Include/slope.h"
#include "Include/solver.h"
#include "Include/fluxes.h"
int main(void)
{
	int i0;
	int j;
	static double zc[1764];
	int k;
	static double wse[1764];
	static const signed char iv0[31] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
	  13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 34, 35, 36, 37, 38, 39, 40, 41 };

	static double U[5292];
	static double h[1764];
	static double u[1764];
	static double v[1764];
	int simtime;
	static double dzcx[1764];
	static double dzcy[1764];
	double dwsex[1764];
	double dwsey[1764];
	double dux[1764];
	double duy[1764];
	double dvx[1764];
	double dvy[1764];
	double sox[1764];
	double soy[1764];
	double sfx[100];
	double sfy[100];
	static int iv1[2] = { 10, 10 };

	double wsep[1764];
	double up[1764];
	double vp[1764];
	signed char fileid;
	static double unusedExpr[5292];

	//
	grav = 9.806;

	//  ManN=0.00;
	//  minimum depth of water (threshold for dry bed)(1cm)
	//  Cell size(dx=dy)m
	//  box size(200*200m)
	//  (n=40) number of cells (nx=ny)
	//  number of edges   (will be used to compute fluxes at the interface)                     
	//  x=0:cellsize:L;                         % array for edge coordinates in x    
	//  y=0:cellsize:L;                         % array for edge coordinates in y
	//  x coordinate of cell center
	//  y coordinate of cell center
	//  number of time steps
	//  plotting interval
	//  time steps
	//  ratio of dt/dx or dt/dy
	//  Define bed elevation
	for (i0 = 0; i0 < 1764; i0++) {
		zc[i0] = 1.0;
	}

	for (j = 0; j < 42; j++) {
		for (k = 0; k < 42; k++) {
			zc[42 * k] = 99999.0;
			zc[41 + 42 * k] = 99999.0;
			zc[j] = 99999.0;
			zc[1722 + j] = 99999.0;
			for (i0 = 0; i0 < 23; i0++) {
				zc[20 + 42 * i0] = 99999.0;
			}

			for (i0 = 0; i0 < 8; i0++) {
				zc[20 + 42 * (34 + i0)] = 99999.0;
			}

			for (i0 = 0; i0 < 23; i0++) {
				zc[21 + 42 * i0] = 99999.0;
			}

			for (i0 = 0; i0 < 8; i0++) {
				zc[21 + 42 * (34 + i0)] = 99999.0;
			}
		}
	}

	//  Set up initial conditions
	//  wse=1.01*ones(n,n);                           % Initial water surface elevation  
	for (i0 = 0; i0 < 1764; i0++) {
		wse[i0] = 6.0;
	}

	for (i0 = 0; i0 < 42; i0++) {
		for (j = 0; j < 21; j++) {
			wse[j + 42 * i0] = 11.0;
		}
	}

	//  0.5*n = 20 or higher water level
	for (i0 = 0; i0 < 31; i0++) {
		for (j = 0; j < 2; j++) {
			wse[(j + 42 * iv0[i0]) + 20] = 99999.0;
		}
	}

	//  0.5*n = 20 or higher water level
	//  water depth to water surface elevatin relation ship
	for (j = 0; j < 42; j++) {
		for (k = 0; k < 42; k++) {
			h[j + 42 * k] = wse[j + 42 * k] - zc[j + 42 * k];
			h[42 * k] = 0.0;
			h[41 + 42 * k] = 0.0;
			h[j] = 0.0;
			h[1722 + j] = 0.0;
			if (h[j + 42 * k] < 0.0) {
				h[j + 42 * k] = 0.0;
			}
		}
	}

	// Intialize the arrays
	memset(&U[0], 0, 5292U * sizeof(double));

	//  fluxes in the $x$ direction
	//  fluxes in the $y$ direction
	//  friction slope Sf
	//  Bed slope So
	//  variables on the new time level
	for (i0 = 0; i0 < 42; i0++) {
		for (j = 0; j < 42; j++) {
			U[j + 42 * i0] = h[j + 42 * i0];
			U[1764 + (j + 42 * i0)] = 0.0 * h[j + 42 * i0];
			U[3528 + (j + 42 * i0)] = 0.0 * h[j + 42 * i0];
			u[j + 42 * i0] = U[1764 + (j + 42 * i0)] / (U[j + 42 * i0] + 0.1);
			v[j + 42 * i0] = U[3528 + (j + 42 * i0)] / (U[j + 42 * i0] + 0.1);
		}
	}

	//  Bed slope along X and Y
	//  writes the outputs for the sensors
	//  hnorm=zeros(size(num));     hsens_1=zeros(m,3); hsens_2=zeros(m,3); hsens_3=zeros(m,3); 
	cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor1.txt",
		"wb");
	cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor2.txt",
		"wb");
	cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor3.txt",
		"wb");

	//  ******************************************************************************************************************************************          
	//  tcntr=0;
	for (j = 0; j < 20000; j++) {
		simtime = 1 + j;
		limiter(zc, dzcx, dzcy);
		limiter(wse, dwsex, dwsey);
		limiter(u, dux, duy);
		limiter(v, dvx, dvy);

		//  Slope calculation
		slope(h, u, v, dzcx, dzcy, sox, soy, sfx, sfy);

		//  predictor step (estimate the values at half timestep)
		predictor(wse, h, u, v, dwsex, dwsey, dux, duy, dvx, dvy, dzcx, dzcy, zc,
			sox, sfx, iv1, soy, sfy, iv1, wsep, up, vp);

		//    Compute fluxes at the interfaces
		// [F, G, amax] = fluxes(UP,n,dwsex,dwsey,dux,duy,dvx,dvy,hextra,zc);
		//  Estimate the flux vectors on the next time step
		corrector(U, sox, sfx, iv1, soy, sfy, iv1, grav, unusedExpr);

		//  U(:,:,:)=Unew(:,:,:);
		//  h(:,:)= Unew(:,:,1);                    % computed water depth (water level)      
		//  if Unew(:,:,1)< epsilon
		memset(&u[0], 0, 1764U * sizeof(double));
		memset(&v[0], 0, 1764U * sizeof(double));

		// else
		// u(:,:)= Unew(:,:,2)./(Unew(:,:,1)+hextra);
		// v(:,:)= Unew(:,:,3)./(Unew(:,:,1)+hextra);
	}

	// wse = Unew(:,:,1)+zc;                  % compute the new free surface height 
	// if wse < 0
	//    wse = 0;
	// end
	//      tcntr=tcntr+dt;
	//    dt= cr*cellsize/amax;
	//   cr=amax*dt/cellsize;
	//  Writes the output at each counters seconds
	if (b_mod((double)simtime * 0.05) == 0.0) {
		fileid = cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hOut_",
			"wb");

		//    fprintf(fid1,'%.2f\r\n',h(:,:));
		//   dlmwrite(outfile,x,'delimiter','\t','precision',12)
		b_fclose((double)fileid);

		//  % %   outputs the 2D-water depth plot
		//  ctrs=sprintf('L%s',ctrs);
		// saveas(gca,output);
	}

	if (c_mod((double)simtime * 0.05) == 0.0) {
		fileid = cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor1.txt",
			"ab");
		b_fclose((double)fileid);

		// dlmwrite(outfile2,hsens_1,'-append','newline','pc','delimiter','\t','precision',12) 
		// ctrs2=num2str(ctrs2);
	}

	if (c_mod((double)simtime * 0.05) == 0.0) {
		fileid = cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor2.txt",
			"ab");

		// dlmwrite(outfile3,hsens_2,'-append','newline','pc','delimiter','\t','precision',12) 
		b_fclose((double)fileid);

		//  ctrs3=num2str(ctrs3);
	}

	if (c_mod((double)simtime * 0.05) == 0.0) {
		fileid = cfopen("C:\\Users\\raj\\source\\repos\\Flood2dOutput\\hsensor3.txt",
			"ab");

		// dlmwrite(outfile4,hsens_3,'-append','newline','pc','delimiter','\t','precision',12) 
		b_fclose((double)fileid);

		// ctrs4=num2str(ctrs4);
	}

	//
	//      k=tstep*dt;
	//      if k==7.50,
	//         break;
	//      end
	//      fprintf('Time Step = %d, Courant counter tcntr = %g \n',dt,tcntr)
	//    fprintf('Courant value = %d, Time Step = %g \n',cr,t)
	//      if (mod(tstep,ntplot) == 0),
	//          surf(xc,yc,h')
	//          meshgrid(cellsize/2:cellsize:L, cellsize/2:cellsize:L, 5:0.5:10);
	//          xlabel('xc (m)'); ylabel('yc (m)'); zlabel('h (m)')
	//          title(['Flood wave propagation downstreams of the Labscale dam at ',num2str(t,'%02g'),' seconds ',date]); 
	//  % %         view(-150,150)
	//          view(0,90)
	//          hcolr=colorbar;
	//          q2=get(hcolr,'Title'); titlesg='Water Depth'; set(q2,'string',titlesg); 
	//          pause(0.1)
	//      end
}

//
// File trailer for fvm2d.cpp
//
// [EOF]
//
	return 0;
}