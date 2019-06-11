#include "solver.cpp"
#include <stdlib.h>
#include "fluxes_emxAPI.cpp"
class fluxes
{

public:
	fluxes()
	{

	}
	// Function Definitions

	//
	// Compute fluxes in X-direction
	// Arguments    : const double UP[5292]
	//                double n
	//                const double dwsex[1764]
	//                const double dwsey[1764]
	//                const double dux[1764]
	//                const double duy[1764]
	//                const double dvx[1764]
	//                const double dvy[1764]
	//                double hextra
	//                const double zc[1764]
	//                emxArray_real_T *F
	//                emxArray_real_T *G
	//                double *amax
	// Return Type  : void
	//
	void ffluxes(const double UP[5292], double n, const double dwsex[1764], const
		double dwsey[1764], const double dux[1764], const double duy[1764],
		const double dvx[1764], const double dvy[1764], double hextra, const
		double zc[1764], emxArray_real_T* F, emxArray_real_T* G, double
		* amax)
	{
		double hr;
		double hl;
		int i0;
		int loop_ub;
		double ul;
		double vl;
		int k;
		emxArray_int32_T* r0;
		double ur;
		int j;
		double vr;
		double dv0[3];
		double zbc;
		*amax = 0.0;
		hr = 0.0;
		hl = 0.0;
		i0 = F->size[0] * F->size[1] * F->size[2];
		F->size[0] = (int)n;
		F->size[1] = (int)n;
		F->size[2] = (int)n;
		//emxEnsureCapacity_real_T(F, i0);
		loop_ub = (int)n * (int)n * (int)n;
		for (i0 = 0; i0 < loop_ub; i0++) {
			F->data[i0] = 0.0;
		}

		i0 = G->size[0] * G->size[1] * G->size[2];
		G->size[0] = (int)n;
		G->size[1] = (int)n;
		G->size[2] = (int)n;
		//emxEnsureCapacity_real_T(G, i0);
		loop_ub = (int)n * (int)n * (int)n;
		for (i0 = 0; i0 < loop_ub; i0++) {
			G->data[i0] = 0.0;
		}

		ul = 0.0;

		//  added
		vl = 0.0;

		//  added
		//  Enforce wall boundary on left side of box (west)
		k = 1;
	//	emxInit_int32_T(&r0, 1);
		while (k - 1 <= (int)((n - 1.0) + -1.0) - 1) {
			hr = UP[1 + 42 * k] - zc[1 + 42 * k];
			if (hr < 0.0) {
				hr = 0.0;
			}

			if (hr == 0.0) {
				ur = 0.0;
				vr = 0.0;
			}
			else {
				ur = UP[1765 + 42 * k] / (hr + hextra);
				vr = UP[3529 + 42 * k] / (hr + hextra);
			}

			//////fsolver(hr, hr, -ur, ur, vr, vr, 0.0, 1.0, hextra, dv0, &zbc);
			loop_ub = F->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//	emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				F->data[F->size[0] * k + F->size[0] * F->size[1] * r0->data[i0]] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}

			k++;
		}

		//  Enforce wall boundary on second rows along the X-axis
		for (k = 1; k - 1 < (int)((n - 1.0) + -1.0); k++) {
			hl = UP[1 + 42 * k] - zc[1 + 42 * k];
			if (hl < 0.0) {
				hl = 0.0;
			}

			if (hl == 0.0) {
				ul = 0.0;
				vl = 0.0;
			}
			else {
				ul = UP[1764 + 42 * k] / (hl + hextra);
				vl = UP[3528 + 42 * k] / (hl + hextra);
			}

			////fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
			loop_ub = F->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				F->data[(F->size[0] * k + F->size[0] * F->size[1] * r0->data[i0]) + 1] =
					dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		// Compute the fluxe in the X-direction on the domain
		for (j = 2; j - 2 < (int)((n - 1.0) + -2.0); j++) {
			for (k = 0; k < (int)n; k++) {
				if ((zc[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1] > zc[j + 42 *
					k]) || isnan(zc[j + 42 * k])) {
					zbc = zc[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1];
				}
				else {
					zbc = zc[j + 42 * k];
				}

				if (UP[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1] - zbc > 0.0) {
					hl = UP[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1] - zbc;
				}
				else {
					if (UP[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1] - zbc <= 0.0)
					{
						hl = 0.0;
					}
				}

				if (UP[j + 42 * k] - zbc > 0.0) {
					hr = UP[j + 42 * k] - zbc;
				}
				else {
					if (UP[j + 42 * k] - zbc <= 0.0) {
						hr = 0.0;
					}
				}

				if (hl > 0.0) {
					hl += 0.5 * dwsex[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) - 1];
					if (hl < 0.0) {
						hl = 0.0;
					}
				}

				if (hr > 0.0) {
					hr -= 0.5 * dwsex[j + 42 * k];
					if (hr < 0.0) {
						hr = 0.0;
					}
				}

				if (hl == 0.0) {
					ul = 0.0;
					vl = 0.0;
				}
				else {
					ul = UP[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) + 1763] / (hl +
						hextra) + 0.5 * dux[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) -
						1];
					vl = UP[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) + 3527] / (hl +
						hextra) + 0.5 * dvx[((int)((3.0 + (double)(j - 2)) - 1.0) + 42 * k) -
						1];
				}

				if (hr == 0.0) {
					ur = 0.0;
					vr = 0.0;
				}
				else {
					ur = UP[1764 + (j + 42 * k)] / (hr + hextra) - 0.5 * dux[j + 42 * k];
					vr = UP[3528 + (j + 42 * k)] / (hr + hextra) - 0.5 * dvx[j + 42 * k];
				}

				////fsolver(hl, hr, ul, ur, vl, vr, 0.0, 1.0, hextra, dv0, &zbc);
				loop_ub = F->size[2];
				i0 = r0->size[0];
				r0->size[0] = loop_ub;
				//emxEnsureCapacity_int32_T(r0, i0);
				for (i0 = 0; i0 < loop_ub; i0++) {
					r0->data[i0] = i0;
				}

				loop_ub = r0->size[0];
				for (i0 = 0; i0 < loop_ub; i0++) {
					F->data[(j + F->size[0] * k) + F->size[0] * F->size[1] * r0->data[i0]] =
						dv0[i0];
				}

				if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
				}
				else {
					*amax = zbc;
				}
			}
		}

		//  From 20th to 21st and 22nd rows
		for (j = 0; j < 2; j++) {
			for (k = 0; k < (int)(n - 1.0); k++) {
				zbc = zc[(j + 42 * k) + 19];
				if (UP[(j + 42 * k) + 19] - zbc > 0.0) {
					hl = UP[(j + 42 * k) + 19] - zbc;
					if (hl < 0.0) {
						hl = 0.0;
					}

					ul = UP[(j + 42 * k) + 1783] / (hl + hextra);
					vl = UP[(j + 42 * k) + 3547] / (hl + hextra);
				}
				else {
					if (UP[(j + 42 * k) + 19] - zbc <= 0.0) {
						hl = 0.0;
						ul = 0.0;
						vl = 0.0;
					}
				}

				////fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
				loop_ub = F->size[2];
				i0 = r0->size[0];
				r0->size[0] = loop_ub;
				//emxEnsureCapacity_int32_T(r0, i0);
				for (i0 = 0; i0 < loop_ub; i0++) {
					r0->data[i0] = i0;
				}

				loop_ub = r0->size[0];
				for (i0 = 0; i0 < loop_ub; i0++) {
					F->data[((j + F->size[0] * k) + F->size[0] * F->size[1] * r0->data[i0])
						+ 20] = dv0[i0];
				}

				if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
				}
				else {
					*amax = zbc;
				}
			}
		}

		//  % For the 23rd row
		// Compute the fluxe in the X-direction on the domain
		for (k = 0; k < (int)n; k++) {
			zbc = zc[22 + 42 * k];
			if (UP[22 + 42 * k] - zbc > 0.0) {
				hl = UP[22 + 42 * k] - zbc;
				hr = UP[22 + 42 * k] - zbc;
			}
			else {
				if (UP[22 + 42 * k] - zbc <= 0.0) {
					hl = 0.0;
				}

				if (UP[22 + 42 * k] - zbc <= 0.0) {
					hr = 0.0;
				}
			}

			if (hr > 0.0) {
				hr -= 0.5 * (UP[22 + 42 * k] - UP[22 + 42 * k]);

				//  wse difference bn 25 & 26
				if (hr < 0.0) {
					hr = 0.0;
				}
			}

			if (hl == 0.0) {
				ul = 0.0;
				vl = 0.0;
			}
			else {
				ul = UP[1786 + 42 * k] / (hl + hextra);
				vl = UP[3550 + 42 * k] / (hl + hextra);
			}

			if (hr == 0.0) {
				ur = 0.0;
				vr = 0.0;
			}
			else {
				ur = UP[1786 + 42 * k] / (hr + hextra);
				vr = UP[3550 + 42 * k] / (hr + hextra);
			}

			//fsolver(hl, hr, ul, ur, vl, vr, 0.0, 1.0, hextra, dv0, &zbc);
			loop_ub = F->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				F->data[(F->size[0] * k + F->size[0] * F->size[1] * r0->data[i0]) + 22] =
					dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//  Enforce wall boundary on right side of box (East)
		for (k = 0; k < (int)n; k++) {
			hl = UP[((int)(n - 1.0) + 42 * k) - 1] - zc[((int)(n - 1.0) + 42 * k) - 1];
			if (hl < 0.0) {
				hl = 0.0;
			}

			if (hl == 0.0) {
				ul = 0.0;
				vl = 0.0;
			}
			else {
				ul = UP[((int)(n - 1.0) + 42 * k) + 1763] / (hl + hextra);
				vl = UP[((int)(n - 1.0) + 42 * k) + 3527] / (hl + hextra);
			}

			//fsolver(hl, hl, ul, -ul, vl, vl, 0.0, 1.0, hextra, dv0, &zbc);
			loop_ub = F->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				F->data[(((int)n + F->size[0] * k) + F->size[0] * F->size[1] * r0->data[i0])
					- 1] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//
		// Compute fluxes in y-direction
		//  Enforce wall boundary on bottom of box (South)
		for (j = 1; j - 1 < (int)((n - 1.0) + -1.0); j++) {
			hl = UP[42 + j] - zc[42 + j];
			if (hl < 0.0) {
				hl = 0.0;
			}

			if (hl == 0.0) {
				ul = 0.0;
				vl = 0.0;
			}
			else {
				ul = UP[1764 + j] / (hr + hextra);
				vl = UP[3528 + j] / (hr + hextra);
			}

			//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
			loop_ub = G->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				G->data[j + G->size[0] * G->size[1] * r0->data[i0]] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		// Enforce wall boundary on second column along Y
		for (j = 1; j - 1 < (int)((n - 1.0) + -1.0); j++) {
			hl = UP[42 + j] - zc[42 + j];
			if (hl < 0.0) {
				hl = 0.0;
			}

			if (hl == 0.0) {
				ul = 0.0;
				vl = 0.0;
			}
			else {
				ul = UP[1764 + j] / (hl + hextra);
				vl = UP[3528 + j] / (hl + hextra);
			}

			//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
			loop_ub = G->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				G->data[(j + G->size[0]) + G->size[0] * G->size[1] * r0->data[i0]] =
					dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//  Compute fluxes along Y
		for (k = 2; k - 2 < (int)((n - 1.0) + -2.0); k++) {
			for (j = 0; j < (int)n; j++) {
				if ((zc[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)] > zc[j + 42 *
					k]) || isnan(zc[j + 42 * k])) {
					zbc = zc[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)];
				}
				else {
					zbc = zc[j + 42 * k];
				}

				if (UP[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)] - zbc > 0.0) {
					hl = UP[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)] - zbc;
				}
				else {
					if (UP[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)] - zbc <= 0.0)
					{
						hl = 0.0;
					}
				}

				if (UP[j + 42 * k] - zbc > 0.0) {
					hr = UP[j + 42 * k] - zbc;
				}
				else {
					if (UP[j + 42 * k] - zbc <= 0.0) {
						hr = 0.0;
					}
				}

				if (hl > 0.0) {
					hl += 0.5 * dwsey[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)];
					if (hl < 0.0) {
						hl = 0.0;
					}
				}

				if (hr > 0.0) {
					hr -= 0.5 * dwsey[j + 42 * k];
					if (hr < 0.0) {
						hr = 0.0;
					}
				}

				if (hl == 0.0) {
					ul = 0.0;
					vl = 0.0;
				}
				else {
					ul = (UP[1764 + (j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1))] +
						0.5 * duy[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)]) /
						(hl + hextra);
					vl = (UP[3528 + (j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1))] +
						0.5 * dvy[j + 42 * ((int)((3.0 + (double)(k - 2)) - 1.0) - 1)]) /
						(hl + hextra);
				}

				if (hr == 0.0) {
					ur = 0.0;
					vr = 0.0;
				}
				else {
					ur = UP[1764 + (j + 42 * k)] / (hr + hextra) - 0.5 * duy[j + 42 * k];
					vr = UP[3528 + (j + 42 * k)] / (hr + hextra) - 0.5 * dvy[j + 42 * k];
				}

				//fsolver(hl, hr, ul, ur, vl, vr, 1.0, 0.0, hextra, dv0, &zbc);
				loop_ub = G->size[2];
				i0 = r0->size[0];
				r0->size[0] = loop_ub;
				//emxEnsureCapacity_int32_T(r0, i0);
				for (i0 = 0; i0 < loop_ub; i0++) {
					r0->data[i0] = i0;
				}

				loop_ub = r0->size[0];
				for (i0 = 0; i0 < loop_ub; i0++) {
					G->data[(j + G->size[0] * k) + G->size[0] * G->size[1] * r0->data[i0]] =
						dv0[i0];
				}

				if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
				}
				else {
					*amax = zbc;
				}
			}
		}

		//  for the 23rd rows /downstream of the dam
		for (k = 0; k < 22; k++) {
			for (j = 0; j < 2; j++) {
				hl = UP[(j + 42 * k) + 20] - zc[(j + 42 * k) + 20];
				if (hl < 0.0) {
					hl = 0.0;
				}

				ul = UP[(j + 42 * k) + 1784] / (hl + hextra);
				vl = UP[(j + 42 * k) + 3548] / (hl + hextra);
				//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
				loop_ub = G->size[2];
				i0 = r0->size[0];
				r0->size[0] = loop_ub;
				//emxEnsureCapacity_int32_T(r0, i0);
				for (i0 = 0; i0 < loop_ub; i0++) {
					r0->data[i0] = i0;
				}

				loop_ub = r0->size[0];
				for (i0 = 0; i0 < loop_ub; i0++) {
					G->data[((j + G->size[0] * (k + 1)) + G->size[0] * G->size[1] * r0->
						data[i0]) + 20] = dv0[i0];
				}

				if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
				}
				else {
					*amax = zbc;
				}
			}
		}

		//  At the 24th cell at the corner of the dam opening
		for (j = 0; j < 2; j++) {
			zbc = zc[j + 986];
			if (UP[j + 986] - zbc > 0.0) {
				hl = UP[j + 986] - zbc;
				if (hl < 0.0) {
					hl = 0.0;
				}
			}
			else {
				if (UP[j + 986] - zbc <= 0.0) {
					hl = 0.0;
					ul = 0.0;
					vl = 0.0;
				}
			}

			//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
			loop_ub = G->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				G->data[((j + G->size[0] * 23) + G->size[0] * G->size[1] * r0->data[i0]) +
					20] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//  At the 35th cell at the left wing of the dam opening
		for (j = 0; j < 2; j++) {
			zbc = zc[j + 1406];
			if (UP[j + 1406] - zbc > 0.0) {
				hl = UP[j + 1406] - zbc;
				if (hl < 0.0) {
					hl = 0.0;
				}

				ul = UP[j + 3170] / (hl + hextra);
				vl = UP[j + 4934] / (hl + hextra);
			}
			else {
				if (UP[j + 1406] - zbc <= 0.0) {
					hl = 0.0;
					ul = 0.0;
					vl = 0.0;
				}
			}

			//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
			loop_ub = G->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				G->data[((j + G->size[0] * 34) + G->size[0] * G->size[1] * r0->data[i0]) +
					20] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//  for the left part of the dam
		for (k = 0; k < (int)((n - 1.0) + -35.0); k++) {
			for (j = 0; j < 2; j++) {
				hl = UP[(j + 42 * ((int)((36.0 + (double)k) - 1.0) - 1)) + 20] - zc[(j +
					42 * ((int)((36.0 + (double)k) - 1.0) - 1)) + 20];
				if (hl < 0.0) {
					hl = 0.0;
				}

				ul = UP[(j + 42 * ((int)((36.0 + (double)k) - 1.0) - 1)) + 1784] / (hl +
					hextra);
				vl = UP[(j + 42 * ((int)((36.0 + (double)k) - 1.0) - 1)) + 3548] / (hl +
					hextra);
				//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
				loop_ub = G->size[2];
				i0 = r0->size[0];
				r0->size[0] = loop_ub;
				//emxEnsureCapacity_int32_T(r0, i0);
				for (i0 = 0; i0 < loop_ub; i0++) {
					r0->data[i0] = i0;
				}

				loop_ub = r0->size[0];
				for (i0 = 0; i0 < loop_ub; i0++) {
					G->data[((j + G->size[0] * (k + 35)) + G->size[0] * G->size[1] *
						r0->data[i0]) + 20] = dv0[i0];
				}

				if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
				}
				else {
					*amax = zbc;
				}
			}
		}

		// Enforce wall boundary on top of box (north)
		for (j = 0; j < (int)n; j++) {
			hl = UP[j + 42 * ((int)(n - 1.0) - 1)] - zc[j + 42 * ((int)(n - 1.0) - 1)];
			if (hl < 0.0) {
				hl = 0.0;
			}

			ul = UP[1764 + (j + 42 * ((int)(n - 1.0) - 1))] / (hl + hextra);
			vl = UP[3528 + (j + 42 * ((int)(n - 1.0) - 1))] / (hl + hextra);
			//fsolver(hl, hl, ul, ul, vl, -vl, 1.0, 0.0, hextra, dv0, &zbc);
			loop_ub = G->size[2];
			i0 = r0->size[0];
			r0->size[0] = loop_ub;
			//emxEnsureCapacity_int32_T(r0, i0);
			for (i0 = 0; i0 < loop_ub; i0++) {
				r0->data[i0] = i0;
			}

			loop_ub = r0->size[0];
			for (i0 = 0; i0 < loop_ub; i0++) {
				G->data[(j + G->size[0] * ((int)n - 1)) + G->size[0] * G->size[1] *
					r0->data[i0]] = dv0[i0];
			}

			if ((zbc < *amax) || (isnan(zbc) && (!isnan(*amax)))) {
			}
			else {
				*amax = zbc;
			}
		}

		//emxFree_int32_T(&r0);
	}

	//
	// File trailer for fluxes.cpp
	//
	// [EOF]
	//
};