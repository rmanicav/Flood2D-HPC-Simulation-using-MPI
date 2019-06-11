#include<string>
#include<stdlib.h>
// Type Definitions
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T

struct emxArray_int32_T
{
	int* data;
	int* size;
	int allocatedSize;
	int numDimensions;
	bool canFreeData;
};

#endif                                 //struct_emxArray_int32_T

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
	double* data;
	int* size;
	int allocatedSize;
	int numDimensions;
	bool canFreeData;
};

#endif                                 //struct_emxArray_real_T

class fluxes_emxAPI
{
	// Function Definitions

//
// Arguments    : emxArray_int32_T *emxArray
//                int oldNumel
// Return Type  : void
//
	void emxEnsureCapacity_int32_T(emxArray_int32_T* emxArray, int oldNumel)
	{
		int newNumel;
		int i;
		void* newData;
		if (oldNumel < 0) {
			oldNumel = 0;
		}

		newNumel = 1;
		for (i = 0; i < emxArray->numDimensions; i++) {
			newNumel *= emxArray->size[i];
		}

		if (newNumel > emxArray->allocatedSize) {
			i = emxArray->allocatedSize;
			if (i < 16) {
				i = 16;
			}

			while (i < newNumel) {
				if (i > 1073741823) {
					i = INT_MAX;
				}
				else {
					i <<= 1;
				}
			}

			newData = calloc((unsigned int)i, sizeof(int));
			if (emxArray->data != NULL) {
				memcpy(newData, (void*)emxArray->data, sizeof(int) * oldNumel);
				if (emxArray->canFreeData) {
					free((void*)emxArray->data);
				}
			}

			emxArray->data = (int*)newData;
			emxArray->allocatedSize = i;
			emxArray->canFreeData = true;
		}
	}

	//
	// Arguments    : emxArray_real_T *emxArray
	//                int oldNumel
	// Return Type  : void
	//
	void emxEnsureCapacity_real_T(emxArray_real_T* emxArray, int oldNumel)
	{
		int newNumel;
		int i;
		void* newData;
		if (oldNumel < 0) {
			oldNumel = 0;
		}

		newNumel = 1;
		for (i = 0; i < emxArray->numDimensions; i++) {
			newNumel *= emxArray->size[i];
		}

		if (newNumel > emxArray->allocatedSize) {
			i = emxArray->allocatedSize;
			if (i < 16) {
				i = 16;
			}

			while (i < newNumel) {
				if (i > 1073741823) {
					i = INT_MAX;
				}
				else {
					i <<= 1;
				}
			}

			newData = calloc((unsigned int)i, sizeof(double));
			if (emxArray->data != NULL) {
				memcpy(newData, (void*)emxArray->data, sizeof(double) * oldNumel);
				if (emxArray->canFreeData) {
					free((void*)emxArray->data);
				}
			}

			emxArray->data = (double*)newData;
			emxArray->allocatedSize = i;
			emxArray->canFreeData = true;
		}
	}

	//
	// Arguments    : emxArray_int32_T **pEmxArray
	// Return Type  : void
	//
	void emxFree_int32_T(emxArray_int32_T** pEmxArray)
	{
		if (*pEmxArray != (emxArray_int32_T*)NULL) {
			if (((*pEmxArray)->data != (int*)NULL) && (*pEmxArray)->canFreeData) {
				free((void*)(*pEmxArray)->data);
			}

			free((void*)(*pEmxArray)->size);
			free((void*)* pEmxArray);
			*pEmxArray = (emxArray_int32_T*)NULL;
		}
	}

	//
	// Arguments    : emxArray_real_T **pEmxArray
	// Return Type  : void
	//
	void emxFree_real_T(emxArray_real_T** pEmxArray)
	{
		if (*pEmxArray != (emxArray_real_T*)NULL) {
			if (((*pEmxArray)->data != (double*)NULL) && (*pEmxArray)->canFreeData) {
				free((void*)(*pEmxArray)->data);
			}

			free((void*)(*pEmxArray)->size);
			free((void*)* pEmxArray);
			*pEmxArray = (emxArray_real_T*)NULL;
		}
	}

	//
	// Arguments    : emxArray_int32_T **pEmxArray
	//                int b_numDimensions
	// Return Type  : void
	//
	void emxInit_int32_T(emxArray_int32_T** pEmxArray, int b_numDimensions)
	{
		emxArray_int32_T* emxArray;
		int i;
		*pEmxArray = (emxArray_int32_T*)malloc(sizeof(emxArray_int32_T));
		emxArray = *pEmxArray;
		emxArray->data = (int*)NULL;
		emxArray->numDimensions = b_numDimensions;
		emxArray->size = (int*)malloc(sizeof(int) * b_numDimensions);
		emxArray->allocatedSize = 0;
		emxArray->canFreeData = true;
		for (i = 0; i < b_numDimensions; i++) {
			emxArray->size[i] = 0;
		}
	}

	//
	// Arguments    : emxArray_real_T **pEmxArray
	//                int b_numDimensions
	// Return Type  : void
	//
	void emxInit_real_T(emxArray_real_T** pEmxArray, int b_numDimensions)
	{
		emxArray_real_T* emxArray;
		int i;
		*pEmxArray = (emxArray_real_T*)malloc(sizeof(emxArray_real_T));
		emxArray = *pEmxArray;
		emxArray->data = (double*)NULL;
		emxArray->numDimensions = b_numDimensions;
		emxArray->size = (int*)malloc(sizeof(int) * b_numDimensions);
		emxArray->allocatedSize = 0;
		emxArray->canFreeData = true;
		for (i = 0; i < b_numDimensions; i++) {
			emxArray->size[i] = 0;
		}
	}

	//
	// File trailer for fluxes_emxutil.cpp
	//
	// [EOF]
	//

	// Function Definitions

	//
	// Arguments    : int b_numDimensions
	//                int *b_size
	// Return Type  : emxArray_real_T *
	//
	emxArray_real_T* emxCreateND_real_T(int b_numDimensions, int* b_size)
	{
		emxArray_real_T* emx;
		int numEl;
		int i;
		emxInit_real_T(&emx, b_numDimensions);
		numEl = 1;
		for (i = 0; i < b_numDimensions; i++) {
			numEl *= b_size[i];
			emx->size[i] = b_size[i];
		}

		emx->data = (double*)calloc((unsigned int)numEl, sizeof(double));
		emx->numDimensions = b_numDimensions;
		emx->allocatedSize = numEl;
		return emx;
	}

	//
	// Arguments    : double *b_data
	//                int b_numDimensions
	//                int *b_size
	// Return Type  : emxArray_real_T *
	//
	emxArray_real_T* emxCreateWrapperND_real_T(double* b_data, int b_numDimensions,
		int* b_size)
	{
		emxArray_real_T* emx;
		int numEl;
		int i;
		emxInit_real_T(&emx, b_numDimensions);
		numEl = 1;
		for (i = 0; i < b_numDimensions; i++) {
			numEl *= b_size[i];
			emx->size[i] = b_size[i];
		}

		emx->data = b_data;
		emx->numDimensions = b_numDimensions;
		emx->allocatedSize = numEl;
		emx->canFreeData = false;
		return emx;
	}

	//
	// Arguments    : double *b_data
	//                int rows
	//                int cols
	// Return Type  : emxArray_real_T *
	//
	emxArray_real_T* emxCreateWrapper_real_T(double* b_data, int rows, int cols)
	{
		emxArray_real_T* emx;
		int b_size[2];
		int numEl;
		int i;
		b_size[0] = rows;
		b_size[1] = cols;
		emxInit_real_T(&emx, 2);
		numEl = 1;
		for (i = 0; i < 2; i++) {
			numEl *= b_size[i];
			emx->size[i] = b_size[i];
		}

		emx->data = b_data;
		emx->numDimensions = 2;
		emx->allocatedSize = numEl;
		emx->canFreeData = false;
		return emx;
	}

	//
	// Arguments    : int rows
	//                int cols
	// Return Type  : emxArray_real_T *
	//
	emxArray_real_T* emxCreate_real_T(int rows, int cols)
	{
		emxArray_real_T* emx;
		int b_size[2];
		int numEl;
		int i;
		b_size[0] = rows;
		b_size[1] = cols;
		emxInit_real_T(&emx, 2);
		numEl = 1;
		for (i = 0; i < 2; i++) {
			numEl *= b_size[i];
			emx->size[i] = b_size[i];
		}

		emx->data = (double*)calloc((unsigned int)numEl, sizeof(double));
		emx->numDimensions = 2;
		emx->allocatedSize = numEl;
		return emx;
	}

	//
	// Arguments    : emxArray_real_T *emxArray
	// Return Type  : void
	//
	void emxDestroyArray_real_T(emxArray_real_T* emxArray)
	{
		emxFree_real_T(&emxArray);
	}

	//
	// Arguments    : emxArray_real_T **pEmxArray
	//                int b_numDimensions
	// Return Type  : void
	//
	void emxInitArray_real_T(emxArray_real_T** pEmxArray, int b_numDimensions)
	{
		emxInit_real_T(pEmxArray, b_numDimensions);
	}

	//
	// File trailer for fluxes_emxAPI.cpp
	//
	// [EOF]
	//
};