#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<vector>
#include "mpi.h"


#undef VERSION_STRING
#define VERSION_STRING "MPIMatrix 0.1"

#define DIM 2
#define GHOST_VAL 0.0
#define MPI_VALUE_T MPI_FLOAT

typedef std::initializer_list<int> init_list;
typedef std::pair<int, int> dims_t; // first=rows, second=cols
typedef double value_t;

dims_t get_local_dims(int globalrows, int globalcols, int rank, int size);

/*
  supports partitioning into a 1D or 2D Cartesian grid
*/
struct partition_data_t {
	int  size, rows, cols, sub_rows, sub_cols;
	int  cart_dims[2];
	std::vector<dims_t > part_dims;

	partition_data_t() : size(0), rows(0), cols(0) {}

	partition_data_t(size_t s, size_t r, size_t c) :
		size(s), rows(r), cols(c)
	{
		part_dims.assign(size, dims_t(0, 0));

		// block row partitioning
		for (int p = 0; p < size; ++p) {
			part_dims[p] = get_local_dims(rows, cols, p, size);
		}
	}

};

class MPIXMatrix
{

public:



	//---------------------------------------------------------------------------

	double* scatter_exchange(double* global, partition_data_t pd,
		const int rank, MPI_Comm cartcomm, int n)
	{
		dims_t lgrid = get_local_dims(pd.rows, pd.cols, rank, pd.size);
		int
			lrows = lgrid.first,
			lcols = lgrid.second,
			subsize = lrows * lcols;
		double* subgrid;
		subgrid = (double*)malloc(sizeof(double) * n * n);
		if (rank == 0) {

			memcpy(&subgrid[0], &global[0], subsize * sizeof(value_t));
			//fprintf(stderr, "Rank: %d Send: %d\n",rank,subsize);

			for (int p = 1; p < pd.size; p++) {
				if (p == pd.size - 1) {
					int rem = (pd.rows % pd.size);
					int local_row = lrows;
					if (rem > 0) {
						local_row += rem;
					}
					int local_subsize = local_row * lcols;
					MPI_Send(&global[(lrows - 2) * p * lcols], local_subsize, MPI_VALUE_T, p, 0, cartcomm);
					//fprintf(stderr, "Rank: %d Send: %d\n",p,local_subsize);
				}
				else {
					MPI_Send(&global[(lrows - 2) * p * lcols], subsize, MPI_VALUE_T, p, 0, cartcomm);
					//fprintf(stderr, "Rank: %d Send: %d\n",p,subsize);
				}
			}
		}
		else {
			MPI_Recv(&subgrid[0], subsize, MPI_VALUE_T, 0, 0, cartcomm, MPI_STATUS_IGNORE);
			//fprintf(stderr, "Rank: %d Received: %d\n",rank,subsize);
		}

		MPI_Barrier(cartcomm);

		return subgrid;
	}

	//---------------------------------------------------------------------------

	// --------------------------------------------------------------------------


	/*  std::ostream& operator<<(std::ostream& out, MPIXMatrix<T>& M)
	  {
		unsigned m = M.num_rows() ;
		unsigned n = M.num_cols() ;
		int off = 1;
		for (unsigned i = off; i < m - off; i++)
		{
		  for (unsigned j = off; j < n - off; j++)
		  {
			out << (T)M(i, j) ;

			if (j < (n - off) - off)
			{
			  out << " " ;
			}
		  }
		  out << std::endl ;
		}

		return out ;
	  }*/

};

// ----------------------------------- prototypes ---------------------------
dims_t get_local_dims(int globalrows, int globalcols, int rank, int size)
{
	int halo = 1;
	dims_t localgrid(halo * 2, halo * 2);

	localgrid.first += globalrows / size;
	localgrid.second += globalcols;

	int rem = (globalrows % size);

	if (rank == size - 1 && rem > 0) {
		localgrid.first += rem;
	}

	return localgrid;
}
