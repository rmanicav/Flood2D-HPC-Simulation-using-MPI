//#include <unistd.h>
#include <cmath>
//#include <cstdlib>
//#include <ctime>
#include <string>
//#include <vector>
#include <iostream>
/** ------------------------------------------------------------------------/
 
  Compare
  
  Test out compare functions

  Sept, 2016
  R. Marshall

  requires C++ 11 and the following classes:

    IOMatrix
    |-FlexMatrix
    |---BaseMatrix
    
    stringutils


  compile:
    g++ -I../include --std=c++0x -o 
 
 
  usage: 
    ./

------------------------------------------------------------------------- */
#include "IOMatrix.h"

#include "stringutils.h"

typedef float value_t ;
typedef IOMatrix<value_t> domain_t;

// default values
#define RNROWS 6          // grid rows to make
#define RNCOLS 6          // grid cols to make
#define NUM_PARTITIONS 4  // number of subgrids to be joined


using namespace std;


/** ------------------------------------------------------------------------/
  read an int from command line args

  nidx    the index of argv where the arg is located
  dfault  default value if bad data

------------------------------------------------------------------------- */
int get_int_arg(int nidx, int dfault, int argc, char* argv[])
{
  return (
    (argc > nidx) && 
    (stringutils::is_numeric(std::string(argv[nidx]))) ? 
      atoi(argv[nidx]) : dfault
  ) ;

} // ------------------------------------------------------------------------

/** ------------------------------------------------------------------------/
  read a float from command line args

  nidx    the index of argv where the arg is located
  dfault  default value if bad data

------------------------------------------------------------------------- */
float get_float_arg(int nidx, float dfault, int argc, char* argv[])
{
  return (
    (argc > nidx) && 
    (stringutils::is_numeric(std::string(argv[nidx]))) ? 
      atof(argv[nidx]) : dfault
  ) ;

} // ------------------------------------------------------------------------






/**

  compute a[x] - b[x] for every x
*/
IOMatrix<value_t> ddiff(IOMatrix<value_t>* a, IOMatrix<value_t>* b, value_t eps) 
{
  if (a->num_rows() != b->num_rows() || a->num_cols() != b->num_cols())
  {
    throw idmex;
  }

  domain_t D(a->num_rows(), a->num_cols());

  for (unsigned i=0; i<a->num_rows(); ++i )
    for(unsigned j=0; j<a->num_cols(); ++j)
    {
      value_t val = (a->get_val(i,j) - b->get_val(i,j));
      D.set_val( i, j, (abs(val) < eps) ? 0 : val);
    }



    return D;
}

//---------------------------------------------------------------------------
//              main
//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
  // get command line args
  unsigned 
    num_rows = get_int_arg(1, RNROWS, argc, argv),
    num_cols = get_int_arg(2, RNCOLS, argc, argv);

    std::string af, bf;

    if (argc > 4) {
      af = std::string(argv[3]);
      bf = std::string(argv[4]);

    }
    else{
      cerr << "usage ./compare nrows ncols fileA fileB" << endl;
      exit(1);
    }


    domain_t A(num_rows, num_cols);
    domain_t B(num_rows, num_cols);  

    A.load_from_file(num_rows, num_cols, af);  
    B.load_from_file(num_rows, num_cols, bf);  

    domain_t D = ddiff(&A, &B, 0.1);

    cout << D << endl;

  return 0;
}

