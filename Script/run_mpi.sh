# Example MPI execution script for Flood2D simulation
module load gcc/6.3.0
#g++ corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp 
mpic++ -o main  corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp -lm -fopenmp
mpirun -n 2 ./main
