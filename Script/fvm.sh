#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=6000M
#SBATCH --time=00:10:00
#SBATCH --mail-user=rmanicava42@students.tntech.edu --mail-type=ALL
#SBATCH --error=preprocesserror.err


module load gcc/6.3.0
#g++ corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp 
mpic++ -o main  corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp -lm -fopenmp
mpirun -n 2 ./main
