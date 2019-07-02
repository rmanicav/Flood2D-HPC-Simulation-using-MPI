#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=6000M
#SBATCH --time=00:10:00
#SBATCH --mail-user=rmanicava42@students.tntech.edu --mail-type=ALL
#SBATCH --error=preprocesserror.err

g++ corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp 
./a.out
