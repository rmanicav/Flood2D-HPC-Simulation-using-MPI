#!/bin/bash
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --mem=6000M
#SBATCH --time=00:10:00
#SBATCH --mail-user=rmanicava42@students.tntech.edu --mail-type=ALL
#SBATCH --error=preprocesserror.err

module load cuda91/toolkit/9.1.85
nvcc -c cuda_slope.cu 
g++ corrector.cpp fluxes.cpp limiter.cpp main.cpp predictor.cpp slope.cpp solver.cpp cuda_slope.o -L/usr/local/lib64 -lcuda -lcudart 
./a.out
