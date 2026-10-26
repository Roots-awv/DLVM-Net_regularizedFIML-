#!/bin/bash

#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=150gb

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

module load math/R/4.3.3
R CMD BATCH --vanilla SAMdata12250-df17var7l1-FIML-lasso6.R
