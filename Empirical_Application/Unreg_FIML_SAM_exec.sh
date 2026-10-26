#!/bin/bash

#SBATCH --partition=cpu-single
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=24:00:00
#SBATCH --mem=150gb

module load math/R/4.4.2
R CMD BATCH --vanilla SAMdata12250-df17var7l1-unregFIML.R

