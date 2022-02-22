#!/bin/bash
#

#BATCH --job-name=scdali
#SBATCH --time=12:00:00
#SBATCH --mem=25g
#SBATCH -N 1
#SBATCH -n 23
#SBATCH --mail-user=wancen@live.unc.edu
#SBATCH --mail-type=ALL

module load r/4.0.3
Rscript comparison.R
