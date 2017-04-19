#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH -N 8
#SBATCH --job-name=ViaSim26

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0325ParaViabilitySimulation.R 