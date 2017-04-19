#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH --nodes=8
#SBATCH --job-name=ViaSim64

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0327ParaViabilitySimulation_updatedPai2gRho.R 