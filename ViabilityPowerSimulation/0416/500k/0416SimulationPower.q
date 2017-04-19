#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH --nodes=32
#SBATCH -c 8
#SBATCH --mem=8G
#SBATCH --job-name=416Sim

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0416SimulationPower.R 