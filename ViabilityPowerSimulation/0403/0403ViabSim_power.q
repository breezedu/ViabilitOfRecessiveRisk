#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH --nodes=16
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH --job-name=ViaSimPower

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0403ViabSim_power.R 