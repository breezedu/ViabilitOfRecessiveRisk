#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH --nodes=16
#SBATCH --job-name=ViaSimPower

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0330ViabSim_power.R 