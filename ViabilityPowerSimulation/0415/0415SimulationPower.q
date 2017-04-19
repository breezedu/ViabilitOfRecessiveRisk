#!/bin/sh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jeff.du@duke.edu
#SBATCH --nodes=16
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH --job-name=413Sim

export PATH=/opt/apps/R-2.15.2/bin:$PATH

R CMD BATCH 0415SimulationPower.R 