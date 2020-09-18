#!/bin/bash
#SBATCH --job-name AttractorSimulation
#SBATCH --time 2:0:0
#SBATCH --cpus 2
#SBATCH --mem 1GB
#SBATCH --mail-type NONE
#SBATCH --output /home/mpib/kamp/LNDG/AttractorModel/logs/slurm-%j.out

cd $HOME/LNDG/AttractorModel/ppattractor
module load conda 
conda activate py2attractor

ipython colored_noise.py

