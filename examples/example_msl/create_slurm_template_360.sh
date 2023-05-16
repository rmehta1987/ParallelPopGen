#!/bin/bash
#SBATCH --job-name=360_samp.job
#SBATCH --output=360.out
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=pi-jjberg
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

module load cuda
bash run_lots_of_sims_template.sh samples/360_set_of_samples 170
