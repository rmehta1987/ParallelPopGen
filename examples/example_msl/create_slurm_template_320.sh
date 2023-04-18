#!/bin/bash
#SBATCH --job-name=320_samp.job
#SBATCH --output=320.out
#SBATCH --time=12:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --account=pi-jjberg
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2

module load cuda
bash run_lots_of_sims_template.sh samples/320_set_of_samples 170
