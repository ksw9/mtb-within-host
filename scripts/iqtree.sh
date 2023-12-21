#!/bin/bash
#SBATCH --time=10:00:00 # Walltime
#SBATCH --nodes=1 # Use 1 Node (Unless code is multi-node parallelized)
#SBATCH --ntasks=1 # We only run one R instance = 1 task
#SBATCH --cpus-per-task=8 # number of threads we want to run on
#SBATCH --account=owner-guest
#SBATCH --partition=notchpeak-guest
#SBATCH -o slurm-%j.out-%N
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@utah.edu # Your email address
#SBATCH --job-name=iqtree

module load iqtree/2.2.2.4

fasta=$1
prefix=$2
model=$3

# Run IQtree
iqtree2 -s $fasta --prefix $prefix -m $model