#!/bin/bash
#SBATCH -A dfc13_mri
#SBATCH -p mgc-mri
#SBATCH --gres=gpu:rtx_6000:1
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=4:00:00
#SBATCH --job-name=NEID_GRASS_all_lines
#SBATCH --chdir=/storage/home/efg5335/work/GRASS-E
#SBATCH --output=/storage/home/efg5335/work/GRASS-E/NEID_GRASS_all_lines.out

echo "About to start: $SLURM_JOB_NAME"
date
echo "Job id: $SLURM_JOBID"
echo "About to change into $SLURM_SUBMIT_DIR"
cd $SLURM_SUBMIT_DIR
echo "About to start Julia"
julia NEID_GRASS_all_lines.jl
echo "Julia exited"
date
