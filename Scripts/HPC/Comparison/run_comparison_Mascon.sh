#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=15
#SBATCH --cpus-per-task=2
#SBATCH --account=ucb387_asc1
#SBATCH --time=03:30:00
#SBATCH --output=SlurmFiles/Comparison/comparison-Mascon-%A-%a.out
#--dependency=afterok:3964878_0
#--array=12-19

module purge

echo "== Load Anaconda =="

module load slurm/alpine
module load anaconda

echo "== Activate Env =="

conda activate research

echo "== Run data generation =="

srun python /projects/joma5012/GravNN/Scripts/Gen-III/Comparison/train_HPC.py $SLURM_ARRAY_TASK_ID

wait # Necessary to wait for all processes to finish
echo "== End of Job =="
exit 0
