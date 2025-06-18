#!/bin/bash
#SBATCH --job-name=gray_scott_lsrk
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --array=0-15

NAMES=("LSRK" "ERK")
METHODS=("ERK" "ERK")
LOW_STORAGE=("true" "false")
TOLS=("1e-2" "1e-3" "1e-4" "1e-5" "1e-6" "1e-7" "1e-8" "1e-9")

TOL_INDEX=$((SLURM_ARRAY_TASK_ID / ${#METHODS[@]}))
METHOD_INDEX=$((SLURM_ARRAY_TASK_ID % ${#METHODS[@]}))

mkdir -p data

../builddir/gray_scott \
  grid_pts_1d 1024 \
  method ${METHODS[$METHOD_INDEX]} \
  low_storage ${LOW_STORAGE[$METHOD_INDEX]} \
  rel_tol ${TOLS[$TOL_INDEX]} \
  abs_tol 1e-13 \
  threads $SLURM_CPUS_PER_TASK \
  out_file data/${NAMES[$METHOD_INDEX]}_${TOLS[$TOL_INDEX]}.txt
