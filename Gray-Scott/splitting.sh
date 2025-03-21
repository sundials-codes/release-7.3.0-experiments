#!/bin/bash
#SBATCH --job-name=gray_scott_splitting
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --array=0-39

ORDERS=("1" "2" "3" "4" "6")

DT_INDEX=$((SLURM_ARRAY_TASK_ID / ${#ORDERS[@]}))
ORDER_INDEX=$((SLURM_ARRAY_TASK_ID % ${#ORDERS[@]}))

DT=$(echo "2^(-$DT_INDEX)" | bc -l)
ORDER=${ORDERS[$ORDER_INDEX]}

mkdir -p data

./build/gray_scott \
  grid_pts_1d 128 \
  method Splitting \
  order $ORDER \
  dt $DT \
  threads $SLURM_CPUS_PER_TASK \
  out_file data/splitting_${ORDER}_${DT}.txt
