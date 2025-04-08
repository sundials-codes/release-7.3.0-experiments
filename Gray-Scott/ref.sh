#!/bin/bash
#SBATCH --job-name=gray_scott_ref
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --array=0-3

MESH_SIZES=("128" "256" "512" "1024")
MESH_SIZE=${MESH_SIZES[$SLURM_ARRAY_TASK_ID]}

mkdir -p data

./build/gray_scott \
  grid_pts_1d $MESH_SIZE \
  method ERK \
  order 5 \
  rel_tol 1e-13 \
  abs_tol 1e-13 \
  threads $SLURM_CPUS_PER_TASK \
  out_file data/ref_${MESH_SIZE}.txt
