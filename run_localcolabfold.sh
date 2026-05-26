#!/bin/bash
#SBATCH --job-name=localcolabfold
#SBATCH --partition=gpu              # GPU partition (per your docs)
#SBATCH --gres=gpu:1                 # Request 1 GPU
#SBATCH --cpus-per-task=8            # CPU cores for MSA/relax steps
#SBATCH --mem=32G                    # System RAM
#SBATCH --time=04:00:00              # Wall time, adjust to need
#SBATCH --output=localcolabfold_%j.log    # %j = job ID

# Usage: sbatch run_localcolabfold.sh <input.fasta> <output_dir>
# Example: sbatch run_localcolabfold.sh virprot1.fasta virprot1_output

INPUT="${1:?Usage: sbatch run_localcolabfold.sh <input.fasta> <output_dir>}"
OUTPUT="${2:?Usage: sbatch run_localcolabfold.sh <input.fasta> <output_dir>}"

# Activate the conda env
source ~/miniforge3/etc/profile.d/conda.sh
conda activate colabfold

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT"

# Redirect weight cache of localcolabfold to /fast storage
export XDG_CACHE_HOME=/fast/sunny


colabfold_batch \
  --num-recycle 3 \
  --amber \
  --templates \
  --use-gpu-relax \
  --num-models 5 \
  --model-order 1,2,3,4,5 \
  --random-seed 0 \
  ${INPUT} \
  ${OUTPUT}
