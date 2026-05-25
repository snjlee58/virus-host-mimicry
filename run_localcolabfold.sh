#!/bin/bash
#SBATCH --job-name=localcolabfold
#SBATCH --partition=gpu              # GPU partition (per your docs)
#SBATCH --gres=gpu:1                 # Request 1 GPU
#SBATCH --cpus-per-task=8            # CPU cores for MSA/relax steps
#SBATCH --mem=32G                    # System RAM
#SBATCH --time=04:00:00              # Wall time, adjust to need
#SBATCH --output=localcolabfold_%j.log    # %j = job ID

# Activate the conda env
source ~/miniforge3/etc/profile.d/conda.sh
conda activate colabfold

# Optional: redirect weight cache off /home if quota is tight
export XDG_CACHE_HOME=/fast/sunny

INPUT="virprot1.fasta"
OUTPUT="virprot1_output"

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