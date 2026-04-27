#!/bin/bash
#SBATCH --job-name=hrf_m10
#SBATCH --time=06:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=16
#SBATCH --account=PUOM0008
#SBATCH --output=logs/m10_%j.out
#SBATCH --error=logs/m10_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0

export HRF_ROOT=/users/PUOM0008/crsfaaron/HRF
export OUTPUT_DIR=/users/PUOM0008/crsfaaron/HRF/output
export FIELD_DIR=/users/PUOM0008/crsfaaron/HRF/data/field_inventory
export NORM_DIR=/users/PUOM0008/crsfaaron/HRF/output/als_2025_processing

Rscript --vanilla scripts/method10_decimated_kohek.R
