#!/bin/bash
#SBATCH --job-name=hrf_m11
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=logs/m11_%j.out
#SBATCH --error=logs/m11_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0

export HRF_ROOT=/users/PUOM0008/crsfaaron/HRF
export OUTPUT_DIR=/users/PUOM0008/crsfaaron/HRF/output
export FIG_DIR=/users/PUOM0008/crsfaaron/HRF/output/figures
export FIELD_DIR=/users/PUOM0008/crsfaaron/HRF/data/field_inventory

Rscript --vanilla scripts/method11_itc_anchored.R
