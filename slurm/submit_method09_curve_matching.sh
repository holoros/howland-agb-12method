#!/bin/bash
#SBATCH --job-name=hrf_curve_match
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=logs/curve_match_%j.out
#SBATCH --error=logs/curve_match_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF

module load gcc/12.3.0
module load gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0

export HRF_ROOT=/users/PUOM0008/crsfaaron/HRF
export OUTPUT_DIR=/users/PUOM0008/crsfaaron/HRF/output
export FIG_DIR=/users/PUOM0008/crsfaaron/HRF/output/figures
export DATA_DIR=/users/PUOM0008/crsfaaron/HRF/data

Rscript --vanilla scripts/11_howland_curve_matching.R
