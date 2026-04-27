#!/bin/bash
#SBATCH --job-name=hrf_m12
#SBATCH --time=01:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=logs/m12_%j.out
#SBATCH --error=logs/m12_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0

export HRF_ROOT=/users/PUOM0008/crsfaaron/HRF
export OUTPUT_DIR=/users/PUOM0008/crsfaaron/HRF/output
export FIELD_DIR=/users/PUOM0008/crsfaaron/HRF/data/field_inventory
export TLS_DIR=/users/PUOM0008/crsfaaron/HRF/data/tls_2025

Rscript --vanilla scripts/method12_itc_tls.R
