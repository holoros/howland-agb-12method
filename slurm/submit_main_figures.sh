#!/bin/bash
#SBATCH --job-name=hrf_figs2
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=logs/figs2_%j.out
#SBATCH --error=logs/figs2_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
Rscript --vanilla scripts/regenerate_figures_v45_v2.R
