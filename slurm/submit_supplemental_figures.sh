#!/bin/bash
#SBATCH --job-name=hrf_suppfig
#SBATCH --time=00:15:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=logs/suppfig_%j.out
#SBATCH --error=logs/suppfig_%j.err

set -euo pipefail
cd /users/PUOM0008/crsfaaron/HRF
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
Rscript --vanilla scripts/supplemental_figures.R
