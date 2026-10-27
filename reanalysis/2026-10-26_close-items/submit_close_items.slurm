#!/bin/bash
#SBATCH --job-name=howland_close
#SBATCH --account=PUOM0008
#SBATCH --time=01:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/HRF_v2maps/close_%j.out
cd /fs/scratch/PUOM0008/crsfaaron/HRF_v2maps
python3 -c "import rasterio, shapely, scipy, pandas" || exit 1
python3 howland_close_items.py . csi_check/csi_trajectories.csv close_out_20261026 > close_out_20261026.log 2>&1
echo DONE $?

