#!/bin/bash
#SBATCH --job-name=aggregateNDVI
#SBATCH --nodes=1
#SBATCH --array=1981-2022
#SBATCH --error=/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/outfiles/cdo.err
#SBATCH --output=/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/outfiles/cdo.out
#SBATCH --time=5:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=120GB
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=aminaly@stanford.edu
#SBATCH -p diffenbaugh

ml devel physics netcdf cdo/2.1.1

set year=$SLURM_ARRAY_TASK_ID
cd $OAK/group_members/aminaly/biodiversity_cmip6/raw_data/NDVI/www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/$SLURM_ARRAY_TASK_ID

cdo ensmean *.nc ../../../../../../../processed_data/NDVI/${SLURM_ARRAY_TASK_ID}_mean.nc
cdo ensmax *.nc ../../../../../../../processed_data/NDVI/${SLURM_ARRAY_TASK_ID}_max.nc

