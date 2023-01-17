#!/bin/bash
#SBATCH --job-name=extractCMIP
#SBATCH --nodes=1
#SBATCH --array=1982-1983
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
cd $OAK/aminaly/biodiversity_cmip6/raw_data/NDVI/www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/$SLURM_ARRAY_TASK_ID

cdo ensmean *.nc /oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/processed_data/NDVI/%year%_mean.nc
cdo ensmax *.nc /oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/processed_data/NDVI/%year%_max.nc
cdo ensmin *.nc /oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/processed_data/NDVI/%year%_min.nc

