#!/bin/bash
#SBATCH --job-name=extractCMIP
#SBATCH --nodes=3
#SBATCH --array=1-50
#SBATCH --error=/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/outfiles/extractCMIP.err
#SBATCH --output=/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6/outfiles/extractCMIP.out
#SBATCH --time=48:00:00
#SBATCH --ntasks-per-node=2
#SBATCH --mem-per-cpu=95GB
#SBATCH --mail-type=END,FAIL 
#SBATCH --mail-user=aminaly@stanford.edu
#SBATCH -p diffenbaugh

ml system math devel sqlite gcc
ml physics proj geos gdal udunits curl netcdf R/4.1.2;
ml gcc/9.1.0

cd $OAK/group_members/aminaly/biodiversity_cmip6
let buffer=$SLURM_ARRAY_TASK_ID
Rscript ./create_dataset.R $buffer
