#### File to run model on historical, test performance, and return figures

## load in dependencies
library(caret)

## set working directory
ifelse(dir.exists("~/Box Sync/biodiversity_cmip6"),
       setwd("~/Box Sync/biodiversity_cmip6"),
       setwd("/oak/stanford/groups/omramom/group_members/aminaly/biodiversity_cmip6"))

## source util to get reusable functions
source(paste0(getwd(), "/analysis/util.R"))

## global variables (most changes only need to happen here)
CMIP6_MODEL_NAME = "tas_Amon_historical_CMCC-ESM2_r1i1p1f1_ann_mean_2pt5degree.nc"
ML_MODEL = "earth"
  
## allow the rest of this to run
##get the data
data <- get_input_hist(CMIP6_MODEL_NAME)

## index for training and test data
data$split <- createDataPartition(data$WDPAID, p = .75, list = FALSE)

## get model
model_mars = train(max_ndvi ~ , data=trainData, method=ML_MODEL)





