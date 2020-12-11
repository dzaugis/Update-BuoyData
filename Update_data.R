library(rerddap)
library(tidyverse)
library(lubridate)
library(foreach)
library(ncdf4)
source("Update_buoy_functions.R")


# buoy (can be a vector)
# destfolder: place to deposit files
buoys <- c("B01", "E01", "F01", "I01", "M01", "N01")
destfolder <- "buoy_dataNC"

DownloadBuoyCTD(buoy = buoys, destfolder = destfolder)

files <- list.files("buoy_dataNC", full.names = TRUE)
variables <- c("temperature", "salinity", "sigma_t")

buoy_data <- loadNERACOOSnetcdf(filename = files, varname = variables)

buoy_data2 <- surfbotcorr(buoy_data)

buoy_data_daily <- makeDailyMean(buoy_data)

