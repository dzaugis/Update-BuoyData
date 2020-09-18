loadNERACOOSnetcdf<-function(filename,varname){
  
  # open the netcdf file, extract the variable, and save as a vector
  ncid<-ncdf4::nc_open(filename)
  V<-ncdf4::ncvar_get(ncid,varid= varname)
  V<-as.numeric(V)
  
  # look for any flagged values, replace with NaN
  # values >0 indicate values out of range, broken sensor, or invalid input 
  Q<-ncdf4::ncvar_get(ncid,varid=paste(varname,'_qc',sep=''))
  I<-which(Q>0)
  V[I]<-NaN
  
  # load the time values
  time<-ncdf4::ncvar_get(ncid,varid='time')
  depth <- ncdf4::ncvar_get(ncid, varid = "depth")
  
  # close the netcdf file and combine time and variable output into a dataframe
  ncdf4::nc_close(ncid)
  
  time <- as.POSIXct(time, origin = "1970-01-01T00:00:00Z")
  time <- data.frame(time)
  
  df <- data.frame(cbind(time, V, depth))
  
  # Convert decimal day to Year-month-day, calculate daily means
  df <- df %>% mutate(day = day(time)) %>% 
    mutate(month = month(time)) %>% mutate(year = year(time))
  df <- df %>% group_by(year, month, day, depth) %>% summarise(daily_mean = mean(V, na.rm = TRUE)) %>%
    mutate(Date = as.Date(paste(year, month, day, sep = "-")))
  
  # Create a tibble with all of the dates to make a full time series including missing values
  start_date <- head(as.Date(df$Date), n=1)
  end_date <- tail(as.Date(df$Date), n=1)
  ts_date <- tibble(Date = seq(as_date(start_date), as_date(end_date), by = "day"))
  df <- left_join(ts_date, df[,c("Date", "daily_mean", "depth")], by = "Date", all.x =TRUE)  
}

surfbotcorr<-function(buoydat){
  # extract the surface (1m) amd bottom (50m) temperature data from each buoy
  bottomT<-buoydat %>% filter(depth == 50)
  surfT<-buoydat %>% filter(depth == 1)
  
  # match rows with the same dates to make a table of 1m and 50m temperature data
  Table <- left_join(surfT, bottomT, by = "Date", all.x = TRUE)
  
  #create a column of ordinal days and year
  Table <- Table %>% mutate(yr = year(Date)) %>% mutate(ord = yday(Date))
  
  # Linear model: 50 meter temperature dependent on 1 meter temperature grouped by ordinal day
  lm1 <-  Table %>% group_by(ord) %>% do(broom::tidy(lm(daily_mean.y ~ daily_mean.x, data = .))) %>% 
    select(ord, term, estimate) %>% spread(key = term, value = estimate) 
  
  # combine the temperature data with the slope and intercept of the lm
  Table <- merge(Table, lm1[,c("ord", "(Intercept)", "daily_mean.x")], by = "ord", all.x = TRUE)
  Table <- Table %>% group_by(ord) %>% mutate(sum = sum(is.na(daily_mean.y)))
  end <- length(Table$Date)
  
  # using the lm coefficients, calculate the 50m temperature based on the 1m temperature for each missing 50m value
  for(i in 1:end){
    if(is.na(Table$daily_mean.y[i]) & !is.na(Table$daily_mean.x.x[i]) & Table$sum[i] <10)
    {Table$daily_mean.y[i]<-Table$'(Intercept)'[i]+Table$daily_mean.x.y[i]*Table$daily_mean.x.x[i]}
  }
  Table <- ungroup(Table)
  Table <- Table %>% rename(daily_mean = daily_mean.y)
  bottomT <- Table %>% select(Date, daily_mean, "depth" = depth.y) %>% bind_rows(surfT)
  return(bottomT)
} 

library_check<- function(libraries) {
  lapply(libraries, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}

DownloadBuoy <- function(){
  buoy <- c("B01", "E01", "F01", "I01", "M01", "N01")
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores-1)
  foreach(i = 1:length(buoy)) %dopar% {
    path<-"http://www.neracoos.org/erddap/tabledap/"
    file<-"_sbe37_all.nc?station%2Ctime%2Ctemperature%2Ctemperature_qc%2Csalinity%2Csalinity_qc%2Csigma_t%2Csigma_t_qc%2Clongitude%2Clatitude%2Cdepth"
    destfile <- paste(buoy[[i]], ".nc", sep="")
    filename <- paste0(here::here(),"/buoy_dataNC/",destfile)
    
    download.file(paste(path, buoy[[i]], file, sep=""), destfile = filename, mode="wb")
  }
}


UpdateBuoy<-function(buoy,smoothing,surfbot){
  
  # Select buoy: "B01", "E01", "F01", "I01", "M01", or "N01" named in argument 1
  buoy<- buoy

  destfile <- paste(buoy, ".nc", sep="")
  filename <- paste0(here::here(),"/buoy_dataNC/",destfile)
  

  # TEMPERATURE
  temp <- loadNERACOOSnetcdf(filename = filename, "temperature")
  
  temp <- temp %>% rename("temp" = daily_mean)
  
  # SALINITY
  salt <- loadNERACOOSnetcdf(filename = filename, "salinity")
  salt <- salt %>% rename("salt" = daily_mean)
  
  # DENSITY
  den <- loadNERACOOSnetcdf(filename = filename, "sigma_t")
  den <- den %>% rename("den" = daily_mean)
  
  # To smooth or not to smooth
  if(smoothing==TRUE){
    
    # 8-day mean loop
    rolling_8 <- rollify(mean, window = 8)
    temp <- temp %>% group_by(depth) %>% mutate(rolling_mean = rolling_8(daily_mean)) 
    salt <- salt %>% group_by(depth) %>% mutate(rolling_mean = rolling_8(daily_mean))
    den <- den %>% group_by(depth) %>% mutate(rolling_mean = rolling_8(daily_mean))
  }
  
  if(surfbot == TRUE){
    temp <- surfbotcorr(temp)
  }
  
  if(buoy %in% c("B01", "E01", "F01", "I01")){
    depth_1m <- list("temp" = select(filter(temp, depth == 1), Date, "daily_mean" = temp), 
                     "sal" = select(filter(salt, depth == 1), Date, "daily_mean" = salt),
                     "density" = select(filter(den, depth == 1), Date, "daily_mean" = den))
    
    depth_20m <- list("temp" = select(filter(temp, depth == 20), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 20), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 20), Date, "daily_mean" = den))
    
    depth_50m <- list("temp" = select(filter(temp, depth == 50), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 50), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 50), Date, "daily_mean" = den)) 
    B <- list(depth_1m=depth_1m, depth_20m=depth_20m, depth_50m=depth_50m)}
  
  if(buoy == "M01"){
    depth_1m <- list("temp" = select(filter(temp, depth == 1), Date, "daily_mean" = temp), 
                     "sal" = select(filter(salt, depth == 1), Date, "daily_mean" = salt),
                     "density" = select(filter(den, depth == 1), Date, "daily_mean" = den))
    
    depth_20m <- list("temp" = select(filter(temp, depth == 20), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 20), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 20), Date, "daily_mean" = den))
    
    depth_50m <- list("temp" = select(filter(temp, depth == 50), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 50), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 50), Date, "daily_mean" = den))
    
    depth_100m <- list("temp" = select(filter(temp, depth == 100), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 100), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 100), Date, "daily_mean" = den))
    
    depth_150m <- list("temp" = select(filter(temp, depth == 150), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 150), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 150), Date, "daily_mean" = den))
    
    depth_200m <- list("temp" = select(filter(temp, depth == 200), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 200), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 200), Date, "daily_mean" = den)) 
    
    depth_250m <- list("temp" = select(filter(temp, depth == 250), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 250), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 250), Date, "daily_mean" = den))
    
    B <- list(depth_1m=depth_1m, depth_20m=depth_20m, depth_50m=depth_50m, depth_100m=depth_100m, depth_150m=depth_150m, depth_200m=depth_200m, depth_250m=depth_250m)
  }
  
  
  if(buoy == "N01"){
    depth_1m <- list("temp" = select(filter(temp, depth == 1), Date, "daily_mean" = temp), 
                     "sal" = select(filter(salt, depth == 1), Date, "daily_mean" = salt),
                     "density" = select(filter(den, depth == 1), Date, "daily_mean" = den))
    
    depth_20m <- list("temp" = select(filter(temp, depth == 20), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 20), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 20), Date, "daily_mean" = den))
    
    depth_50m <- list("temp" = select(filter(temp, depth == 50), Date, "daily_mean" = temp), 
                      "sal" = select(filter(salt, depth == 50), Date, "daily_mean" = salt),
                      "density" = select(filter(den, depth == 50), Date, "daily_mean" = den))
    
    depth_100m <- list("temp" = select(filter(temp, depth == 100), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 100), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 100), Date, "daily_mean" = den))
    
    depth_150m <- list("temp" = select(filter(temp, depth == 150), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 150), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 150), Date, "daily_mean" = den))
    
    depth_180m <- list("temp" = select(filter(temp, depth == 180), Date, "daily_mean" = temp), 
                       "sal" = select(filter(salt, depth == 180), Date, "daily_mean" = salt),
                       "density" = select(filter(den, depth == 180), Date, "daily_mean" = den)) 
    B <- list(depth_1m=depth_1m, depth_20m=depth_20m, depth_50m=depth_50m, depth_100m=depth_100m, depth_150m=depth_150m, depth_180m=depth_180m)
    
  }
  
  return(B)
  
}



update_buoy_cron <- function(){
  buoy <- c("B01", "E01", "F01", "I01", "M01", "N01")
  numCores <- parallel::detectCores()
  doParallel::registerDoParallel(numCores-1)
  Buoys <- foreach(i = 1:length(buoy)) %dopar% {
    BB<-UpdateBuoy(buoy[[i]],smoothing=FALSE,surfbot=FALSE)
  }
  Buoy_names <-c("Buoy_B", "Buoy_E", "Buoy_F", "Buoy_I", "Buoy_M", "Buoy_N")
  names(Buoys) <- Buoy_names
  rm(Buoy_names, buoy, numCores)
  save.image(file = "Buoys.RData")
  Buoys <- reshape2::melt(Buoys, id.vars = "Date") %>% 
    rename("Type" = variable, "Variable" = L3, "Depth" = L2, "name" = L1) %>% 
    mutate(Depth = as.double(parse_number(Depth)), 
           name = str_replace_all(name, "_", " "), 
           Type = paste("raw"),
           Date = as.Date(Date),
           dayz = day(Date),
           mon = month(Date),
           yr = year(Date)) %>% 
    group_by(Variable, Depth, name, dayz, mon) %>%
    mutate(clim = mean(value, na.rm = TRUE)) %>% ungroup() %>%
    mutate(Anomaly = value - clim) %>% dplyr::select(name, Date, Depth, Variable, Anomaly, "raw" = value) %>% 
    pivot_longer(., cols = c(Anomaly, raw), names_to = "Type", values_to = "Values")
  
  data.table::fwrite(Buoys, "Buoy_data.csv")
}


