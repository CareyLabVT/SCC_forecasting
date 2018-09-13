if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"lubridate" %in% installed.packages()) install.packages("lubridate")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
if (!"RCurl" %in% installed.packages()) install.packages('RCurl')
library(mvtnorm)
library(glmtools)
library(ncdf4)
library(lubridate)
library(RCurl)

sim_name <- 'FCR_betaV2'
Folder <- '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/'
forecast_location <- '/Users/quinn/Dropbox/Research/SSC_forecasting/test_forecast/' 
data_location <- '/Users/quinn/Dropbox/Research/SSC_forecasting/SCC_data/' 
start_day <- '2018-07-10 00:00:00'
forecast_start_day <- '2018-09-01 00:00:00'
spin_up_days <- 5
num_forecast_days <- NA  #Set to NA if running into future
init_restart_file <- NA
init_run <- FALSE
wait_time <- 60*60*2.5
push_to_git <- FALSE



source(paste0(Folder,'/','Rscripts/EnKF_GLM_wNOAAens_V2.R'))
source(paste0(Folder,'/','Rscripts/evaluate_forecast.R'))


if(!init_run){
  hist_days <- as.numeric(difftime(as.POSIXct(forecast_start_day, format = "%Y-%m-%d %H:%M:%S"), as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S")))
  
  #FIRST DAY
  out <- run_forecast(
    first_day = start_day,
    sim_name = sim_name, 
    hist_days = hist_days-1,
    forecast_days = 0,
    spin_up_days = spin_up_days,
    restart_file = NA,
    Folder = Folder,
    forecast_location = forecast_location,
    push_to_git=push_to_git,
    data_location = data_location
  )
  
  #ADVANCE TO NEXT DAY
  start_day <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S") + days(hist_days) - days(1)
  restart_file <- unlist(out)[1]
}else{
  start_day <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S")
  restart_file <- init_restart_file
}

forecast_day_count <- 1
#ALL SUBSEQUENT DAYS
repeat{
  
  startTime <- Sys.time()
  
  
  #LOOP TO KEEP CHECKING FOR A NOAA FORECAST
  forecast_avialable = FALSE
  while(forecast_avialable == FALSE){
    forecast_start_time <- start_day + days(1)
    if(day(forecast_start_time) < 10){
      forecast_day <- paste0('0',day(forecast_start_time))
    }else{
      forecast_day <- paste0(day(forecast_start_time))
    }
    if(month(forecast_start_time) < 10){
      forecast_month <- paste0('0',month(forecast_start_time))
    }else{
      forecast_month <- paste0(month(forecast_start_time))
    }
    forecast_base_name <- paste0(year(forecast_start_time),forecast_month,forecast_day,'gep_all_00z.csv')
    
    noaa_location <- paste0(data_location,'/','noaa-data')
    setwd(noaa_location)
    system(paste0('git pull'))
    
    
    #tmp <-getURL(paste0('https://github.com/CareyLabVT/SCCData/raw/noaa-data/',forecast_base_name))
    #tmp <- unlist(strsplit(tmp, '<'))
    #if(tmp[2] == "!DOCTYPE html>\n"){
    if(!file.exists(paste0(noaa_location,'/',forecast_base_name))){
      print('Waiting for NOAA forecast')
      Sys.sleep(wait_time)
    }else{
      forecast_avialable = TRUE
    }
  }
  
  start_day <- paste0(strftime(start_day,format = "%Y-%m-%d",usetz = FALSE)," 00:00:00")
  
  out <- run_forecast(
    first_day= start_day,
    sim_name = sim_name, 
    hist_days = 1,
    forecast_days = 15,
    spin_up_days = 0,
    restart_file = restart_file,
    Folder = Folder,
    forecast_location = forecast_location,
    push_to_git=push_to_git,
    data_location = data_location
  )
  forecast_day_count <- forecast_day_count + 1
  
  restart_file <- unlist(out)[1]
  
  #ADVANCE TO NEXT DAY
  start_day <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S") + days(1)
  if(!is.na(num_forecast_days)){
    if(forecast_day_count > num_forecast_days){
      break
    }
  }
  
}