if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"lubridate" %in% installed.packages()) install.packages("lubridate")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)
library(lubridate)

Folder <- '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/'
forecast_location <- '/Users/quinn/Dropbox/Research/SSC_forecasting/test_forecast/' 
start_day <- '2018-08-30 00:00:00'
forecast_start_day <- '2018-08-31 00:00:00'

hist_days <- as.numeric(difftime(as.POSIXct(forecast_start_day, format = "%Y-%m-%d %H:%M:%S"), as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S")))

num_days <- 2
wait_time <- 60*60*2.5

push_to_git <- FALSE

source(paste0(Folder,'/','Rscripts/EnKF_GLM_wNOAAens_V2.R'))
source(paste0(Folder,'/','Rscripts/evaluate_forecast.R'))

#FIRST DAY
out <- run_forecast(
  first_day = start_day,
  sim_name = NA, 
  hist_days = hist_days,
  forecast_days = 0,
  restart_file = NA,
  Folder = Folder,
  forecast_location = forecast_location,
  push_to_git=push_to_git
)

day_count <- 0
#ALL SUBSEQUENT DAYS
repeat {
  
  startTime <- Sys.time()
  
  #ADVANCE TO NEXT DAY
  tmp <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S") + days(1)
  start_day <- paste0(strftime(tmp,format = "%Y-%m-%d",usetz = FALSE)," 00:00:00")
  
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
    
    tmp <-getURL(paste0('https://github.com/CareyLabVT/SCCData/raw/noaa-data/',forecast_base_name))
    tmp <- unlist(strsplit(tmp, '<'))
    if(tmp[2] == "!DOCTYPE html>\n"){
      print('Waiting for NOAA forecast')
      Sys.sleep(wait_time)
    }else{
      forecast_avialable = TRUE
    }
  }
  
  out <- run_forecast(
    first_day= start_day,
    sim_name = NA, 
    hist_days = 1,
    forecast_days = 10,
    restart_file = paste0(forecast_location,'/',unlist(out)[3],'/',unlist(out)[1]),
    Folder = Folder,
    forecast_location = forecast_location,
    push_to_git=push_to_git
  )
  day_count <- day_count + 1

}