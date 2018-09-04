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
start_day <- '2018-07-15 00:00:00'
forecast_start_day <- '2018-09-04 00:00:00'
num_forecast_days <- NA  #Set to NA if 

hist_days <- as.numeric(difftime(as.POSIXct(forecast_start_day, format = "%Y-%m-%d %H:%M:%S"), as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S")))


wait_time <- 60*60*2.5

push_to_git <- TRUE

source(paste0(Folder,'/','Rscripts/EnKF_GLM_wNOAAens_V2.R'))
source(paste0(Folder,'/','Rscripts/evaluate_forecast.R'))

#FIRST DAY
out <- run_forecast(
  first_day = start_day,
  sim_name = 'test_historical', 
  hist_days = hist_days-1,
  forecast_days = 0,
  restart_file = NA,
  Folder = Folder,
  forecast_location = forecast_location,
  push_to_git=push_to_git
)


#ADVANCE TO NEXT DAY
start_day <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S") + days(hist_days) - days(1)
forecast_day_count <- 1
#ALL SUBSEQUENT DAYS
repeat {
  
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
    
    tmp <-getURL(paste0('https://github.com/CareyLabVT/SCCData/raw/noaa-data/',forecast_base_name))
    tmp <- unlist(strsplit(tmp, '<'))
    if(tmp[2] == "!DOCTYPE html>\n"){
      print('Waiting for NOAA forecast')
      Sys.sleep(wait_time)
    }else{
      forecast_avialable = TRUE
    }
  }
  
  start_day <- paste0(strftime(start_day,format = "%Y-%m-%d",usetz = FALSE)," 00:00:00")
  
  out <- run_forecast(
    first_day= start_day,
    sim_name = NA, 
    hist_days = 1,
    forecast_days = 15,
    restart_file = paste0(forecast_location,'/',unlist(out)[3],'/',unlist(out)[1]),
    Folder = Folder,
    forecast_location = forecast_location,
    push_to_git=push_to_git
  )
  forecast_day_count <- forecast_day_count + 1
  
  #ADVANCE TO NEXT DAY
  start_day <- as.POSIXct(start_day, format = "%Y-%m-%d %H:%M:%S") + days(1)
  if(!is.na(num_forecast_days)){
    if(forecast_day_count > num_forecast_days){
      break
    }
  }
  
}