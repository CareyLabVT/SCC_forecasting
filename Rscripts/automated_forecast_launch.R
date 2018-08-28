if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"lubridate" %in% installed.packages()) install.packages("lubridate")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)
library(lubridate)

Folder <- '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/'
forecast_location <- '/Users/quinn/Dropbox/Research/SSC_forecasting/temp_forecast/' 
current_day <- '2018-08-10 00:00:00'
num_days <- 17
#--------
#   Use launch_mode <- 1 if you want to schedule the forecast to run at a certain time
#   Use launch_mode <-2 if you want to the forecast to lauch at a given time interval
#-------
launch_mode <- 2
#--------
#   if launch_mode == 1: launch_time corresponds to the time of launch
#   if launch_mode == 2: launch_time corresponds to the time between launches
#-------
launch_time <- 60*60*2.5  #"12:00:00"

source(paste0(Folder,'/','Rscripts/EnKF_GLM_wNOAAens_V2.R'))
source(paste0(Folder,'/','Rscripts/evaluate_forecast.R'))

#FIRST DAY
out <- run_forecast(
  first_day = current_day,
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 0,
  restart_file = NA,
  Folder = Folder,
  forecast_location = forecast_location
)

day_count <- 0
#ALL SUBSEQUENT DAYS
repeat {
  startTime <- Sys.time()
  if((as.POSIXct(current_day, format = "%Y-%m-%d %H:%M:%S") + days(1)  > Sys.time()) | day_count > num_days){
    #There won't be a NOAA forecast avaialable so end the forecasting
    break
  }
  #ADVANCE TO NEXT DAY
  tmp <- as.POSIXct(current_day, format = "%Y-%m-%d %H:%M:%S") + days(1)
  current_day <- paste0(strftime(tmp,format = "%Y-%m-%d",usetz = FALSE)," 00:00:00")
  out <- run_forecast(
    first_day= current_day,
    sim_name = NA, 
    hist_days = 1,
    forecast_days = 10,
    restart_file = paste0(forecast_location,'/',unlist(out)[3],'/',unlist(out)[1]),
    Folder = Folder,
    forecast_location = forecast_location
  )
  day_count <- day_count + 1
  if(launch_mode == 2){
    sleepTime <- startTime + launch_time - Sys.time()
  }else if(launch_mode == 1){
    tmp <- as.POSIXct(Sys.time(), format = "%Y-%m-%d %H:%M:%S") + days(1)
    next_launch <- paste0(strftime(tmp,format = "%Y-%m-%d",usetz = FALSE)," ",launch_time)
    sleepTime <- as.POSIXct(next_launch) - Sys.time()
  }
  sleepTime <- startTime + days(1) - Sys.time()
  if (sleepTime > 0){
    Sys.sleep(sleepTime)
  }
}