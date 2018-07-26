if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)
library(lubridate)

Folder = '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/'

## EXAMPLE LAUCHING A FORECAST
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/EnKF_GLM_wNOAAens_V2.R')
out <- run_forecast(
  first_day = '2018-07-06 00:00:00',
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 1,
  restart_file = NA,
  Folder = Folder,
  machine = 'mac'
  )

## EXAMPLE EVALUATING FORECAST AFTER TIME HAS PAST
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/evaluate_forecast.R')
evaluate_forecast(
  forecast_folder = 'forecast_2018_7_6_2018726_12_9',
  Folder = Folder,
  sim_name = '2018_7_6'
)

## EXAMPLE EVALUATING FORECAST AFTER TIME HAS PAST USING THE WHAT IS RETURNED FROM RUN_FORECAST.R
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/evaluate_forecast.R')
evaluate_forecast(
  forecast_folder = unlist(out)[3],
  Folder = Folder,
  sim_name = unlist(out)[2]
)

## EXAMPLE OF LAUCHING A FORECAST FROM A PREVIOUS STEP THROUGH THE ENKF

#INIITIAL LAUNCH
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/EnKF_GLM_wNOAAens_V2.R')
out <- run_forecast(
  first_day= '2018-07-06 00:00:00',
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 1,
  restart_file = NA,
  Folder = Folder,
  machine = 'mac'
)

#SUBSEQUENT DAYS LAUNCH
restart_file_name <- run_forecast(first_day= '2018-07-07 00:00:00',
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 8,
  restart_file = paste0(Folder,'/Forecasts/',unlist(out)[3],'/',unlist(out)[1]),
  Folder = Folder,
  machine = 'mac'
)

