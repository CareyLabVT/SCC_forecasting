## EXAMPLE LAUCHING A FORECAST
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/EnKF_GLM_wNOAAens_V2.R')
out <- run_forecast(
  first_day = '2018-07-06 00:00:00',
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 1,
  restart_file = NA
  )

## EXAMPLE EVALUATING FORECAST AFTER TIME HAS PAST
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/evaluate_forecast.R')
evaluate_forecast(
  forecast_folder = 'forecast_2018_7_6_2018726_12_9',
  Folder = '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/',
  sim_name = '2018_7_6'
)

## EXAMPLE EVALUATING FORECAST AFTER TIME HAS PAST USING THE WHAT IS RETURNED FROM RUN_FORECAST.R
source('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Rscripts/evaluate_forecast.R')
evaluate_forecast(
  forecast_folder = unlist(out)[3],
  Folder = '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/',
  sim_name = unlist(out)[2]
)

## EXAMPLE OF LAUCHING A FORECAST FROM A PREVIOUS STEP THROUGH THE ENKF

#INIITIAL LAUNCH
out <- run_forecast(
  first_day= '2018-07-06 00:00:00',
  sim_name = NA, 
  hist_days = 1,
  forecast_days = 1,
  restart_file = NA
)

#SUBSEQUENT DAYS LAUNCH
restart_file_name <- run_forecast(first_day= '2018-07-07 00:00:00',
             sim_name = NA, 
             hist_days = 1,
             forecast_days = 1,
             restart_file = paste0('/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/Forecasts/',unlist(out)[3],'/',unlist(out)[1])
)

