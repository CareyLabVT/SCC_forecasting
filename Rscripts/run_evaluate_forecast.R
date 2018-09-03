source(paste0(Folder,'/','Rscripts/evaluate_forecast.R'))
## EXAMPLE EVALUATING FORECAST AFTER TIME HAS PAST USING THE WHAT IS RETURNED FROM RUN_FORECAST.R
evaluate_forecast(
  forecast_folder = 'forecast_2018_8_30_2018831_15_42',
  Folder = '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting',
  sim_name = '2018_8_30',
  forecast_location = '/Users/quinn/Dropbox/Research/SSC_forecasting/test_forecast/'
)
