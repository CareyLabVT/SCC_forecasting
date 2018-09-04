archive_forecast <- function(workingGLM,Folder,forecast_base_name,full_time,forecast_location,push_to_git,save_file_name){
  ###ARCHIVE AND CLEAN UP FORECAST
  unlink(paste0(workingGLM,'/','FCRmet.csv'),recursive = FALSE)
  unlink(paste0(workingGLM,'/','Catwalk.csv'),recursive = FALSE)
  unlink(paste0(workingGLM,'/',forecast_base_name,'.csv'),recursive = FALSE)
  time_of_forecast <- paste0(year(Sys.time()),month(Sys.time()),day(Sys.time()),'_',hour(Sys.time()),'_',(minute(Sys.time())))
  forecast_archive_dir_name <- paste0(save_file_name,'_',time_of_forecast)
  #forecast_archive_dir_name <- paste0(sim_name,forecast_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_',time_of_forecast)
  


  if(is.na(forecast_location)){
    forecast_archive_dir <- paste0(Folder,'/','Forecasts','/',forecast_archive_dir_name)
  }else{
    forecast_archive_dir <- paste0(forecast_location,'/',forecast_archive_dir_name)
  }
  
  dir.create(forecast_archive_dir)
  files <- list.files(paste0(workingGLM))
  tmp <- file.copy(files, forecast_archive_dir)
  zip(zipfile = forecast_archive_dir,files = files)
  if(push_to_git){
    setwd(forecast_location)
    system(paste0('git add ',forecast_archive_dir_name,'.zip'))
    system('git commit -m "forecast"')
    system('git push')
  }
  return(list(forecast_archive_dir_name))
}