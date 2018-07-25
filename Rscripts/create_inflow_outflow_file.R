create_inflow_outflow_file <- function(full_time_day){

  full_time_day_past <- as.POSIXct(full_time_day) - 365*24*60*60
  
  inflow = read.csv(paste0(workingGLM,'FCR_weir_inflow_2013_2017_20180716.csv'))
  spillway = read.csv(paste0(workingGLM,'FCR_spillway_outflow_2013_2017_20180716.csv'))
  
  inflow_new =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_past),]
  inflow_new$time =  full_time_day
  
  spillway_new = spillway[which(as.POSIXct(spillway$time) %in% full_time_day_past),,]
  spillway_new$time =  full_time_day
  
  write.csv(inflow_new,file = paste(out_directory,'FCR_inflow.csv',sep=''),row.names = FALSE,quote = FALSE)
  write.csv(spillway_new,file = paste(out_directory,'FCR_spillway_outflow.csv',sep=''),row.names = FALSE,quote = FALSE)
}