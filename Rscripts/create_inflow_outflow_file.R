create_inflow_outflow_file <- function(full_time_day,workingGLM){

  full_time_day_2017 <- as.POSIXct(full_time_day) - 365*24*60*60
  full_time_day_2016 <- as.POSIXct(full_time_day) - 2*365*24*60*60
  full_time_day_2015 <- as.POSIXct(full_time_day) - 3*365*24*60*60
  full_time_day_2014 <- as.POSIXct(full_time_day) - 4*365*24*60*60
  full_time_day_2015 <- as.POSIXct(full_time_day) - 5*365*24*60*60
  
  inflow = read.csv(paste0(workingGLM,'FCR_weir_inflow_2013_2017_20180716.csv'))
  spillway = read.csv(paste0(workingGLM,'FCR_spillway_outflow_2013_2017_20180716.csv'))
  
  
  inflow_2017 =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_2017),]
  inflow_2016 =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_2016),]
  inflow_2015 =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_2015),]
  inflow_2014 =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_2014),]
  inflow_2013 =  inflow[which(as.POSIXct(inflow$time) %in% full_time_day_2015),]
  inflow_new = inflow_2017
  for(i in 2:8){
  inflow_new[,i] = rowMeans(cbind(inflow_2017[,i],inflow_2016[,i],inflow_2015[,i],inflow_2014[,i],inflow_2014[,i],inflow_2013[,i]))
}
  inflow_new$time =  full_time_day
  
  spillway_2017 =  spillway[which(as.POSIXct(spillway$time) %in% full_time_day_2017),]
  spillway_2016 =  spillway[which(as.POSIXct(spillway$time) %in% full_time_day_2016),]
  spillway_2015 =  spillway[which(as.POSIXct(spillway$time) %in% full_time_day_2015),]
  spillway_2014 =  spillway[which(as.POSIXct(spillway$time) %in% full_time_day_2014),]
  spillway_2013 =  spillway[which(as.POSIXct(spillway$time) %in% full_time_day_2015),]
  spillway_new = spillway_2017
  spillway_new[,2] = rowMeans(cbind(spillway_2017[,2],spillway_2016[,2],spillway_2015[,2],spillway_2014[,2],spillway_2014[,2],spillway_2013[,2]))
  spillway_new$time =  full_time_day
  
  write.csv(inflow_new,file = paste(workingGLM,'FCR_inflow.csv',sep=''),row.names = FALSE,quote = FALSE)
  write.csv(spillway_new,file = paste(workingGLM,'FCR_spillway_outflow.csv',sep=''),row.names = FALSE,quote = FALSE)
}