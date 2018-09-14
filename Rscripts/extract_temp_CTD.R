extract_temp_CTD <- function(fname,full_time_day,depths){
  
  d <- read.csv(fname)
  d_fcr <- d[which(d$Reservoir == 'FCR' & d$Site == '50'),]
  d_fcr_day <- as.POSIXct(strftime(d_fcr$Date, format="%Y-%m-%d"))
  obs <- array(NA,dim=c(length(full_time),length(depths)))
  
  for(i in 1:length(full_time_day)){
    index = which(d_fcr_day==full_time_day[i])
    if(length(index)>0){
      curr_day <- d_fcr[index,]
      for(m in 1:length(depths)){
        index <- which.min(abs(curr_day$Depth_m -depths[m]))
        obs[i,m] <- curr_day$Temp_C[index]
        print(c(i,m,obs[i,m]))
      }
    }
  }
  
  return(list(obs = obs, depths = depths))
}