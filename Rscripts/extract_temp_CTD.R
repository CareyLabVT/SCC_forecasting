extract_temp_CTD <- function(fname,full_time_day,depths){
  
  d <- read.csv(fname)
  d_fcr <- d[which(d$Reservoir == 'FCR' & d$Site == '50'),]
  d_fcr_day <- as.POSIXct(strftime(d_fcr$Date, format="%Y-%m-%d"))
  obs <- array(NA,dim=c(length(full_time_day),length(depths)))
  
  for(i in 1:length(full_time_day)){
    index1 = which(d_fcr_day==full_time_day[i])
    if(length(index1)>0){
      curr_day <- d_fcr[index1,]
      for(j in 1:length(depths)){
        index2 <- which.min(abs(curr_day$Depth_m -depths[j]))
        obs[i,j] <- curr_day$Temp_C[index2]
      }
    }
  }
  
  return(list(obs = obs, depths = depths))
}