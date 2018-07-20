extract_temp_CTD <- function(fname){
  d <- read.csv(fname)
  obs <- array(NA,dim=c(length(full_time),10))
  depths <- c(0.1,1,2,3,4,5,6,7,8,9)
  for(i in 1:length(depths)){
      index = which.min(abs(d$Depth_m -depths[i]))
      obs[1,i] <- d$Temp_C[index]
  }
  return(list(obs = obs, depths = depths))
}

