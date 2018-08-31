
extract_temp_chain <- function(fname,full_time){
  d <- read.csv(fname, skip =4, na.strings = 'NAN')
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  obs <- array(NA,dim=c(length(full_time),10))
  depths <- c(0.1,1,2,3,4,5,6,7,8,9)
  
  d$TIMESTAMP <- as.POSIXct(d$TIMESTAMP)
  full_time <- as.POSIXct(full_time)
  for(i in 1:length(full_time)){
    index = which(d$TIMESTAMP==full_time[i])
    if(length(index)>0){
      obs[i,] <- unlist(d[index,5:14])
      if(is.na(obs[i,2]) & !is.na(d[index,23])){
        obs[i,2] <- d[index,23]
      }
      if(is.na(obs[i,6]) & !is.na(d[index,17])){
        obs[i,6] <- d[index,17]
      }
      if(is.na(obs[i,10]) & !is.na(d[index,20])){ 
        obs[i,10] <- d[index,20]
      }
    }
  }
  
  return(list(obs = obs, depths = depths))
}


extract_do_chain <- function(fname = catwalk_fname,full_time){
  d <- read.csv(fname, skip =3)
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  obs <- array(NA,dim=c(length(full_time),3))
  depths <- c(1,5,9)
  
  d$TIMESTAMP <- as.POSIXct(d$TIMESTAMP)
  full_time <- as.POSIXct(full_time)
  for(i in 1:length(full_time)){
    index = which(abs(as.POSIXct(d$TIMESTAMP)-as.POSIXct(full_time[i])) == min(abs(as.POSIXct(d$TIMESTAMP) - as.POSIXct(full_time[i]))))
    index = which(as.POSIXct(d$TIMESTAMP)==as.POSIXct(full_time[i]))
    if(length(index)>0){
      obs[i,1] <- d$doobs_5[index]
      obs[i,2] <- d$doobs_5[index]
      obs[i,3] <- d$doobs_9[index]
    }
  }
  return(list(obs = obs, depths = depths))
}