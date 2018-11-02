
extract_temp_chain <- function(fname,full_time,input_tz, output_tz){
  d <- read.csv(fname, skip =4, na.strings = 'NAN')
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  obs <- array(NA,dim=c(length(full_time),10))
  depths <- c(0.1,1,2,3,4,5,6,7,8,9)
  
  TIMESTAMP_in <- as.POSIXct(d$TIMESTAMP,origin = '1970-01-01 00:00.00 UTC',tz = input_tz)
  d$TIMESTAMP <- as.POSIXct(TIMESTAMP_in,tz = output_tz)
  full_time <- as.POSIXct(full_time,tz = output_tz)
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


extract_do_chain <- function(fname = catwalk_fname,full_time,input_tz = 'EST5EDT', output_tz = 'GMT'){
  d <- read.csv(fname, skip =3, na.strings = 'NAN')
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  obs <- array(NA,dim=c(length(full_time),3))
  depths <- c(1,5,9)
  
  TIMESTAMP_in <- as.POSIXct(d$TIMESTAMP,origin = '1970-01-01 00:00.00 UTC',tz = input_tz)
  d$TIMESTAMP <- as.POSIXct(TIMESTAMP_in,tz = output_tz)
  full_time <- as.POSIXct(full_time,tz = output_tz)
  for(i in 1:length(full_time)){
    index = which(abs(as.POSIXct(d$TIMESTAMP)-as.POSIXct(full_time[i])) == min(abs(as.POSIXct(d$TIMESTAMP) - as.POSIXct(full_time[i]))))
    index = which(as.POSIXct(d$TIMESTAMP)==as.POSIXct(full_time[i]))
    if(length(index)>0){
      obs[i,1] <- max(c(d$doobs_1[index],0.0))
      obs[i,2] <- max(c(d$doobs_5[index],0.0))
      obs[i,3] <- max(c(d$doobs_9[index],0.0))
    }
  }
  return(list(obs = obs, depths = depths))
}

extract_chla_chain <- function(fname = catwalk_fname,full_time,input_tz = 'EST5EDT', output_tz = 'GMT'){
  d <- read.csv(fname, skip =3, na.strings = 'NAN')
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  Chla_obs <- array(NA,dim=c(length(full_time)))
  BGAPC_obs <- array(NA,dim=c(length(full_time)))
  fDOM_obs <- array(NA,dim=c(length(full_time)))
  
  depths <- c(1)
  
  TIMESTAMP_in <- as.POSIXct(d$TIMESTAMP,origin = '1970-01-01 00:00.00 UTC',tz = input_tz)
  d$TIMESTAMP <- as.POSIXct(TIMESTAMP_in,tz = output_tz)
  full_time <- as.POSIXct(full_time,tz = output_tz)
  for(i in 1:length(full_time)){
    index = which(abs(as.POSIXct(d$TIMESTAMP)-as.POSIXct(full_time[i])) == min(abs(as.POSIXct(d$TIMESTAMP) - as.POSIXct(full_time[i]))))
    index = which(as.POSIXct(d$TIMESTAMP)==as.POSIXct(full_time[i]))
    if(length(index)>0){
      Chla_obs[i] <- d$Chla_1[index]
      BGAPC_obs[i] <- d$BGAPC_1[index]
      fDOM_obs[i] <- d$fDOM_QSU_1[index]
    }
  }
  return(list(Chla_obs = Chla_obs, BGAPC_obs = BGAPC_obs, fDOM_obs = fDOM_obs, depths = depths))
}