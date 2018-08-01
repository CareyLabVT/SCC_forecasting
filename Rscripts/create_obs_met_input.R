create_obs_met_input <- function(fname,outfile,full_time_hour_obs){
  d <- read.csv(fname, skip =3)
  d_names <- read.csv(fname, skip =1)
  names(d) <- names(d_names)
  
  d <- d[-85572,]
  
  #d$day <- day(as.POSIXct(d$TIMESTAMP))
  #d$month <- month(as.POSIXct(d$TIMESTAMP))
  #d$hour <- hour(as.POSIXct(d$TIMESTAMP))
  #d$year <- year(as.POSIXct(d$TIMESTAMP))
  
  #full_time_day <- day(as.POSIXct(full_time_hour_obs))
  #full_time_month <- month(as.POSIXct(full_time_hour_obs))
  #full_time_year<- year(as.POSIXct(full_time_hour_obs))
  #full_time_hour<- hour(as.POSIXct(full_time_hour_obs))
  
  ShortWave = rep(NA,length(full_time_hour_obs)-1)
  LongWave = rep(NA,length(full_time_hour_obs)-1)
  AirTemp = rep(NA,length(full_time_hour_obs)-1)
  RelHum =rep(NA,length(full_time_hour_obs)-1)
  WindSpeed= rep(NA,length(full_time_hour_obs)-1)
  Rain = rep(NA,length(full_time_hour_obs)-1)
  Snow = rep(NA,length(full_time_hour_obs)-1)
  
  d_time_tmp <- as.POSIXct(d$TIMESTAMP)
  full_time_tmp <- as.POSIXct(full_time_hour_obs)
  
  if(length(which(d_time_tmp==full_time_tmp[1]))>0){
  
  for(i in 1:(length(full_time_hour_obs)-1)){
    index = which(d_time_tmp==full_time_tmp[i])
    index_2 = which(d_time_tmp==full_time_tmp[i+1])
    
    #index = which(d$day==full_time_year[i] & d$day==full_time_hour[i] & d$hour==full_time_hour[i] & d$hour==full_time_hour[i])
    #index_2 = which(d$day==full_time_year[i+1] & d$day==full_time_hour[i+1] & d$hour==full_time_hour[i+1])
    if(length(index)>0 & length(index_2)>0){
    ShortWave[i] = max(mean(d$SR01Up_Avg[index:index_2]),0.0)
    LongWave[i] = mean(d$IR01UpCo_Avg[index:index_2])
    AirTemp[i] = d$AirTC_Avg[index]
    RelHum[i] = d$RH[index]
    if(is.na(RelHum[i])){
      RelHum[i] = mean(c(d$RH[index-1],d$RH[index+1]))
    }
    WindSpeed[i] = mean(d$WS_ms_Avg[index:index_2])
    Rain[i] = (sum(d$Rain_mm_Tot[index:index_2])*24)/1000
    Snow[i] = 0
    }
  }
    
  remove_hours <- which(!is.na(AirTemp))
  
  ShortWave = ShortWave[remove_hours]
  LongWave = LongWave[remove_hours]
  AirTemp = AirTemp[remove_hours]
  RelHum = RelHum[remove_hours]
  WindSpeed= WindSpeed[remove_hours]
  Rain = Rain[remove_hours]
  Snow = Snow[remove_hours]
  full_time_hour_obs = full_time_hour_obs[remove_hours]
  
  
  #GMT to EST conversion
  #force_tz(full_time_hour_obs, tzone = "EST")
  #Save in GLM Format
  full_time_hour_obs = strftime(full_time_hour_obs, format="%Y-%m-%d %H:%M")
  GLM_climate = data.frame(full_time_hour_obs,ShortWave,LongWave,AirTemp,RelHum,WindSpeed,Rain,Snow)
  n= noquote(c('time','ShortWave','LongWave','AirTemp','RelHum','WindSpeed','Rain','Snow'))
  colnames(GLM_climate) = noquote(c('time','ShortWave','LongWave','AirTemp','RelHum','WindSpeed','Rain','Snow'))
  write.csv(GLM_climate,file = outfile,row.names = FALSE,quote = FALSE)
  }
}