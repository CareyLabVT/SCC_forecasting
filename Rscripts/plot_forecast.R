plot_forecast <- function(workingGLM,sim_name){
  ###LOAD FORECAST FOR ANALYSIS
  load(file = paste0(workingGLM,sim_name,'_EnKF_output.Rdata'))

  ###PLOT FORECAST
  pdf(paste0(workingGLM,sim_name,'_forecast.pdf'))
  par(mfrow=c(4,3))
  
  z = z_obs
  nmembers = dim(x)[2]
  z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol =length(full_time)))
  nMETmembers =21
  
  for(i in 1:nlayers_init){
    model = i
    if(length(which(z_states[1,] == i) > 0)){
      obs = which(z_states[1,] == i)
    }else{
      obs = NA
    }
    #if(!is.na(obs)){
    ylim = range(c(x[,,],c(z[,])),na.rm = TRUE)
    as.POSIXct(full_time)
    plot(as.POSIXct(full_time_day),x[,1,model],type='l',ylab='water temperature (celsius)',xlab='time step (day)',main = paste('depth: ',the_depths_init[i],' m',sep=''),ylim=ylim)
    if(nmembers > 1){
      for(m in 2:nmembers){
        points(as.POSIXct(full_time_day),x[,m,model],type='l')
      }
    }
    if(!is.na(obs)){
      tmp = z[,obs]
      tmp[is.na(tmp)] = -999
      points(as.POSIXct(full_time_day),tmp,col='red',pch=19,cex=1.0)
    }
    
    abline(v = as.POSIXct(full_time_day[1+hist_days]))
    #}
  }
  
  ###PLOT HISTOGRAMS OF FORECAST
  par(mfrow=c(2,3))
  if(forecast_days > 6){
    xlim<- range(c(x[1+hist_days+7,,obs_index[1]],z[1+hist_days+7,1]),na.rm = TRUE)
    hist(x[1+hist_days+7,,obs_index[1]],main='0.1m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+7,1],col='red')
    xlim<- range(c(x[1+hist_days+7,,obs_index[5]],z[1+hist_days+7,5]),na.rm = TRUE)
    hist(x[1+hist_days+7,,obs_index[5]],main='4m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+7,5],col='red')
    xlim<- range(c(x[1+hist_days+7,,obs_index[10]],z[1+hist_days+7,10]),na.rm = TRUE)
    hist(x[1+hist_days+7,,obs_index[10]],main='9m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+7,10],col='red')
  }
  if(forecast_days > 13){
    xlim<- range(c(x[1+hist_days+14,,obs_index[1]],z[1+hist_days+14,1]),na.rm = TRUE)
    hist(x[1+hist_days+14,,obs_index[1]],main='0.1m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+14,1],col='red')
    xlim<- range(c(x[1+hist_days+14,,obs_index[5]],z[1+hist_days+14,5]),na.rm = TRUE)
    hist(x[1+hist_days+14,,obs_index[5]],main='4m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+14,5],col='red')
    xlim<- range(c(x[1+hist_days+14,,obs_index[10]],z[1+hist_days+14,10]),na.rm = TRUE)
    hist(x[1+hist_days+14,,obs_index[10]],main='9m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[1+hist_days+14,10],col='red')
  }
  
  ###PLOT NOAA MET TO VIEWING 
  d = read.csv(file.path(workingGLM,met_file_names[1]))
  air_temp = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  ShortWave = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  LongWave = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  RelHum = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  WindSpeed = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  Rain = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  
  y <- as.POSIXct(d$time)
  
  for(ens in 1:nMETmembers){
    d = read.csv(paste0(workingGLM,met_file_names[ens]))
    air_temp[ens,] = d$AirTemp
    ShortWave[ens,] = d$ShortWave
    LongWave[ens,] = d$LongWave
    RelHum[ens,] = d$RelHum
    WindSpeed[ens,] = d$WindSpeed
    Rain[ens,] = d$Rain
  }
  par(mfrow=c(2,3))
  ylim = range(c(air_temp))
  plot(y,air_temp[1,],type='l',ylab='Air Temp',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,air_temp[m,],type='l')
    }
  }
  
  ylim = range(c(ShortWave),na.rm = TRUE)
  plot(y,ShortWave[1,],type='l',ylab='Shortwave',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,ShortWave[m,],type='l')
    }
  }
  
  ylim = range(c(LongWave),na.rm = TRUE)
  plot(y,LongWave[1,],type='l',ylab='Longwave',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,LongWave[m,],type='l')
    }
  }
  
  ylim = range(c(RelHum))
  plot(y,RelHum[1,],type='l',ylab='Rel Hum',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,RelHum[m,],type='l')
    }
  }
  
  ylim = range(c(WindSpeed))
  plot(y,WindSpeed[1,],type='l',ylab='Wind Speed',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,WindSpeed[m,],type='l')
    }
  }
  
  ylim = range(c(Rain),na.rm = TRUE)
  plot(y,Rain[1,],type='l',ylab='Rain',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points(y,Rain[m,],type='l')
    }
  }
  
 
  
  dev.off()
}