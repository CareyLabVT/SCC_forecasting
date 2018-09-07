plot_forecast_netcdf <- function(pdf_file_name,output_file,catwalk_fname,include_wq,code_location){
  
  code_location <- '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SSC_forecasting/Rscripts'
  source(paste0(code_location,'/extract_temp_chain.R'))
  
  output_file <-'/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/test_forecast/FCR_hist_2018_7_15_201897_15_32.nc'
  include_wq <- FALSE
  data_location <- '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SCC_data'
  mia_location <- paste0(data_location,'/','mia-data')
  setwd(mia_location)
  system(paste0('git pull'))

  catwalk_fname <- paste0(mia_location,'/','Catwalk.csv')
  
  plot_summaries <- TRUE
  
  nc <- nc_open(output_file)
  
  t <- ncvar_get(nc,'time')
  full_time <- as.POSIXct(t, origin = '1970-01-01 00:00.00 UTC')
  full_time_day <- strftime(full_time, format="%Y-%m-%d")
  temp_mean <- ncvar_get(nc,'temp_mean')
  temp <- ncvar_get(nc,'temp')
    temp_upper <- ncvar_get(nc,'temp_upperCI')
    temp_lower  <- ncvar_get(nc,'temp_lowerCI')
    depths <- ncvar_get(nc,'z')
  Kw <- ncvar_get(nc,'Kw')
  forecasted <- ncvar_get(nc,'forecasted')
  
  nsteps <- length(full_time)
  forecast_index <- which(forecasted == 1)[1]
  nlayers <- length(depths)
  
  obs_temp <- extract_temp_chain(fname = catwalk_fname,full_time)
  for(i in 1:length(obs_temp$obs[,1])){
    for(j in 1:length(obs_temp$obs[1,])){
      if(obs_temp$obs[i,j] == 0 | is.na(obs_temp$obs[i,j]) | is.nan(obs_temp$obs[i,j])){
        obs_temp$obs[i,j] = NA
      } 
    }
  }
  TempObservedDepths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8,9)
  DoObservedDepths <- c(1,5,9)
  
  
  if(include_wq){
    nobs <- length(TempObservedDepths) + length(DoObservedDepths)
  }else{
    nobs <- length(TempObservedDepths)
  }
  
  #Observations for each observed state at each time step
  #an observation with at least 1 observation but without an observation in a time-step gets assigned an NA
  z <- t(matrix(rep(NA,nobs), nrow = nobs, ncol = nsteps))
  
  if(include_wq){
    z <- cbind(obs_temp$obs,obs_do$obs)
  }else{
    z <- cbind(obs_temp$obs) 
  }
  
  #FIGURE OUT WHICH DEPTHS HAVE OBSERVATIONS
  if(include_wq){
    obs_index <- rep(NA,length(TempObservedDepths)+length(DoObservedDepths))
    for(i in 1:length(TempObservedDepths)){
      obs_index[i] <- which.min(abs(depths - TempObservedDepths[i]))
    }
    for(i in 1:length(DoObservedDepths)){
      obs_index[length(TempObservedDepths)+i] <- length(depths) + which.min(abs(depths - DoObservedDepths[i]))
    }
  }else{
    obs_index <- rep(NA,length(TempObservedDepths))
    for(i in 1:length(TempObservedDepths)){
      obs_index[i] <- which.min(abs(depths - TempObservedDepths[i]))
    } 
  }
  
  #Matrix for knowing which state the observation corresponds to
  z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol = nsteps))
  
  nMETmembers =21
  pdf(paste0(workingGLM,'/',pdf_file_name))
  par(mfrow=c(4,3))
  
  for(i in 1:nlayers){
    model = i
    if(length(which(z_states[1,] == i) > 0)){
      obs = which(z_states[1,] == i)
    }else{
      obs = NA
    }
    
    ylim = range(c(temp_mean,temp_lower,temp_upper,c(z[,])),na.rm = TRUE)
    if(plot_summaries){
      plot(as.POSIXct(full_time_day),temp_mean[,model],type='l',ylab='water temperature (celsius)',xlab='time step (day)',main = paste('depth: ',depths[i],' m',sep=''),ylim=ylim)
      points(as.POSIXct(full_time_day),temp_upper[,model],type='l',lty='dashed')
      points(as.POSIXct(full_time_day),temp_lower[,model],type='l',lty='dashed')
    }else{
      plot(as.POSIXct(full_time_day),temp[,1,model],type='l',ylab='water temperature (celsius)',xlab='time step (day)',main = paste('depth: ',depths[i],' m',sep=''),ylim=ylim)
      if(nmembers > 1){
        for(m in 2:nmembers){
          points(as.POSIXct(full_time_day),temp[,m,model],type='l')
        }
      }
    }
    if(!is.na(obs)){
      tmp = z[,obs]
      tmp[is.na(tmp)] = -999
      points(as.POSIXct(full_time_day),tmp,col='red',pch=19,cex=1.0)
    }
    
    abline(v = as.POSIXct(full_time_day[forecast_index]))
    #}
  }
  
  ###PLOT OF PARAMETERS IF FIT
  plot(rowMeans(Kw[,]),xlab ='time step (day)',ylab = 'Kw parameter')
  
  
  ###PLOT HISTOGRAMS OF FORECAST
  par(mfrow=c(2,3))
  if(!is.na(forecast_index)){
  if(length(which(forecast_index == 1)) > 6){
    xlim<- range(c(temp[forecast_index+7,,obs_index[1]],z[forecast_index,1]),na.rm = TRUE)
    hist(temp[forecast_index+7,,obs_index[1]],main='0.1m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[forecast_index+7,1],col='red')
    xlim<- range(c(temp[forecast_index+7,,obs_index[5]],z[forecast_index+7,5]),na.rm = TRUE)
    hist(temp[forecast_index+7,,obs_index[5]],main='4m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[forecast_index+7,5],col='red')
    xlim<- range(c(temp[forecast_index+7,,obs_index[10]],z[forecast_index+7,10]),na.rm = TRUE)
    hist(temp[forecast_index+7,,obs_index[10]],main='9m temp. 7 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[forecast_index+7,10],col='red')
  }
  if(length(which(forecast_index == 1)) > 13){
    xlim<- range(c(temp[forecast_index+14,,obs_index[1]],z[forecast_index+14,1]),na.rm = TRUE)
    hist(temp[forecast_index+14,,obs_index[1]],main='0.1m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[forecast_index+14,1],col='red')
    xlim<- range(c(temp[forecast_index+14,,obs_index[5]],z[forecast_index+14,5]),na.rm = TRUE)
    hist(temp[forecast_index+14,,obs_index[5]],main='4m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= temp[forecast_index+14,5],col='red')
    xlim<- range(c(temp[forecast_index+14,,obs_index[10]],z[forecast_index+14,10]),na.rm = TRUE)
    hist(temp[forecast_index+14,,obs_index[10]],main='9m temp. 14 days forecast',xlab='Temperature',xlim=xlim)
    abline(v= z[forecast_index+14,10],col='red')
  }
  }
}

