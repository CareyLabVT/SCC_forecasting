plot_forecast_management <- function(pdf_file_name,output_file,catwalk_fname,include_wq,code_location,save_location,data_location,plot_summaries,PRE_SCC){
  library(ncdf4)
  library(lubridate)
  
  
  #code_location <- '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SSC_forecasting/Rscripts'
  #output_file <-'/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/test_forecast/FCR_betaV2_hist_2018_9_15_forecast_2018_9_16_2018916_9_47.nc'
  #include_wq <- FALSE
  #data_location <- '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SCC_data'
  source(paste0(code_location,'/extract_temp_chain.R'))
  source(paste0(code_location,'/extract_temp_CTD.R'))
  
  nc <- nc_open(output_file)
  
  t <- ncvar_get(nc,'time')
  depths <- ncvar_get(nc,'z')
  
  
  full_time <- as.POSIXct(t, origin = '1970-01-01 00:00.00 UTC', tz = 'EST5EDT')
  full_time_day <- strftime(full_time, format="%Y-%m-%d")
  
  full_time_past <- seq(full_time[1] - days(5), full_time[2], by = "1 day") # grid
  
  full_time_combined <- seq(full_time_past[1], full_time[length(full_time)], by = "1 day")
  
  full_time_plotting <- seq(full_time_past[1]-days(3), full_time[length(full_time)]+days(4), by = "1 day")
  
  
  mia_location <- paste0(data_location,'/','mia-data')
  setwd(mia_location)
  system(paste0('git pull'))
  
  catwalk_fname <- paste0(mia_location,'/','Catwalk.csv')
  
  obs_temp_past <- extract_temp_chain(fname = catwalk_fname,full_time_past)
  for(i in 1:length(obs_temp_past$obs[,1])){
    for(j in 1:length(obs_temp_past$obs[1,])){
      if(obs_temp_past$obs[i,j] == 0 | is.na(obs_temp_past$obs[i,j]) | is.nan(obs_temp_past$obs[i,j])){
        obs_temp_past$obs[i,j] = NA
      } 
    }
  }
  TempObservedDepths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8,9)
  
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
  
  temp_mean <- ncvar_get(nc,'temp_mean')
  temp <- ncvar_get(nc,'temp')
  temp_upper <- ncvar_get(nc,'temp_upperCI')
  temp_lower  <- ncvar_get(nc,'temp_lowerCI')
  depths <- ncvar_get(nc,'z')
  zone1temp <- ncvar_get(nc,'zone1temp')
  zone2temp <- ncvar_get(nc,'zone2temp')
  forecasted <- ncvar_get(nc,'forecasted')
  
  nsteps <- length(full_time)
  forecast_index <- which(forecasted == 1)[1]
  nlayers <- length(depths)
  
  focal_depths <- c(4,16,25)
  pdf(paste0(save_location,'/',pdf_file_name),width = 6, height = 10)
  par(mfrow=c(2,1))
  
  #PLOT OF TURNOVER PROBABILITY
  prob_zero <- rep(NA,length(seq(3,17,1)))
  for(i in 3:17){
    prob_zero[i-3] = 100*length(which(temp[i,,obs_index[1]] - temp[i,,obs_index[9]] < 1))/length((temp[i,,obs_index[1]]))
  }
  
  plot(full_time_plotting,rep(-99,length(full_time_plotting)),ylim=c(0,100),xlab = 'date',ylab = '% chance')
  title(paste0('Falling Creek Reservior\n',month(tmp_day),'/',day(tmp_day),'/',year(tmp_day), '\n\nTurnover forecast'),cex.main=0.9)
  
  points(full_time[3:17],prob_zero,type='o',ylim=c(0,100),xlab = 'date',ylab = 'Probablity of turnover')
  axis(1, at=full_time_plotting,las=2, cex.axis=0.7, tck=-0.01,labels=FALSE)
  abline(v = full_time_past[length(full_time_past)])
  text(full_time_past[length(full_time_past)-2],80,'past')
  text(full_time[4],80,'future')
  #HISTORICAL AND FUTURE TEMPERATURE
  depth_colors <- c("firebrick4","firebrick1","DarkOrange1","gold","greenyellow","medium sea green","sea green","DeepSkyBlue4","blue2","blue4")

  plot(full_time_plotting,rep(-99,length(full_time_plotting)),ylim=c(5,35),xlab = 'date',ylab = expression(~degree~C))
  title(paste0('Water temperature forecast'),cex.main=0.9)
  tmp_day <- full_time[-1][1]
  axis(1, at=full_time_plotting,las=2, cex.axis=0.7, tck=-0.01,labels=FALSE)
  
  for(i in 1:length(obs_index)){
    points(full_time_past, obs_temp_past$obs[,i],type='l',col=depth_colors[i],lwd=1.5)
    index <- which(obs_index[i]  == focal_depths)
    if(length(index) == 1){
      points(full_time[-1], temp_mean[-1,obs_index[i]],type='l',lty='dashed',col=depth_colors[i],lwd=1.5) 
      points(full_time[-1], temp_upper[-1,obs_index[i]],type='l',lty='dotted',col=depth_colors[i],lwd=1.5)
      points(full_time[-1], temp_lower[-1,obs_index[i]],type='l',lty='dotted',col=depth_colors[i],lwd=1.5)
    }
  }
  
  abline(v = full_time_past[length(full_time_past)])
  text(full_time_past[length(full_time_past)-2],30,'past')
  text(full_time[4],30.1,'future')
  if(temp_mean[length(temp_mean[,focal_depths[1]]),focal_depths[1]] == temp_mean[length(temp_mean[,focal_depths[2]]),focal_depths[2]] |
     temp_mean[length(temp_mean[,focal_depths[1]]),focal_depths[1]] == temp_mean[length(temp_mean[,focal_depths[3]]),focal_depths[3]] |
     temp_mean[length(temp_mean[,focal_depths[2]]),focal_depths[2]] == temp_mean[length(temp_mean[,focal_depths[3]]),focal_depths[3]]){
  text(full_time_plotting[length(full_time_plotting)-3], temp_mean[length(temp_mean[,focal_depths[1]]),focal_depths[1]], '1m', col='firebrick1')
  text(full_time_plotting[length(full_time_plotting)-2], temp_mean[length(temp_mean[,focal_depths[2]]),focal_depths[2]], '5m', col='medium sea green')
  text(full_time_plotting[length(full_time_plotting)-1], temp_mean[length(temp_mean[,focal_depths[3]]),focal_depths[3]], '8m', col='blue2')
  }else{
    text(full_time_plotting[length(full_time_plotting)-2], temp_mean[length(temp_mean[,focal_depths[1]]),focal_depths[1]], '1m', col='firebrick1')
    text(full_time_plotting[length(full_time_plotting)-2], temp_mean[length(temp_mean[,focal_depths[2]]),focal_depths[2]], '5m', col='medium sea green')
    text(full_time_plotting[length(full_time_plotting)-2], temp_mean[length(temp_mean[,focal_depths[3]]),focal_depths[3]], '8m', col='blue2') 
  }
  
    legend("left",c("0.1m","1m", "2m", "3m", "4m", "5m", "6m", "7m","8m", "9m"),
         text.col=c("firebrick4", "firebrick1", "DarkOrange1", "gold", "greenyellow", "medium sea green", "sea green",
                    "DeepSkyBlue4", "blue2", "blue4"), cex=1, y.intersp=1, x.intersp=0.001, inset=c(0,0), xpd=T, bty='n')
  legend('topright', c('mean','confidence bounds'), lwd=1.5, lty=c('dashed','dotted'),bty='n',cex = 1)

  
  dev.off()
}