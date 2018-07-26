evaluate_forecast <- function(){
  ###LOAD FORECAST FOR ANALYSIS
  load(file = paste0(workingGLM,sim_name,'_EnKF_output.Rdata'))
  
  full_time_hour_obs <- seq(as.POSIXct(full_time[1]), as.POSIXct(full_time[length(full_time)]), by = "1 hour") # grid
  
  ###LOAD SHARE R FUNCTIONS
  source(paste0(Folder,'/Rscripts/mcmc_enkf_shared_functions.R'))
  source(paste0(Folder,'/Rscripts/create_obs_met_input.R'))
  source(paste0(Folder,'/Rscripts/extract_temp_chain.R'))
  source(paste0(Folder,'/Rscripts/process_GEFS2GLM_v2.R'))
  source(paste0(Folder,'/Rscripts/extract_temp_CTD.R'))
  source(paste0(Folder,'/Rscripts/create_inflow_outflow_file.R'))
  
  ###SET FILE NAMES
  forecast_base_name <- paste0(year(forecast_start_time),forecast_month,forecast_day,'gep_all_00z',sep='')
  catwalk_fname <-  paste0(workingGLM,'Catwalk.csv')
  met_obs_fname <-paste0(workingGLM,'FCRmet.csv')
  ctd_fname <- '/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/test_data/070218_fcr50.csv' 
  met_base_file_name <- paste0('met_hourly_',forecast_base_name,'_ens')
  if(is.na(sim_name)){
    sim_name <- paste0('historical_start_',year(first_day),'_',month(first_day),'_',day(first_day),'_forecast_start_',paste0(year(forecast_start_time),forecast_month,forecast_day))
  }
  
  ###DOWNLOAD FILES TO WORKING DIRECTORY
  download.file('https://github.com/CareyLabVT/SCCData/raw/carina-data/FCRmet.csv',paste0(workingGLM,'FCRmet.csv'))
  download.file('https://github.com/CareyLabVT/SCCData/raw/mia-data/Catwalk.csv',paste0(workingGLM,'Catwalk.csv'))
  download.file(paste0('https://github.com/CareyLabVT/SCCData/raw/noaa-data/',forecast_base_name,'.csv'),paste0(workingGLM,forecast_base_name,'.csv'))
  
  ###CREATE HISTORICAL MET FILE
  obs_met_outfile <- paste0(workingGLM,'GLM_met_eval.csv')
  create_obs_met_input(fname = met_obs_fname,outfile=obs_met_outfile,full_time_hour_obs)
  
  
  obs_temp <- extract_temp_chain(fname = catwalk_fname,full_time)
  
  #mg/L (obs) -> mol/m3 * 31.25
  obs_do <- extract_do_chain(fname = catwalk_fname,full_time)
  
  #KLUDGE TO GET WORKING
  TempObservedDepths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8,9)
  init_temps1 <- obs_temp$obs[1,]
  
  DoObservedDepths <- c(1,5,9)
  
  temp_inter <- approxfun(TempObservedDepths,init_temps1,rule=2)
  
  #SET UP INITIAL CONDITIONS
  if(USE_OBS_DEPTHS){
    nlayers_init <- length(TempObservedDepths)
    the_depths_init <- TempObservedDepths
    the_temps_init <- init_temps1
  }else{
    the_depths_init <- c(0.1, 0.33, 0.66, 1.00, 1.33,1.66,2.00,2.33,2.66,3.0,3.33,3.66,4.0,4.33,4.66,5.0,5.33,5.66,6.0,6.33,6.66,7.00,7.33,7.66,8.0,8.33,8.66,9.00,9.33)
    nlayers_init <- length(the_depths_init)
    the_temps_init <- temp_inter(the_depths_init)
    do_init <- rep(NA,length(the_depths_init))
    do_init[1:13] <- obs_do$obs[1,1]
    do_init[14:23] <- obs_do$obs[1,2]
    do_init[24:29] <- obs_do$obs[1,3]
  }
  
  temp_start <- 1
  temp_end <- length(the_depths_init)
  do_start <- temp_end+1
  do_end <- temp_end + (length(the_depths_init))
  
  #TEMPORARY TO ADD PHYTOS
  CYANOPCH1_init_depth <- c(rep(1,nlayers_init))
  CYANONPCH2_init_depth <- c(rep(1,nlayers_init)) 
  CHLOROPCH3_init_depth <- c(rep(1,nlayers_init)) 
  DIATOMPCH4_init_depth <- c(rep(1,nlayers_init))
  GREENCH5_init_depth <- c(rep(1,nlayers_init))
  
  #UPDATE NML WITH PARAMETERS AND INITIAL CONDITIONS
  wq_init_vals <- c(rep(OGM_doc_init,nlayers_init),do_init,rep(CAR_dic_init,nlayers_init),rep(NIT_amm_init,nlayers_init),rep(NIT_nit_init,nlayers_init),rep(PHS_frp_init,nlayers_init),rep(CAR_ch4_init,nlayers_init),CYANOPCH1_init_depth)
  update_var(wq_init_vals,origNML,'wq_init_vals')
  update_var(rep(the_sals_init,nlayers_init),origNML,'the_sals')
  update_var(lake_depth_init,origNML,'lake_depth')
  update_var(nlayers_init,origNML,'num_depths')
  update_var(the_temps_init,origNML,'the_temps')
  update_var(the_depths_init,origNML,'the_depths')
  
  update_var(Kw,origNML,'Kw')
  update_var(coef_mix_conv,origNML,'coef_mix_conv')
  update_var(coef_wind_stir,origNML,'coef_wind_stir')
  update_var(coef_mix_shear,origNML,'coef_mix_shear')
  update_var(coef_mix_turb,origNML,'coef_mix_turb')
  update_var(coef_mix_KH,origNML,'coef_mix_KH')
  update_var(coef_mix_hyp,origNML,'coef_mix_hyp')
  update_var(wind_factor,origNML,'wind_factor')
  update_var(sw_factor,origNML,'sw_factor')
  update_var(lw_factor,origNML,'lw_factor')
  update_var(at_factor,origNML,'at_factor')
  update_var(rh_factor,origNML,'rh_factor')
  update_var(rain_factor,origNML,'rain_factor')
  update_var(cd,origNML,'cd')
  update_var(ce,origNML,'ce')
  update_var(ch,origNML,'ch')
  
  #NUMBER OF STATE SIMULATED = SPECIFIED DEPTHS
  if(include_wq){
    nstates <- nlayers_init*(1+num_wq)
  }else{
    nstates <- nlayers_init 
  }
  
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
  
  z_obs <- z
  if(!USE_OBS_CONTRAINT){
    z[,] <- NA
  }
  
  #FIGURE OUT WHICH DEPTHS HAVE OBSERVATIONS
  if(include_wq){
    obs_index <- rep(NA,length(TempObservedDepths)+length(DoObservedDepths))
    for(i in 1:length(TempObservedDepths)){
      obs_index[i] <- which.min(abs(the_depths_init - TempObservedDepths[i]))
    }
    for(i in 1:length(DoObservedDepths)){
      obs_index[length(TempObservedDepths)+i] <- length(the_depths_init) + which.min(abs(the_depths_init - DoObservedDepths[i]))
    }
  }else{
    obs_index <- rep(NA,length(TempObservedDepths))
    for(i in 1:length(TempObservedDepths)){
      obs_index[i] <- which.min(abs(the_depths_init - TempObservedDepths[i]))
    } 
  }
  
  #Matrix for knowing which state the observation corresponds to
  z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol = nsteps))
  
  
  ###
  update_var(sw_factor,origNML,'sw_factor')
  update_var(wq_init_vals,origNML,'wq_init_vals')
  update_var(rep(the_sals_init,nlayers_init),origNML,'the_sals')
  update_var(lake_depth_init,origNML,'lake_depth')
  update_var(nlayers_init,origNML,'num_depths')
  update_var(the_temps_init,origNML,'the_temps')
  update_var(the_depths_init,origNML,'the_depths')
  update_time(start_value  = full_time[1], stop_value = full_time[length(full_time)],origNML)
  update_var('GLM_met_eval.csv',origNML,'meteo_fl')
  update_var(paste0('FCR_inflow.csv'),origNML,'inflow_fl')
  update_var(paste0('FCR_spillway_outflow.csv'),origNML,'outflow_fl')
  update_var(length(full_time),origNML,'num_days')
  system(paste0(workingGLM,"/glm"))
  
  glm_prediction <- get_temp(file = "output.nc", reference = "surface", z_out = the_depths_init)
  
  ###PLOT FORECAST
  pdf(paste0(workingGLM,sim_name,'_EnKF_output.pdf'))
  par(mfrow=c(4,3))
  
  z = z_obs
  
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
    points(as.POSIXct(glm_prediction$DateTime),glm_prediction[,1+i],col='lightblue',type='o')
    #}
  }
  
  
  ###PLOT NOAA MET TO VIEWING 
  d = read.csv(paste0(workingGLM,'met_hourly_',forecast_base_name,'_ens1.csv'))
  air_temp = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  ShortWave = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  LongWave = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  RelHum = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  WindSpeed = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  Rain = array(NA,dim=c(nMETmembers,length(d$AirTemp)))
  
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
  plot((1:ncol(air_temp))/24,air_temp[1,],type='l',ylab='Air Temp',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(air_temp))/24,air_temp[m,],type='l')
    }
  }
  
  ylim = range(c(ShortWave),na.rm = TRUE)
  plot((1:ncol(ShortWave))/24,ShortWave[1,],type='l',ylab='Shortwave',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(ShortWave))/24,ShortWave[m,],type='l')
    }
  }
  
  ylim = range(c(LongWave),na.rm = TRUE)
  plot((1:ncol(LongWave))/24,LongWave[1,],type='l',ylab='Longwave',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(LongWave))/24,LongWave[m,],type='l')
    }
  }
  
  ylim = range(c(RelHum))
  plot((1:ncol(RelHum))/24,RelHum[1,],type='l',ylab='Rel Hum',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(RelHum))/24,RelHum[m,],type='l')
    }
  }
  
  ylim = range(c(WindSpeed))
  plot((1:ncol(WindSpeed))/24,WindSpeed[1,],type='l',ylab='Wind Speed',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(WindSpeed))/24,WindSpeed[m,],type='l')
    }
  }
  
  ylim = range(c(Rain),na.rm = TRUE)
  plot((1:ncol(Rain))/24,Rain[1,],type='l',ylab='Rain',xlab = 'days in future',ylim=ylim)
  if(nMETmembers > 1){
    for(m in 2:nMETmembers){
      points((1:ncol(Rain))/24,Rain[m,],type='l')
    }
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
  
  
  dev.off()
  
  
  
}