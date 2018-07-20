if (!"mvtnorm" %in% installed.packages()) install.packages("mvtnorm")
if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4")
if (!"glmtools" %in% installed.packages()) install.packages('glmtools', repos=c('http://cran.rstudio.com', 'http://owi.usgs.gov/R'))
library(mvtnorm)
library(glmtools)
library(ncdf4)
library(lubridate)

###### OPTIONS ######

run_forecast<-function(first_day= '2018-07-10 00:00:00', sim_name = NA, hist_days = 1,forecast_days = 5){
  
  ###RUN OPTIONS
  Folder <- '/Users/quinn/Dropbox/Research/SSC_forecasting/SSC_forecasting/'
  nEnKFmembers <- 4
  nMETmembers <- 21
  nmembers = nEnKFmembers*nMETmembers
  
  use_CTD <- FALSE
  include_wq <- FALSE
  cold_start <- TRUE
  NO_UNCERT <- FALSE
  ADD_NOISE_TO_OBS <- FALSE
  USE_OBS_DEPTHS <- FALSE
  USE_OBS_CONTRAINT <- TRUE

  ###CREATE TIME VECTOR
  begin_sim  <- as.POSIXct(first_day)
  total_days <- hist_days + forecast_days
  end_sim <- begin_sim + total_days*24*60*60
  start_forecast_step <- hist_days
  forecast_start_time <- begin_sim + (start_forecast_step)*24*60*60
  if(day(forecast_start_time) < 10){
    forecast_day <- paste0('0',day(forecast_start_time))
  }else{
    forecast_day <- paste0(day(forecast_start_time))
  }
  if(month(forecast_start_time) < 10){
    forecast_month <- paste0('0',month(forecast_start_time))
  }else{
    forecast_month <- paste0(month(forecast_start_time))
  }
  full_time <- seq(begin_sim, end_sim, by = "1 day") # grid
  full_time <- strftime(full_time, format="%Y-%m-%d %H:%M")
  full_time_day <- strftime(full_time, format="%Y-%m-%d")
  full_time_hour_obs <- seq(as.POSIXct(full_time[1]), as.POSIXct(full_time[1+hist_days]), by = "1 hour") # grid
  nsteps <- length(full_time)
  
  ###CREATE DIRECTORY PATHS AND STRUCTURE
  workingGLM <- paste0(Folder,'/GLM_working/')  
  unlink(paste0(workingGLM,'*'),recursive = FALSE)    #Clear out temp GLM working directory
  
  ###LOAD SHARE R FUNCTIONS
  source(paste0(Folder,'/Rscripts/mcmc_enkf_shared_functions.R'))
  source(paste0(Folder,'/Rscripts/create_obs_met_input.R'))
  source(paste0(Folder,'/Rscripts/extract_temp_chain.R'))
  source(paste0(Folder,'/Rscripts/process_GEFS2GLM_v2.R'))
  source(paste0(Folder,'/Rscripts/extract_temp_CTD.R'))
  
  ###SHARED GLM LIBRARIES
  #Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH= paste(pathGLM,'/glm_lib_files/',sep=''))
  Sys.setenv(PATH='/opt/local/bin:/opt/local/sbin:/Users/quinn/anaconda2/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/local/bin')
  #system(paste('export DYLD_FALLBACK_LIBRARY_PATH=~',pathGLM,'/glm_lib_files:$DYLD_FALLBACK_LIBRARY_PATH',sep=''))
  
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
  obs_met_outfile <- paste0(workingGLM,'GLM_met.csv')
  create_obs_met_input(fname = met_obs_fname,outfile=obs_met_outfile,full_time_hour_obs)
  
  ###CREATE FUTURE MET FILES
  in_directory <- workingGLM
  out_directory <- workingGLM
  file_name <- forecast_base_name
  process_GEFS2GLM(in_directory,out_directory,file_name)
  met_file_names <- rep(NA,nMETmembers)
  for(i in 1:nMETmembers){
    met_file_names[i] <- paste(met_base_file_name,i,'.csv',sep='')
  }
  spillway_outflow_file_name <- paste0('FCR_spillway_outflow_',forecast_base_name,'.csv')
  inflow_file_name <- paste0('FCR_inflow_',forecast_base_name,'.csv')

  ###MOVE FILES AROUND
  SimFilesFolder <- paste0(Folder,'/sim_files/')
  fl <- c(list.files(SimFilesFolder, full.names = TRUE))
  tmp <- file.copy(from = fl, to = workingGLM,overwrite = TRUE)
  
  ###SET UP RUN
  lake_depth_init <- 10.0
  the_sals_init <- 0.5
  OGM_doc_init <- 47.4
  OXY_oxy_init <- 300.62
  CAR_dic_init <- 59.1
  NIT_amm_init <- 0.69
  NIT_nit_init <- 0.05
  PHS_frp_init <- 0.07
  CAR_ch4_init <- 0.58
  
  #Parameters
  Kw <- 0.86
  coef_mix_conv <- 0.165
  coef_wind_stir <- 0.23
  coef_mix_shear <- 0.2
  coef_mix_turb <- 0.51
  coef_mix_KH <- 0.3
  coef_mix_hyp <- 0.5
  wind_factor <- 1
  sw_factor <- 0.7
  lw_factor <- 1
  at_factor <- 1
  rh_factor <- 1
  rain_factor <- 1
  cd <- 0.000001
  ce <- 0.0013
  ch <- 0.0013
  
  #PROCESS TEMPERATURE OBSERVATIONS
  if(!use_CTD){
    obs_temp <- extract_temp_chain(fname = catwalk_fname)
  }else{
    obs_temp <- extract_temp_CTD(fname = ctd_fname)
  }
  
  #mg/L (obs) -> mol/m3 * 31.25
  obs_do <- extract_do_chain(fname = catwalk_fname)
  
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
  
  #SET UP GLM VARIANCE PARAMETERS
  #RMVNORM USES VARIANCE RATHER THAN STANDARD DEVIATATION
  variance_depths <- c(0.1,1,2,3,4,5,6,7,8,9,9.33)
  variance_values <- rep(0.5,length(variance_depths))
  #variance_values <- c(0.12,0.13,0.14,0.10,0.25,0.43,0.06,0.01,0.03,0.04,0.05)*10
  #variance_values <- rep(0.5^2,length(variance_depths))
  variance_inter <- approxfun(variance_depths,variance_values,rule=2)
  temps_variance <- variance_inter(the_depths_init)
  
  temps_variance_init <- 0.001^2
  temps_variance_init <- 1^2
  
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
    nstates <- nlayers_init*num_wq
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
  
  #Process error 
  cross_var <-0.0
  Qt_init <- diag(temps_variance_init, nstates)
  Qt <- diag(temps_variance, nstates)
  for(s in 1:nstates){
    if(s == 1){
      Qt[s,s+1] <- cross_var
    }else if(s == nstates){
      Qt[s,s-1] <- cross_var
    }else{
      Qt[s,s-1] <- cross_var
      Qt[s,s+1] <- cross_var
    }
  }
  #Measurement error 
  psi <- rep(0.0001,length(obs_index))
  
  ### INITILIZE FIRST TIME STEP
  restart_present <- file.exists(paste0(workingGLM,'restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'.csv'))
  x <- array(NA,dim=c(nsteps,nmembers,nstates))
  #Initial conditions
  if(cold_start){
    if(include_wq){
      x <- array(NA,dim=c(nsteps,nmembers,nstates))
      x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init,do_init), sigma=as.matrix(Qt_init))
      if(NO_UNCERT){
        for(m in 1:nmembers){
          x[1,m,] <- c(the_temps_init,do_init)
        }
      }
    }else{
      
      x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init), sigma=as.matrix(Qt_init))
      if(NO_UNCERT){
        for(m in 1:nmembers){
          x[1,m,] <- the_temps_init
        }
      }  
    }
    if(!restart_present){
      write.csv(x[1,,],paste0(workingGLM,'restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'),row.names = FALSE)
    }
  }
  
  #THIS ALLOWS THE EnKF TO BE RESTARTED FROM YESTERDAY'S RUN
  if(restart_present){
    x_previous <- read.csv(paste0(workingGLM,'restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'.csv'))
  }else{
    x_previous <- read.csv(paste0(workingGLM,'restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'))
  }
  
  if(dim(x[1,,])[1] != dim(x_previous)[1] | dim(x[1,,])[2] != dim(x_previous)[2]){
    print('ERROR: Dimension of the restart file are not correct. State variables or ensemble mismatch?')
    print('Need to fix if the states match but the ensemble number does not')
  }
  
  #Set initial conditions
  x[1,,] <- as.matrix(x_previous)
  
  
  #Matrix to store ensemble specific deviations and innovations
  dit <- array(NA,dim=c(nmembers,nstates))
  #dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation
  surface_height <- array(NA,dim=c(nsteps,nmembers))
  surface_height[1,] <- lake_depth_init

  ###START EnKF
  met_index <- 1
  for(i in 2:nsteps){
    
    #1) Update GLM NML files to match the current day of the simulation
    curr_start <- (full_time[i-1])
    curr_stop <- (full_time[i])
    update_time(start_value  = curr_start, stop_value = curr_stop,origNML)
    origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
    
    setwd(workingGLM)
    
    #Create array to hold GLM predictions for each ensemble
    x_star <- array(NA, dim = c(nmembers,nstates))
    for(m in 1:nmembers){
      
      #2) Use x[i-1,m,] to update GLM NML files for initial temperature at each depth
      tmp <- update_temps(curr_temps = x[i-1,m,temp_start:temp_end],origNML,the_depths_init)
      
      if(include_wq){
        wq_init_vals <- c(rep(OGM_doc_init,nlayers_init),x[i-1,m,do_start:do_end],rep(CAR_dic_init,nlayers_init),rep(NIT_amm_init,nlayers_init),rep(NIT_nit_init,nlayers_init),rep(PHS_frp_init,nlayers_init),rep(CAR_ch4_init,nlayers_init),CYANOPCH1_init_depth)
        update_var(wq_init_vals,origNML,'wq_init_vals')
      }
      
      #ALLOWS THE LOOPING THROUGH NOAA ENSEMBLES
      if(i > (hist_days+1)){
        update_var(met_file_names[met_index],origNML,'meteo_fl')
        origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
        update_var(paste0('FCR_inflow_',forecast_base_name,'.csv'),origNML,'inflow_fl')
        origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
        update_var(paste0('FCR_spillway_outflow_',forecast_base_name,'.csv'),origNML,'outflow_fl')
      }else{
        update_var(obs_met_outfile,origNML,'meteo_fl')
        origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
        update_var(paste0('FCR_inflow_',forecast_base_name,'.csv'),origNML,'inflow_fl')
        origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
        update_var(paste0('FCR_spillway_outflow_',forecast_base_name,'.csv'),origNML,'outflow_fl')
        
      }
      
      if(i == (hist_days+1)){
        write.csv(x[1-1,,],paste0(workingGLM,'restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'.csv'),row.names = FALSE)
      }
      
      #3) Use GLM NML files to run GLM for a day
      origNML <- read_nml(file.path(workingGLM,'glm3.nml'))
      system(paste0(workingGLM,"/glm"))
      
      #4) Fill x_star with temperatures from GLM
      GLMtemps <- get_glm_nc_var(ncFile = 'output.nc',z_out = the_depths_init, var = 'temp')
      x_star[m,temp_start:temp_end] = GLMtemps
      if(include_wq){
        GLMdo <- get_glm_nc_var(ncFile = 'output.nc',z_out = the_depths_init, var = 'OXY_oxy')
        x_star[m,do_start:do_end] <- GLMdo
      }
      
      #INCREMENT THE MET_INDEX TO MOVE TO THE NEXT NOAA ENSEMBLE
      met_index = met_index + 1
      if(met_index > nMETmembers){
        met_index <- 1
      }
    }
    
    #Corruption [nmembers x nstates] 
    NQt <- rmvnorm(n=nmembers, sigma=as.matrix(Qt))
    
    #Matrix Corrupted state estimate [nmembers x nstates]
    x_corr <- x_star + NQt
    
    #Obs for time step
    print(c(i,dim(z)))
    z_index <- which(!is.na(z[i,]))

    
    #if no observations at a time step then just propogate model uncertainity
    if(length(z_index) == 0 | i > (hist_days+1)){
      x[i,,] <- x_corr
      if(NO_UNCERT){
        x[i,,] <- x_star
      }
      
    }else{
      
      #if observation then calucate Kalman adjustment
      zt <- z[i,z_index]
      z_states_t <- z_states[i,z_index]
      yit <- array(NA,dim=c(nmembers,length(zt)))
      
      #Assign which states have obs in the time step
      H <- array(0,dim=c(length(zt),nstates))
      for(j in 1:length(z_index)){
        H[j,z_states_t[j]] <- 1
      }
      
      #Extract the data uncertainity for the data types present during the time-step
      if(length(z_index)>1){
        psi_t <- diag(psi[z_index])
      }else{
        #Special case where there is only one data type during the time-step
        psi_t <- psi[z_index]
      }
      
      #Ensemble mean
      ens_mean <- apply(x_corr, 2, mean)
      #ens_mean_star = apply(x_star, 2, mean) #Adaptive noise estimation
      
      #Loop through ensemble members
      for(m in 1:nmembers){  
        
        #Ensemble specific deviation
        dit[m,] <- x_corr[m,]-ens_mean
        #dit_star[m,] x_star[m,] - ens_mean_star #Adaptive noise estimation
        
        #Observational uncertainity
        N_psi = rmvnorm(n=1,sigma=as.matrix(psi_t))
        
        #Ensemble specific innovations
        yit[m,] <- t(zt - crossprod(t(H),(x_corr[m,])) + t(N_psi))
        
        #Ensemble specific estimate and innovation covariance
        if(m == 1){
          Pit <- tcrossprod(dit[m,])
          #Pit_star = tcrossprod(dit_star[m,])  #Adaptive noise estimation
          Sit <- tcrossprod(yit[m,])
        }else{
          Pit <- tcrossprod(dit[m,]) +  Pit
          #Pit_star = tcrossprod(dit_star[m,]) + Pit_star #Adaptive noise estimation
          Sit <- tcrossprod(yit[m,]) +  Sit
        }
      }
      
      #estimate covariance
      Pt <- Pit/nmembers
      #Pt_star = Pit_star/nmembers  #Adaptive noise estimation
      #Innovations covariance
      St <- Sit/nmembers
      
      #Kalman gain
      Kt <- crossprod(t(crossprod(Pt,t(H))), solve(St, tol=1e-30))
      
      #Adaptive noise estimation
      beta <- 0.55
      #Gammat = (1 - beta)*crossprod(t(crossprod(Pt,t(H))), solve(H, tol=1e-30))
      #Qt_hat = Gammat(St - H*Pt_star*H - t(N_psi))*Gammat
      #Q = alpha*Q+ (1-alpha)*Qt_hat
      #Ensemble specific updated state
      for(m in 1:nmembers){
        x[i,m,] <- t((x_corr[m,]) + crossprod(t(Kt), yit[m,]))
        if(NO_UNCERT){
          x[i,m,] <- x_star[m,]
        }
      }
      
      #print('here')
      #print(curr_start)
      #print(Kt)
      #print(zt[1])
      #print(x_corr[1:10,1])
      #print(x[i,1:10,1])
      #print(mean(x[i,,1]))
      #readline(prompt="Press [enter] to continue")
      if(length(which(is.na(x[i,,]))) > 0){dies = i}
    }
  }
  
  
  ###SAVE FORECAST
  save(x,full_time,the_depths_init,file = paste0(workingGLM,sim_name,'_EnKF_output.Rdata'))
  
  ###LOAD FORECAST FOR ANALYSIS
  load(file = paste0(workingGLM,sim_name,'_EnKF_output.Rdata'))
  
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
  par(mfrow=c(3,2))
  if(forecast_days > 13){
    hist(x[1+hist_days+14,,1],main='surface temperature 14 days in the future',xlab='Temperature')
    abline(v= z[1+hist_days+14,obs_index[1]],col='red')
    #hist(x[hist_days+14,,1]-x[1,,1],main='surface temperature change 14 days in the future',xlab='Temperature')
    hist(x[1+hist_days+14,,10],main='4m temperature 14 days in the future',xlab='Temperature')
  }
  if(forecast_days > 6){
    hist(x[1+hist_days+7,,1],main='surface temperature 7 days in the future',xlab='Temperature')
    abline(v= z[1+hist_days+7,obs_index[1]],col='red')
    #hist(x[hist_days+4,,1]-x[1,,1],main='surface temperature change 7 days in the future',xlab='Temperature')
    #abline(v= z[hist_days+4,obs_index[1]]-21.8,col='red')
    hist(x[1+hist_days+7,,10],main='4m temperature 7 days in the future',xlab='Temperature')
    abline(v= z[1+hist_days+7,4],col='red')
  }
  
  dev.off()
  
  ###ARCHIVE AND CLEAN UP FORECAST
  unlink(paste0(workingGLM,'FCRmet.csv'),recursive = FALSE)
  unlink(paste0(workingGLM,'Catwalk.csv'),recursive = FALSE)
  unlink(paste0(workingGLM,'20180711gep_all_00z.csv'),recursive = FALSE)
  time_of_forecast <- paste0(year(Sys.time()),month(Sys.time()),day(Sys.time()),'_',hour(Sys.time()),'_',(minute(Sys.time())))
  forecast_archive_dir <- paste0(Folder,'/Forecasts/','forecast_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_',time_of_forecast)
  dir.create(forecast_archive_dir)
  files <- list.files(paste0(workingGLM))
  tmp <- file.copy(files, forecast_archive_dir)

}



