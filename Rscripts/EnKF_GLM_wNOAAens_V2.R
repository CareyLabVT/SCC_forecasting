run_forecast<-function(first_day= '2018-07-06 00:00:00', sim_name = NA, hist_days = 1,forecast_days = 15,  spin_up_days = 0,restart_file = NA, Folder, forecast_location = NA,push_to_git=FALSE,data_location = NA, nEnKFmembers = NA){
  
  ###RUN OPTIONS
  include_wq <- FALSE
  num_pars <- 3
  
  USE_QT_MATRIX <- TRUE
  USE_CTD <- FALSE
  PRE_SCC <- FALSE
  
  USE_OBS_CONTRAINT <- TRUE
  NO_UNCERT <- FALSE
  
  
  #Parameters
  sw_factor <- 0.95
  lw_factor <- 0.95
  lake_depth_init <- 9.4  #not a modeled state
  kw_init <- 0.87
  kw_init_Qt <- 0.00000001
  zone2_temp <- 17
  zone1_temp <- 11
  zone1temp_init_Qt <- 0.0001
  zone2temp_init_Qt <- 0.0001
  
  #ERROR TERMS
  obs_error <- 0.0001 #NEED TO FIX
  thermo_depth_error <- 0.15
  temp_error <- 0.5
  b <- 0.5  #0 = error is distributed to observed states, 1 = error is distributed to unobserved states
  alpha <- 0.5
  
  ####################################################
  #### YOU WON'T NEED TO MODIFY ANYTHING BELOW HERE ##
  ####################################################
  
  # SET UP NUMBER OF ENSEMBLE MEMBERS
  nMETmembers <- 21
  if(is.na(nEnKFmembers) & is.na(restart_file)){ 
    nEnKFmembers <- 50
    nmembers <- nEnKFmembers*nMETmembers
  }else if(!is.na(restart_file)){
    nc <- nc_open(restart_file)
    nmembers <- length(ncvar_get(nc,'ens'))
    nc_close(nc)
  }else{
    nmembers <- nEnKFmembers*nMETmembers
  }
  
  ###DETECT THE PLATFORM###
  
  switch(Sys.info() [['sysname']],
         Linux = { machine <- 'unix' },
         Darwin = { machine <- 'mac' })
  
  ###INSTALL PREREQUISITES##
  
  #INSTALL libnetcdf
  if(machine == 'unix') {
    system("if [ $(dpkg-query -W -f='${Status}' libnetcdf-dev 2>/dev/null | grep -c 'ok installed') -eq 0 ]; then sudo apt update && sudo apt install libnetcdf-dev; fi;")
    Sys.setenv(LD_LIBRARY_PATH=paste("../glm/unix/", Sys.getenv("LD_LIBRARY_PATH"),sep=":"))
  }
  
  #--CREATE TIME VECTOR---
  # The simulations are run from 00:00:00 GMT time 
  # so that they directly interface with the NOAA forecast
  # The output is converted back to local time before being saved
  
  begin_sim  <- as.POSIXct(first_day,tz = reference_tzone)
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
  full_time_local <- with_tz(full_time,tzone = 'EST5EDT')
  full_time <- strftime(full_time, format="%Y-%m-%d %H:%M",tz = reference_tzone)
  full_time_local <- strftime(full_time_local, format="%Y-%m-%d %H:%M")
  full_time_day <- strftime(full_time, format="%Y-%m-%d")
  full_time_day_local <- strftime(full_time_local, format="%Y-%m-%d")
  full_time_hour_obs <- seq(as.POSIXct(full_time[1],tz = reference_tzone), as.POSIXct(full_time[length(full_time)],tz = reference_tzone), by = "1 hour") # grid
  nsteps <- length(full_time)
  
  ###CREATE DIRECTORY PATHS AND STRUCTURE
  workingGLM <- paste0(Folder,'/','GLM_working')  
  print(workingGLM)
  unlink(paste0(workingGLM,'/*'),recursive = FALSE)    #Clear out temp GLM working directory
  
  ###LOAD SHARE R FUNCTIONS
  source(paste0(Folder,'/','Rscripts/mcmc_enkf_shared_functions.R'))
  source(paste0(Folder,'/','Rscripts/create_obs_met_input.R'))
  source(paste0(Folder,'/','Rscripts/extract_temp_chain.R'))
  source(paste0(Folder,'/','Rscripts/process_GEFS2GLM_v2.R'))
  source(paste0(Folder,'/','Rscripts/extract_temp_CTD.R'))
  source(paste0(Folder,'/','Rscripts/create_inflow_outflow_file.R'))
  source(paste0(Folder,'/','Rscripts/plot_forecast_netcdf.R'))
  source(paste0(Folder,'/','Rscripts/archive_forecast.R'))
  source(paste0(Folder,'/','Rscripts/write_forecast_netcdf.R')) 
  ###SHARED GLM LIBRARIES
  #Sys.setenv(DYLD_FALLBACK_LIBRARY_PATH= paste(pathGLM,'/glm_lib_files/',sep=''))
  #Sys.setenv(PATH='/opt/local/bin:/opt/local/sbin:/Users/quinn/anaconda2/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/local/bin')
  #system(paste('export DYLD_FALLBACK_LIBRARY_PATH=~',pathGLM,'/glm_lib_files:$DYLD_FALLBACK_LIBRARY_PATH',sep=''))
  
  ###SET FILE NAMES
  forecast_base_name <- paste0(year(forecast_start_time),forecast_month,forecast_day,'gep_all_00z')
  catwalk_fname <-  paste0(workingGLM,'/','Catwalk.csv')
  met_obs_fname <-paste0(workingGLM,'/','FCRmet.csv')
  met_base_file_name <- paste0('met_hourly_',forecast_base_name,'_ens')
  if(is.na(sim_name)){
    sim_name <- paste0(year(full_time_local[1]),'_',month(full_time_local[1]),'_',day(full_time_local[1]))
  }
  
  ###DOWNLOAD FILES TO WORKING DIRECTORY
  mia_location <- paste0(data_location,'/','mia-data')
  setwd(mia_location)
  system(paste0('git pull'))
  carina_location <- paste0(data_location,'/','carina-data')
  setwd(carina_location)
  system(paste0('git pull'))
  noaa_location <- paste0(data_location,'/','noaa-data')
  setwd(noaa_location)
  system(paste0('git pull'))
  
  #download.file('https://github.com/CareyLabVT/SCCData/raw/carina-data/FCRmet.csv',paste0(workingGLM,'/','FCRmet.csv'))
  #download.file('https://github.com/CareyLabVT/SCCData/raw/mia-data/Catwalk.csv',paste0(workingGLM,'/','Catwalk.csv'))
  #download.file(paste0('https://github.com/CareyLabVT/SCCData/raw/noaa-data/',forecast_base_name,'.csv'),paste0(workingGLM,'/',forecast_base_name,'.csv'))
  
  met_obs_fname <- paste0(carina_location,'/FCRmet.csv')
  ###CREATE HISTORICAL MET FILE
  obs_met_outfile <- paste0(workingGLM,'/','GLM_met.csv')
  create_obs_met_input(fname = met_obs_fname,outfile=obs_met_outfile,full_time_hour_obs, input_tz = 'EST5EDT', output_tz = reference_tzone)
  
  ###CREATE FUTURE MET FILES
  if(forecast_days >0 ){
    in_directory <- paste0(noaa_location)
    out_directory <- workingGLM
    file_name <- forecast_base_name
    #NEED TO DOUBLE CHECK THE INPUT_TZ AND WHY IT IS EST
    process_GEFS2GLM(in_directory,out_directory,file_name, input_tz = 'EST5EDT', output_tz = reference_tzone)
    met_file_names <- rep(NA,nMETmembers)
    for(i in 1:nMETmembers){
      met_file_names[i] <- paste0(met_base_file_name,i,'.csv')
    }
  }
  
  ###MOVE FILES AROUND
  SimFilesFolder <- paste0(Folder,'/','sim_files')
  GLM_folder <- paste0(Folder,'/','glm','/',machine) 
  fl <- c(list.files(SimFilesFolder, full.names = TRUE))
  tmp <- file.copy(from = fl, to = workingGLM,overwrite = TRUE)
  fl <- c(list.files(GLM_folder, full.names = TRUE))
  tmp <- file.copy(from = fl, to = workingGLM,overwrite = TRUE)
  if(!is.na(restart_file)){
    tmp <- file.copy(from = restart_file, to = workingGLM,overwrite = TRUE)
  }
  if(PRE_SCC){
    fl <- c(list.files('/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SCC_data/preSCC/', full.names = TRUE))
    tmp <- file.copy(from = fl, to = workingGLM,overwrite = TRUE)
  }
  
  ##CREATE INFLOW AND OUTFILE FILES
  if(!PRE_SCC){
    create_inflow_outflow_file(full_time = full_time_day,workingGLM = workingGLM, input_tz = 'EST5EDT', output_tz = reference_tzone )
  }
  
  if(include_wq){
    file.copy(from = paste0(workingGLM,'/','glm3_wAED.nml'), to = paste0(workingGLM,'/','glm3.nml'),overwrite = TRUE)
  }else{
    if(!PRE_SCC){
      file.copy(from = paste0(workingGLM,'/','glm3_woAED.nml'), to = paste0(workingGLM,'/','glm3.nml'),overwrite = TRUE)
    }else{
      file.copy(from = paste0(workingGLM,'/','glm3_woAED_preSCC.nml'), to = paste0(workingGLM,'/','glm3.nml'),overwrite = TRUE)
    }
  }
  
  ###SET UP RUN
  
  #DEFINE DEPTHS THAT ARE MODELED
  the_depths_init <- c(0.1, 0.33, 0.66, 1.00, 1.33,1.66,2.00,2.33,2.66,3.0,3.33,3.66,4.0,4.33,4.66,5.0,5.33,5.66,6.0,6.33,6.66,7.00,7.33,7.66,8.0,8.33,8.66,9.00,9.33)
  nlayers_init <- length(the_depths_init)
  
  wq_names <- c('OXY_oxy',
                'CAR_pH','CAR_dic','CAR_ch4', 
                'SIL_rsi',
                'NIT_amm', 'NIT_nit',
                'PHS_frp',
                'OGM_doc','OGM_poc','OGM_don','OGM_pon','OGM_dop','OGM_pop',  #'OGM_docr', 'OGM_donr', 'OGM_dopr','OGM_cpom', 
                'PHY_CYANOPCH1','PHY_CYANONPCH2','PHY_CHLOROPCH3','PHY_DIATOMPCH4')
  num_wq_vars <- length(wq_names) 
  glm_output_vars <- c('temp',wq_names)
  
  
  #Initial States
  the_sals_init <- 0.0
  OXY_oxy_init <- 300.62
  CAR_pH_init <- 6.5
  CAR_dic_init <- 59.1
  CAR_ch4_init <- 0.58
  SIL_rsi_init <- 300
  NIT_amm_init <- 0.69
  NIT_nit_init <- 0.05
  PHS_frp_init <- 0.07
  OGM_doc_init <- 47.4
  OGM_poc_init <- 78.5
  OGM_don_init <- 1.3
  OGM_pon_init <- 8.3
  OGM_dop_init <- 1.5
  OGM_pop_init <- 8.3
  OGM_docr_init <- 350.00
  OGM_donr_init <- 13.0
  OGM_dopr_init <- 3.0
  OGM_cpom_init <- 100.00
  PHY_CYANOPCH1_init <- 2.0
  PHY_CYANONPCH2_init <-2.0
  PHY_CHLOROPCH3_init <-2.0
  PHY_DIATOMPCH4_init <- 2.0
  
  if(!PRE_SCC){
    catwalk_fname <- paste0(mia_location,'/','Catwalk.csv')
    #PROCESS TEMPERATURE OBSERVATIONS
    obs_temp <- extract_temp_chain(fname = catwalk_fname,full_time,input_tz = 'EST5EDT', output_tz = reference_tzone)
    for(i in 1:length(obs_temp$obs[,1])){
      for(j in 1:length(obs_temp$obs[1,])){
        if(obs_temp$obs[i,j] == 0 | is.na(obs_temp$obs[i,j]) | is.nan(obs_temp$obs[i,j])){
          obs_temp$obs[i,j] = NA
        } 
      }
    }
    TempObservedDepths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8,9)
    init_temps1 <- obs_temp$obs[1,]
    
    #PROCESS DO OBSERVATIONS
    DoObservedDepths <- c(1,5,9)
    obs_do <- extract_do_chain(fname = catwalk_fname,full_time,input_tz = 'EST5EDT', output_tz = reference_tzone)
    #mg/L (obs) -> mol/m3 * 31.25
    
  }else{
    fname <- paste0('/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SCC_data/preSCC/CTD_Meta_13_17.csv')
    obs_temp <- extract_temp_CTD(fname,full_time_day,depths = the_depths_init,input_tz = 'EST5EDT', output_tz = reference_tzone)
    TempObservedDepths <- the_depths_init
    init_temps1 <- obs_temp$obs[1,]
  }
  
  
  
  #NEED AN ERROR CHECK FOR WHETHER THERE ARE OBSERVED DATA
  if(is.na(restart_file)){
    if((length(which(init_temps1 != 0.0)) == 0) | length(which(is.na(init_temps1))) >0){
      print('Pick another start day or provide an initial condition file: observations not avialable for starting day')
      break
    }
    temp_inter <- approxfun(TempObservedDepths,init_temps1,rule=2)
    the_temps_init <- temp_inter(the_depths_init)
  }
  
  #SET UP INITIAL CONDITIONS
  
  temp_start <- 1
  temp_end <- length(the_depths_init)
  if(num_pars > 0){
    par1 <- temp_end + 1
    par2 <- par1 + 1
    par3 <-  par2 + 1
  }else{
    par1 <- temp_end
    par2 <- temp_end
    par3 <- temp_end
  }
  wq_start <- rep(NA,num_wq_vars)
  wq_end <- rep(NA,num_wq_vars)
  for(wq in 1:num_wq_vars){
    if(wq == 1){
      wq_start[wq] <- par1+1
      wq_end[wq] <- par1 + (length(the_depths_init))
    }else{
      wq_start[wq] <- wq_end[wq-1]+1
      wq_end[wq] <- wq_end[wq-1] + (length(the_depths_init))
    }
  }
  
  #UPDATE NML WITH PARAMETERS AND INITIAL CONDITIONS
  OXY_oxy_init_depth <- rep(OXY_oxy_init,nlayers_init)
  CAR_pH_init_depth <- rep(CAR_pH_init,nlayers_init)
  CAR_dic_init_depth <- rep(CAR_dic_init,nlayers_init)
  CAR_ch4_init_depth <- rep(CAR_ch4_init,nlayers_init)
  SIL_rsi_init_depth <- rep(SIL_rsi_init,nlayers_init)
  NIT_amm_init_depth <- rep(NIT_amm_init,nlayers_init)
  NIT_nit_init_depth <- rep(NIT_nit_init,nlayers_init)
  PHS_frp_init_depth <- rep(PHS_frp_init,nlayers_init)
  OGM_doc_init_depth <- rep(OGM_doc_init,nlayers_init)
  OGM_poc_init_depth <- rep(OGM_poc_init,nlayers_init)
  OGM_don_init_depth <- rep(OGM_don_init,nlayers_init)
  OGM_pon_init_depth <- rep(OGM_pon_init,nlayers_init)
  OGM_dop_init_depth <- rep(OGM_dop_init,nlayers_init)
  OGM_pop_init_depth <- rep(OGM_pop_init,nlayers_init)
  #OGM_docr_init_depth <- rep(OGM_docr_init,nlayers_init)
  #OGM_donr_init_depth <- rep(OGM_donr_init,nlayers_init)
  #OGM_dopr_init_depth <- rep(OGM_dopr_init,nlayers_init)
  #OGM_cpom_init_depth <- rep(OGM_cpom_init,nlayers_init)
  PHY_CYANOPCH1_init_depth <- rep(PHY_CYANOPCH1_init,nlayers_init)
  PHY_CYANONPCH2_init_depth <- rep(PHY_CYANONPCH2_init,nlayers_init)
  PHY_CHLOROPCH3_init_depth <- rep(PHY_CHLOROPCH3_init,nlayers_init)
  PHY_DIATOMPCH4_init_depth <- rep(PHY_DIATOMPCH4_init,nlayers_init)
  
  wq_init_vals <- c(OXY_oxy_init_depth,
                    CAR_pH_init_depth,
                    CAR_dic_init_depth,
                    CAR_ch4_init_depth,
                    SIL_rsi_init_depth,
                    NIT_amm_init_depth,
                    NIT_nit_init_depth,
                    PHS_frp_init_depth,
                    OGM_doc_init_depth,
                    OGM_poc_init_depth,
                    OGM_don_init_depth,
                    OGM_pon_init_depth,
                    OGM_dop_init_depth,
                    OGM_pop_init_depth,
                    #OGM_docr_init_depth,
                    #OGM_donr_init_depth,
                    #OGM_dopr_init_depth,
                    #OGM_cpom_init_depth,
                    PHY_CYANOPCH1_init_depth,
                    PHY_CYANONPCH2_init_depth,
                    PHY_CHLOROPCH3_init_depth,
                    PHY_DIATOMPCH4_init_depth)
  
  #UPDATE NML WITH PARAMETERS AND INITIAL CONDITIONS
  update_var(wq_init_vals,'wq_init_vals',workingGLM)
  update_var(num_wq_vars,'num_wq_vars',workingGLM)
  update_var(nlayers_init,'num_depths',workingGLM)
  update_var(the_depths_init,'the_depths',workingGLM)
  
  update_var(rep(the_sals_init,nlayers_init),'the_sals',workingGLM)
  
  update_var(sw_factor,'sw_factor',workingGLM)
  update_var(lw_factor,'lw_factor',workingGLM)
  
  #NUMBER OF STATE SIMULATED = SPECIFIED DEPTHS
  if(include_wq){
    nstates <- nlayers_init*(1+num_wq) + num_pars
  }else{
    nstates <- nlayers_init + num_pars
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
  
  if(USE_QT_MATRIX){
    Qt <- read.csv(paste0(workingGLM,'/','Qt_cov_matrix.csv'))
  }
  
  if(num_pars >0){
    Qt <- rbind(Qt,rep(0.0,ncol(Qt)))
    Qt <- cbind(Qt,rep(0.0,nrow(Qt)))  
    Qt[ncol(Qt),nrow(Qt)] <- zone1temp_init_Qt
    Qt <- rbind(Qt,rep(0.0,ncol(Qt)))
    Qt <- cbind(Qt,rep(0.0,nrow(Qt))) 
    Qt[ncol(Qt),nrow(Qt)] <- zone2temp_init_Qt
    Qt <- rbind(Qt,rep(0.0,ncol(Qt)))
    Qt <- cbind(Qt,rep(0.0,nrow(Qt))) 
    Qt[ncol(Qt),nrow(Qt)] <- kw_init_Qt
  }
  
  #NEED TO FIX
  psi <- rep(obs_error,length(obs_index))
  
  ### INITILIZE FIRST TIME STEP
  restart_present <- FALSE
  if(!is.na(restart_file)){
    if(file.exists(restart_file)){
      restart_present <- TRUE
    }
  }
  
  x <- array(NA,dim=c(nsteps,nmembers,nstates))
  x_prior <- array(NA,dim=c(nsteps,nmembers,nstates))
  
  #Initial conditions
  if(!restart_present){
    if(include_wq){
      x <- array(NA,dim=c(nsteps,nmembers,nstates))
      x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init,do_init), sigma=as.matrix(Qt))
      if(NO_UNCERT){
        for(m in 1:nmembers){
          x[1,m,] <- c(the_temps_init,do_init)
        }
      }
    }else{
      if(USE_QT_MATRIX){
        if(num_pars){
          x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init,zone1_temp,zone2_temp,kw_init), sigma=as.matrix(Qt))
        }else{
          x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init), sigma=as.matrix(Qt))
          #x[1,,par1] <- rlnorm(n=nmembers, mean = kw_init, sigma = 0.005)
        }
      }else{
        for(m in 1:nmembers){
          corr_temps <- the_temps_init + rnorm(1,0,temp_error) #tmp(the_depths_init)
          corr_depths <- the_depths_init + rnorm(1,0,thermo_depth_error)
          corrupt_profile <- approxfun(corr_depths,corr_temps,rule = 2)
          x[1,m,temp_start:temp_end] <- corrupt_profile(the_depths_init)
        }
      }
    }
    if(NO_UNCERT){
      for(m in 1:nmembers){
        if(num_pars){
          x[1,m,] <- c(the_temps_init,zone1_temp,zone2_temp,kw_init)
        }else{
          x[1,m,] <- c(the_temps_init)
          #x[1,,par1] <- rlnorm(n=nmembers, mean = kw_init, sigma = 0.005)
        }
      }
    }
    if(!restart_present){
      write.csv(x[1,,],paste0(workingGLM,'/','restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'),row.names = FALSE)
    }
  }
  
  #THIS ALLOWS THE EnKF TO BE RESTARTED FROM YESTERDAY'S RUN
  if(restart_present){
    print('Using restart file')
    nc <- nc_open(restart_file)
    x_previous <- ncvar_get(nc, "x_restart")
    Qt <- ncvar_get(nc, "Qt_restart")
  }else{
    x_previous <- read.csv(paste0(workingGLM,'/','restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'))
  }
  
  if(dim(x[1,,])[1] != dim(x_previous)[1] | dim(x[1,,])[2] != dim(x_previous)[2]){
    print('ERROR: Dimension of the restart file are not correct. State variables or ensemble mismatch?')
    print('Need to fix if the states match but the ensemble number does not')
  }
  
  #Set initial conditions
  x[1,,] <- as.matrix(x_previous)
  x_prior[1,,] <- as.matrix(x_previous)
  #parameter_matrix <- array(NA,dim=c(nmembers,3))
  #parameter_matrix[,1] <- rnorm(nmembers,sw_factor,0.05)
  #parameter_matrix[,2] <- rnorm(nmembers,1.0,0.2)
  #parameter_matrix[,3] <- rnorm(nmembers,0.0013,0.0003)
  
  #Matrix to store ensemble specific deviations and innovations
  dit <- array(NA,dim=c(nmembers,nstates))
  dit_star = array(NA,dim=c(nmembers,nstates)) #Adaptive noise estimation
  
  surface_height <- array(NA,dim=c(nsteps,nmembers))
  surface_height[1,] <- lake_depth_init
  
  file.copy(from = paste0(workingGLM,'/','glm3.nml'), to = paste0(workingGLM,'/','glm3_initial.nml'),overwrite = TRUE)
  
  ###START EnKF
  met_index <- 1
  for(i in 2:nsteps){
    
    #1) Update GLM NML files to match the current day of the simulation
    curr_start <- (full_time[i-1])
    curr_stop <- (full_time[i])
    update_time(start_value  = curr_start, stop_value = curr_stop,workingGLM)
    setwd(workingGLM)
    
    #Create array to hold GLM predictions for each ensemble
    x_star <- array(NA, dim = c(nmembers,nstates))
    x_corr <- array(NA, dim = c(nmembers,nstates))
    for(m in 1:nmembers){
      
      tmp <- update_temps(curr_temps = x[i-1,m,temp_start:temp_end],the_depths_init,workingGLM)
      update_var(surface_height[i-1,m],'lake_depth',workingGLM)
      if(num_pars > 0){
        update_var(c(x[i-1,m,par1],x[i-1,m,par2]),'sed_temp_mean',workingGLM)
        update_var(x[i-1,m,par3],'Kw',workingGLM)
      }
      
      if(include_wq){
        wq_init_vals <- c(x[i-1,wq_start[1]:wq_end[num_wq_vars]])
        update_var(wq_init_vals,'wq_init_vals',workingGLM)
      }
      
      
      #ALLOWS THE LOOPING THROUGH NOAA ENSEMBLES
      if(!PRE_SCC){
        if(i > (hist_days+1)){
          update_var(met_file_names[met_index],'meteo_fl',workingGLM)
          update_var(paste0('FCR_inflow.csv'),'inflow_fl',workingGLM)
          update_var(paste0('FCR_spillway_outflow.csv'),'outflow_fl',workingGLM)
        }else{
          update_var(obs_met_outfile,'meteo_fl',workingGLM)
          update_var(paste0('FCR_inflow.csv'),'inflow_fl',workingGLM)
          update_var(paste0('FCR_spillway_outflow.csv'),'outflow_fl',workingGLM)
        }
      }
      
      #3) Use GLM NML files to run GLM for a day
      # Only allow simulations without NaN values in the output to proceed.  Necessary due to random
      # Nan in AED output
      pass <- FALSE
      num_reruns <- 0
      
      while(!pass){
        unlink(paste0(workingGLM,'/output.nc')) 
        system(paste0(workingGLM,'/','glm'))
        
        if(file.exists(paste0(workingGLM,'/output.nc')) & !has_error(nc_open('output.nc'))){
          if(include_wq){
            GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = the_depths_init,vars = glm_output_vars)
            x_star[m,1:(nstates-num_pars)] <- c(GLM_temp_wq_out$output)
          }else{
            GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = the_depths_init,vars = 'temp')
            x_star[m,temp_start:temp_end] <- c(GLM_temp_wq_out$output)
          }
          x_star[m,par1] <- x[i-1,m,par1]
          x_star[m,par2] <- x[i-1,m,par2]
          x_star[m,par3] <- x[i-1,m,par3]
          surface_height[i,m] <- GLM_temp_wq_out$surface_height 
          if(length(which(is.na(x_star[m,])))==0){
            pass = TRUE
          }else{
            num_reruns <- num_reruns + 1
          }
        }else{
          num_reruns <- num_reruns + 1
        }
        
        if(num_reruns > 1000){
          stop(paste0('Too many re-runs (> 1000) due to NaN values in output'))
        }
      }
      
      #INCREMENT THE MET_INDEX TO MOVE TO THE NEXT NOAA ENSEMBLE
      met_index = met_index + 1
      if(met_index > nMETmembers){
        met_index <- 1
      }
      
      #if(!USE_QT_MATRIX){
      #  corr_temps <- x_star[m,temp_start:temp_end] + rnorm(1,0,temp_error) #tmp(the_depths_init)
      #  corr_depths <- the_depths_init + rnorm(1,0,thermo_depth_error)
      #  corrupt_profile <- approxfun(corr_depths,corr_temps,rule = 2)
      #  corrupt_profile <- approxfun(corr_depths,corr_temps,rule = 2)
      #  x_corr[m,temp_start:temp_end] <- corrupt_profile(the_depths_init)
      #}
      
    }
    
    #DEAL WITH ENSEMBLE MEMBERS THAT ARE 'BAD' AND PRODUCE NA VALUES OR HAVE NEGATIVE TEMPERATURES
    # THIS RANDOMLY REPLACES IT WITH A GOOD ENSEMBLE MEMBER
    if(length(which(is.na(c(x_star))))>0){
      good_index <- NULL
      for(m in 1:nmembers){
        if(length(which(is.na(c(x_star[m,])))) == 0 & length(which(c(x_star[m,]) <= 0)) == 0){
          good_index <- c(good_index,m)
        }
      }
      for(m in 1:nmembers){
        if(length(which(is.na(c(x_star[m,])))) > 0 | length(which(c(x_star[m,]) <= 0) > 0)){
          replace_index <- sample(good_index,1)
          x_star[m,] <- x_star[replace_index,]
          surface_height[i,m] <- surface_height[i,replace_index]
        }
      }
    }
    
    #Corruption [nmembers x nstates] 
    if(USE_QT_MATRIX){
      NQt <- rmvnorm(n=nmembers, sigma=as.matrix(Qt))
      zero_par3 <- which(x_star[,par3] + NQt[,par3] <= 0.0)
      if(length(zero_par3) > 0.0){
        for(kk in 1:length(zero_par3)){
          pass <- FALSE
          while(!pass){
            NQt[kk,] <- rmvnorm(n=1, sigma=as.matrix(Qt))
            if(NQt[kk,par3] > 0.0){pass <- TRUE}
          }
        }
      }
      x_corr <- x_star + NQt
    }
    
    x_prior[i,,] <- x_corr
    
    if(i >= spin_up_days+1){
      #Obs for time step
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
        
        N_psi = t(rmvnorm(n=1,mean = zt, sigma=as.matrix(psi_t)))
        D_mat <- t(matrix(rep(N_psi,each=nmembers), nrow = nmembers, ncol=length(N_psi)))
        
        #Loop through ensemble members
        for(m in 1:nmembers){  
          #  
          #  #Ensemble specific deviation
          dit[m,] <- x_corr[m,]-ens_mean
          #dit_star[m,] <- x_star[m,] - ens_mean_star #Adaptive noise estimation
          
          #  #Ensemble specific estimate and innovation covariance
          if(m == 1){
            Pit <- dit[m,] %*% t(dit[m,]) 
            #    Pit_star <- dit_star[m,] %*% t(dit_star[m,]) 
          }else{
            Pit <- dit[m,] %*% t(dit[m,]) +  Pit 
            #     Pit_star <- dit_star[m,] %*% t(dit_star[m,]) +  Pit_star 
          }
        }
        
        #estimate covariance
        Pt <- Pit/nmembers
        
        #Kalman gain
        Kt <- Pt %*% t(H) %*% solve(H%*%Pt%*%t(H)+psi_t)
        
        #Update states array (transposes are necessary to convert between the dims here and the dims in the EnKF formulations)
        x[i,,] <- t(t(x_corr) + Kt%*%(D_mat - H%*%t(x_corr)))
        
        if(length(which(is.na(x[i,,]))) > 0){dies = i}
      }
    }else{
      x[i,,] <- x_star
    }
    
    if(i == (hist_days+1)){
      x_restart <- x[i,,]
      Qt_restart <- Qt
    }
    
  }
  
  if(forecast_days >0){
    save_file_name <- paste0(sim_name,'_hist_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_forecast_',
                             year(full_time[hist_days+1]),'_',month(full_time[hist_days+1]),'_',day(full_time[hist_days+1]))
  }else{
    save_file_name <- paste0(sim_name,'_hist_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]))    
  }
  
  ### SUMMARIZE FORECAST
  time_of_forecast <- Sys.time() #paste0(year(Sys.time()),month(Sys.time()),day(Sys.time()),'_',hour(Sys.time()),'_',(minute(Sys.time())))
  time_of_forecast_string <- paste0(year(Sys.time()),month(Sys.time()),day(Sys.time()),'_',hour(Sys.time()),'_',(minute(Sys.time())))
  
  ###SAVE FORECAST
  write_forecast_netcdf(x = x,
                        full_time = full_time_local,
                        Qt = Qt,
                        the_depths_init = the_depths_init,
                        save_file_name = save_file_name,
                        x_restart=x_restart,
                        Qt_restart = Qt_restart,
                        time_of_forecast = time_of_forecast,
                        hist_days = hist_days,
                        x_prior)
  
  ##ARCHIVE FORECAST
  restart_file_name <- archive_forecast(workingGLM = workingGLM,
                                        Folder = Folder, 
                                        forecast_base_name = forecast_base_name, 
                                        forecast_location = forecast_location,
                                        push_to_git = push_to_git,
                                        save_file_name = save_file_name, 
                                        time_of_forecast = time_of_forecast)
  
  return(list(restart_file_name <- restart_file_name,sim_name <- paste0(save_file_name,'_',time_of_forecast_string)))
}
