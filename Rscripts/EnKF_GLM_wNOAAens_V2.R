run_forecast<-function(start_day= '2018-07-06 00:00:00', 
                       sim_name = NA, 
                       hist_days = 1,
                       forecast_days = 16,  
                       spin_up_days = 0,
                       restart_file = NA,
                       Folder, 
                       forecast_location = NA,
                       push_to_git=FALSE,
                       data_location = NA, 
                       nEnKFmembers = NA,
                       include_wq = FALSE,
                       USE_CTD = USE_CTD,
                       uncert_mode = 1,
                       cov_matrix = NA,
                       alpha = c(0.5,0.5,0.5)){
  
  ###RUN OPTIONS
  npars <- 3
  
  PRE_SCC <- FALSE
  
  if(uncert_mode == 1){
    #All sources of uncertainity and data used to constrain 
    USE_OBS_CONTRAINT <- TRUE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- TRUE
    WEATHER_UNCERTAINITY <- TRUE
    INITIAL_CONDITION_UNCERTAINITY <- TRUE
    PARAMETER_UNCERTAINITY <- TRUE
  }else if(uncert_mode == 2){
    #No sources of uncertainity and no data used to constrain 
    USE_OBS_CONTRAINT <- TRUE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- FALSE
    WEATHER_UNCERTAINITY <- FALSE
    INITIAL_CONDITION_UNCERTAINITY <- FALSE
    PARAMETER_UNCERTAINITY <- FALSE
  }else if(uncert_mode == 3){
    #Only process uncertainity
    USE_OBS_CONTRAINT <- TRUE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- TRUE
    WEATHER_UNCERTAINITY <- FALSE
    INITIAL_CONDITION_UNCERTAINITY <- FALSE
    PARAMETER_UNCERTAINITY <- FALSE
  }else if(uncert_mode == 4){
    #only weather uncertainity
    USE_OBS_CONTRAINT <- TRUE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- FALSE
    WEATHER_UNCERTAINITY <- TRUE
    INITIAL_CONDITION_UNCERTAINITY <- FALSE
    PARAMETER_UNCERTAINITY <- FALSE
  }else if(uncert_mode == 5){
    #only initial condition uncertainity with data constraint
    USE_OBS_CONTRAINT <- TRUE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- FALSE
    WEATHER_UNCERTAINITY <- FALSE
    INITIAL_CONDITION_UNCERTAINITY <- TRUE
    PARAMETER_UNCERTAINITY <- FALSE
  }else if(uncert_mode == 6){
    #only initial condition uncertainity without data constraint
    USE_OBS_CONTRAINT <- FALSE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- FALSE
    WEATHER_UNCERTAINITY <- FALSE
    INITIAL_CONDITION_UNCERTAINITY <- TRUE
    PARAMETER_UNCERTAINITY <- FALSE
  }else if(uncert_mode == 7){
    #only parameter uncertainity
    USE_OBS_CONTRAINT <- FALSE
    #SOURCES OF UNCERTAINITY
    OBSERVATION_UNCERTAINITY <- TRUE
    PROCESS_UNCERTAINITY <- FALSE
    WEATHER_UNCERTAINITY <- FALSE
    INITIAL_CONDITION_UNCERTAINITY <- FALSE
    PARAMETER_UNCERTAINITY <- TRUE
  }
  
  #Parameters
  lake_depth_init <- 9.4  #not a modeled state
  zone2_temp <- 17
  zone1_temp <- 11
  zone1temp_init_Qt <- 0.01 #THIS IS THE VARIANCE, NOT THE SD
  zone2temp_init_Qt <- 0.01 #THIS IS THE VARIANCE, NOT THE SD
  if(include_wq){
    kw_init <- 0.1
    kw_init_Qt <- 0.00000000001
  }else{
    kw_init <- 1.0
    kw_init_Qt <- 0.01^2 #THIS IS THE VARIANCE, NOT THE SD
  }
  
  #ERROR TERMS
  if(OBSERVATION_UNCERTAINITY){
    obs_error <- 0.0001 #NEED TO DOUBLE CHECK
  }else{
    obs_error <- 0.0001 #NEED TO DOUBLE CHECK   
  }
  
  ####################################################
  #### YOU WON'T NEED TO MODIFY ANYTHING BELOW HERE ##
  ####################################################
  
  # SET UP NUMBER OF ENSEMBLE MEMBERS
  nMETmembers <- 21
  nmembers <- nEnKFmembers*nMETmembers
  
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
  
  begin_sim  <- as.POSIXct(start_day,tz = reference_tzone)
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
  full_time_local <- strftime(full_time_local, format="%Y-%m-%d %H:%M",tz = 'EST5EDT')
  full_time_day <- strftime(full_time, format="%Y-%m-%d",tz = reference_tzone)
  full_time_day_local <- strftime(full_time_local, format="%Y-%m-%d",tz = 'EST5EDT')
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
  
  ###CREATE HISTORICAL MET FILE
  met_obs_fname <- paste0(carina_location,'/FCRmet.csv')
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
  the_depths_init <- c(0.1, 0.33, 0.66, 
                       1.00, 1.33,1.66,
                       2.00,2.33,2.66,
                       3.0,3.33,3.66,
                       4.0,4.33,4.66,
                       5.0,5.33,5.66,
                       6.0,6.33,6.66,
                       7.00,7.33,7.66,
                       8.0,8.33,8.66,
                       9.00,9.33)
  
  nlayers_init <- length(the_depths_init)
  
  #DEFINE WATER QUALITY VARIABLES
  wq_names <- c('OXY_oxy',
                'CAR_pH','CAR_dic','CAR_ch4',
                'SIL_rsi',
                'NIT_amm', 'NIT_nit',
                'PHS_frp',
                'OGM_doc','OGM_poc','OGM_don','OGM_pon','OGM_dop','OGM_pop',
                'PHY_CYANOPCH1','PHY_CYANONPCH2','PHY_CHLOROPCH3','PHY_DIATOMPCH4',
                'ZOO_COPEPODS1','ZOO_DAPHNIABIG2','ZOO_DAPHNIASMALL3')
  
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
  PHY_CYANOPCH1_init <- 2.0
  PHY_CYANONPCH2_init <-2.0
  PHY_CHLOROPCH3_init <-2.0
  PHY_DIATOMPCH4_init <- 2.0
  ZOO_COPEPODS1_init <- 2.9
  ZOO_DAPHNIABIG2_init <- 4.3
  ZOO_DAPHNIASMALL3_init <- 40
  
  OXY_oxy_error <- OXY_oxy_init*1.0
  CAR_pH_error <- CAR_pH_init*0.001
  CAR_dic_error <- CAR_dic_init*0.001
  CAR_ch4_error <- CAR_ch4_init*0.001
  SIL_rsi_error <- SIL_rsi_init*0.001
  NIT_amm_error <- NIT_amm_init*0.001
  NIT_nit_error <- NIT_nit_init*0.001
  PHS_frp_error <- PHS_frp_init*0.001
  OGM_doc_error <- OGM_doc_init*0.001
  OGM_poc_error <- OGM_poc_init*0.001
  OGM_don_error <- OGM_don_init*0.001
  OGM_pon_error <- OGM_pon_init*0.001
  OGM_dop_error <- OGM_dop_init*0.001
  OGM_pop_error <- OGM_pop_init*0.001
  PHY_CYANOPCH1_error <- PHY_CYANOPCH1_init*0.01
  PHY_CYANONPCH2_error <-PHY_CYANONPCH2_init*0.01
  PHY_CHLOROPCH3_error <-PHY_CHLOROPCH3_init*0.01
  PHY_DIATOMPCH4_error <- PHY_DIATOMPCH4_init*0.01
  ZOO_COPEPODS1_error <- ZOO_COPEPODS1_init*0.001
  ZOO_DAPHNIABIG2_error <- ZOO_DAPHNIABIG2_init*0.001
  ZOO_DAPHNIASMALL3_error <- ZOO_DAPHNIASMALL3_init*0.001
  
  wq_var_error <- c(OXY_oxy_error,
                    CAR_pH_error,CAR_dic_error,CAR_ch4_error,
                    SIL_rsi_error,
                    NIT_amm_error,NIT_nit_error,
                    PHS_frp_error,
                    OGM_doc_error,OGM_poc_error, OGM_don_error ,OGM_pon_error,OGM_dop_error,OGM_pop_error,
                    PHY_CYANOPCH1_error,PHY_CYANONPCH2_error,PHY_CHLOROPCH3_error,PHY_DIATOMPCH4_error,
                    ZOO_COPEPODS1_error,ZOO_DAPHNIABIG2_error,ZOO_DAPHNIASMALL3_error)
  
  
  #Extract observations
  catwalk_fname <- paste0(mia_location,'/','Catwalk.csv')
  #PROCESS TEMPERATURE OBSERVATIONS
  TempObservedDepths <- c(0.1, 1, 2, 3, 4, 5, 6, 7, 8,9)
  obs_temp <- extract_temp_chain(fname = catwalk_fname,full_time,depths = the_depths_init,TempObservedDepths = TempObservedDepths,
                                 input_tz = 'EST5EDT', output_tz = reference_tzone)
  for(i in 1:length(obs_temp$obs[,1])){
    for(j in 1:length(obs_temp$obs[1,])){
      if(obs_temp$obs[i,j] == 0 | is.na(obs_temp$obs[i,j]) | is.nan(obs_temp$obs[i,j])){
        obs_temp$obs[i,j] = NA
      } 
    }
  }
  
  init_temps <- obs_temp$obs[1,]
  
  #PROCESS DO OBSERVATIONS
  DoObservedDepths <- c(1,5,9)
  obs_do <- extract_do_chain(fname = catwalk_fname,full_time,depths = the_depths_init,DoObservedDepths= DoObservedDepths,
                             input_tz = 'EST5EDT', output_tz = reference_tzone)
  obs_do$obs <- obs_do$obs*1000/32  #mg/L (obs units) -> mmol/m3 (glm units)
  init_do1 <- obs_do$obs[1,]
  
  Chla_fDOM_ObservedDepths <- 1
  obs_chla_fdom <- extract_chla_chain(fname = catwalk_fname,full_time,depths = the_depths_init,Chla_fDOM_ObservedDepths= Chla_fDOM_ObservedDepths,
                                      input_tz = 'EST5EDT', output_tz = reference_tzone)
  
  #Use the CTD observation rather than the sensor string when CTD data is avialable
  if(USE_CTD){
    ## LOOK AT CTD DATA
    fl <- c(list.files('/Users/quinn/Dropbox (VTFRS)/Research/SSC_forecasting/SCC_data/preSCC/', pattern = 'CTD', full.names = TRUE))
    #NEED TO DOUBLE CHECK TIME ZONE
    obs_ctd <- extract_temp_CTD(fname = fl[1],full_time_day_local,depths = the_depths_init,input_tz = 'EST5EDT', output_tz = reference_tzone)
    obs_ctd$obs_do <- obs_ctd$obs_do*1000/32
    for(i in 1:length(full_time_day)){
      if(!is.na(obs_ctd$obs_temp[i,1])){
        obs_temp$obs[i,] <- obs_ctd$obs_temp[i,]
        obs_do$obs[i,] <- obs_ctd$obs_do[i,]
        obs_chla_fdom$Chla_obs[i,] <- obs_ctd$obs_chla[i,]
      }
    }
    init_pH_obs <- obs_ctd$obs_pH[1,which(!is.na(obs_ctd$obs_pH[1,]))]
    init_obs_pH_depths <- the_depths_init[which(!is.na(obs_ctd$obs_pH[1,]))]
    
    init_sal_obs <- obs_ctd$obs_sal[1,which(!is.na(obs_ctd$obs_sal[1,]))]
    init_obs_sal_depths <- the_depths_init[which(!is.na(obs_ctd$obs_sal[1,]))]
  }
  
  init_temps_obs <- obs_temp$obs[1,which(!is.na(obs_temp$obs[1,]))]
  init_obs_temp_depths <- the_depths_init[which(!is.na(obs_temp$obs[1,]))]
  
  init_do_obs <- obs_do$obs[1,which(!is.na(obs_do$obs[1,]))]
  init_obs_do_depths <- the_depths_init[which(!is.na(obs_do$obs[1,]))]
  
  #NEED AN ERROR CHECK FOR WHETHER THERE ARE OBSERVED DATA
  if(is.na(restart_file)){
    if((length(which(init_temps_obs != 0.0)) == 0) | length(which(is.na(init_temps_obs))) >0){
      print('Pick another start day or provide an initial condition file: observations not avialable for starting day')
      break
    }
    temp_inter <- approxfun(init_obs_temp_depths,init_temps_obs,rule=2)
    the_temps_init <- temp_inter(the_depths_init)
    if(include_wq){
      do_inter <- approxfun(init_obs_do_depths,init_do_obs,rule=2)
      do_init <- do_inter(the_depths_init)
      if(length(which(!is.na(init_pH_obs)))>0){
        pH_inter <- approxfun(init_obs_pH_depths,init_pH_obs,rule=2)
        pH_init <- pH_inter(the_depths_init)
      }
    }
  }
  
  #SET UP INITIAL CONDITIONS
  
  if(include_wq){
    temp_start <- 1
    temp_end <- length(the_depths_init)
    wq_start <- rep(NA,num_wq_vars)
    wq_end <- rep(NA,num_wq_vars)
    for(wq in 1:num_wq_vars){
      if(wq == 1){
        wq_start[wq] <- temp_end+1
        wq_end[wq] <- temp_end + (length(the_depths_init))
      }else{
        wq_start[wq] <- wq_end[wq-1]+1
        wq_end[wq] <- wq_end[wq-1] + (length(the_depths_init))
      }
      
      if(npars > 0){ #NEED TO GENERALIZE
        par1 <- wq_end[num_wq_vars] + 1
        par2 <- par1 + 1
        par3 <-  par2 + 1
      }else{
        par1 <- wq_end[num_wq_vars]
        par2 <- wq_end[num_wq_vars]
        par3 <- wq_end[num_wq_vars]
      }
      
    }
  }else{
    temp_start <- 1
    temp_end <- length(the_depths_init)
    if(npars > 0){
      par1 <- temp_end + 1
      par2 <- par1 + 1
      par3 <-  par2 + 1
    }else{
      par1 <- temp_end
      par2 <- temp_end
      par3 <- temp_end
    }
  }
  
  #UPDATE NML WITH PARAMETERS AND INITIAL CONDITIONS
  if(include_wq){
    OXY_oxy_init_depth <- do_init #rep(OXY_oxy_init,nlayers_init)
  }else{
    OXY_oxy_init_depth <- rep(OXY_oxy_init,nlayers_init)    
  }
  if(include_wq & USE_CTD){
    CAR_pH_init_depth <- pH_init
  }else{
    CAR_pH_init_depth <- rep(CAR_pH_init,nlayers_init) 
  }
  
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
  PHY_CYANOPCH1_init_depth <- rep(PHY_CYANOPCH1_init,nlayers_init)
  PHY_CYANONPCH2_init_depth <- rep(PHY_CYANONPCH2_init,nlayers_init)
  PHY_CHLOROPCH3_init_depth <- rep(PHY_CHLOROPCH3_init,nlayers_init)
  PHY_DIATOMPCH4_init_depth <- rep(PHY_DIATOMPCH4_init,nlayers_init)
  ZOO_COPEPODS1_init_depth <- rep(ZOO_COPEPODS1_init,nlayers_init)
  ZOO_DAPHNIABIG2_init_depth <- rep(ZOO_DAPHNIABIG2_init,nlayers_init)
  ZOO_DAPHNIASMALL3_init_depth <- rep(ZOO_DAPHNIASMALL3_init,nlayers_init)
  
  #if(full_time_day_local[1] == as.data.frame.POSIXct('201')
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8, 9, 10)
  #mg/L
  #curr_values <- c(3.764, 3.781, 3.578, 5.156, 5.2735, 5.5165, 5.222, 5.368)
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #CAR_dic_init_depth <- inter(the_depths_init)
  
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8,9,9.5)
  #umol CH4/L
  #curr_values <- c(3.91E-04,0.370572728,0.107597836,0.126096596,0.088502664,0.086276629,0.07256043,0.07249431)
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #CAR_ch4_init_depth <- inter(the_depths_init)
  
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8, 9, 10)
  #curr_values <- c(12.65291714,4.213596723,10.5935375,13.43611258,11.34765394,11.95676704,11.98577285,12.82695814)
  #ug/L
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #NIT_amm_init_depth <- inter(the_depths_init)
  
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8, 9, 10)
  #curr_values <-c(5.68,3.82,4.46,3.71,4.18,5.08,3.01,7.72)
  #ug/L
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #NIT_nit_init_depth <- inter(the_depths_init)
  
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8, 9, 10)
  #ug/L
  #curr_values <- c(8.96,7.66,6.26,6.22,7.72,9.69,7.95,10.5)
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #PHS_frp_init_depth <- inter(the_depths_init)
  
  #curr_depths <- c(0.1,1.6,3.8,5,6.2,8, 9, 10)
  ##mg/L
  #curr_values <- c(4.2315,4.374, 3.2655,2.9705,2.938,2.922,2.773,2.9525)
  #inter <- approxfun(curr_depths,curr_values,rule=2)
  #OGM_poc_init_depth <- inter(the_depths_init)
  
  #}else{
  
  #}
  
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
                    PHY_CYANOPCH1_init_depth,
                    PHY_CYANONPCH2_init_depth,
                    PHY_CHLOROPCH3_init_depth,
                    PHY_DIATOMPCH4_init_depth,
                    ZOO_COPEPODS1_init_depth,
                    ZOO_DAPHNIABIG2_init_depth,
                    ZOO_DAPHNIASMALL3_init_depth)
  
  #UPDATE NML WITH PARAMETERS AND INITIAL CONDITIONS
  if(include_wq){
    update_var(wq_init_vals,'wq_init_vals',workingGLM)
    update_var(num_wq_vars,'num_wq_vars',workingGLM)
  }else{
    update_var(' ','wq_init_vals',workingGLM)
    update_var(0,'num_wq_vars',workingGLM)
  }
  update_var(nlayers_init,'num_depths',workingGLM)
  update_var(the_depths_init,'the_depths',workingGLM)
  update_var(rep(the_sals_init,nlayers_init),'the_sals',workingGLM)
  
  #NUMBER OF STATE SIMULATED = SPECIFIED DEPTHS
  if(include_wq){
    nstates <- nlayers_init*(1+num_wq_vars)
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
  #z <- t(matrix(rep(NA,nobs), nrow = nobs, ncol = nsteps))
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
    obs_index <- rep(NA,length(the_depths_init)*(num_wq_vars+1))
    obs_index[1:length(the_depths_init)] <- seq(1,length(the_depths_init),1)
    for(wq in 1:num_wq_vars){
      obs_index[wq_start[wq]:wq_end[wq]] <- seq(wq_start[wq],wq_end[wq],1)
    }
  }else{
    obs_index <- rep(NA,length(the_depths_init))
    obs_index[1:length(the_depths_init)] <- seq(1,length(the_depths_init),1)
  }
  
  #Matrix for knowing which state the observation corresponds to
  z_states <- t(matrix(obs_index, nrow = length(obs_index), ncol = nsteps))
  
  #Process error 
  if(is.na(cov_matrix)){
    Qt <- read.csv(paste0(workingGLM,'/','Qt_cov_matrix.csv'))
  }else{
    Qt <- read.csv(paste0(workingGLM,'/',cov_matrix))
  }
  
  if(include_wq){
    for(i in 1:num_wq_vars){
      for(j in 1:nlayers_init){
        Qt <- rbind(Qt,rep(0.0,ncol(Qt)))
        Qt <- cbind(Qt,rep(0.0,nrow(Qt)))
        Qt[ncol(Qt),nrow(Qt)] <- wq_var_error[i]
      }
    }
  }
  
  #Covariance matrix for parameters
  Qt_pars <- matrix(data = 0,nrow = npars, ncol = npars)
  diag(Qt_pars) <- c(zone1temp_init_Qt,zone2temp_init_Qt,kw_init_Qt)
  
  psi <- rep(obs_error,length(obs_index))
  
  ### INITILIZE FIRST TIME STEP
  restart_present <- FALSE
  if(!is.na(restart_file)){
    if(file.exists(restart_file)){
      restart_present <- TRUE
    }
  }
  
  x <- array(NA,dim=c(nsteps,nmembers,nstates + npars))
  x_prior <- array(NA,dim=c(nsteps,nmembers,nstates + npars))
  
  #Initial conditions
  if(!restart_present){
    if(include_wq){
      if(npars > 0){
        x[1,,1:nstates] <- rmvnorm(n=nmembers, mean=c(the_temps_init,wq_init_vals), sigma=as.matrix(Qt))
        x[1,,nstates+1:nstates+npars] <- rmvnorm(n=nmembers, mean=c(zone1_temp,zone2_temp,kw_init),sigma = as.matrix(Qt_pars))
        if(INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            x[1,m,] <- c(the_temps_init,wq_init_vals,zone1_temp,zone2_temp,kw_init)
          }
        }
      }else{
        x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init,wq_init_vals), sigma=as.matrix(Qt))
        if(INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            x[1,m,] <- c(the_temps_init,do_init,wq_init_vals)
          }
        }
      }
    }else{
      if(npars > 0){
        x[1,,1:nstates] <- rmvnorm(n=nmembers, mean=c(the_temps_init), sigma=as.matrix(Qt))
        x[1,,(nstates+1):(nstates+npars)] <- rmvnorm(n=nmembers, mean=c(zone1_temp,zone2_temp,kw_init),sigma = as.matrix(Qt_pars))
        if(INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            x[1,m,] <- c(the_temps_init,zone1_temp,zone2_temp,kw_init)
          }
        }
      }else{
        x[1,,] <- rmvnorm(n=nmembers, mean=c(the_temps_init), sigma=as.matrix(Qt))
        if(INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            if(npars > 0){
              x[1,m,] <- c(the_temps_init)
            }
          }
        }
      }
    }
    if(include_wq){
      for(m in 1:nmembers){
        for(wq in 1:num_wq_vars){
          index <- which(x[1,m,] < 0.0)
          index <- index[which(index > wq_start[1])]
          x[1,m,index] <- 0.0
        }
      }
    }
    write.csv(x[1,,],paste0(workingGLM,'/','restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'),row.names = FALSE)
  }
  
  #THIS ALLOWS THE EnKF TO BE RESTARTED FROM YESTERDAY'S RUN
  if(restart_present){
    print('Using restart file')
    nc <- nc_open(restart_file)
    restart_nmembers <- length(ncvar_get(nc,'ens'))
    if(restart_nmembers > nmembers){
      #sample restart_nmembers
      sampled_nmembers <- sample(seq(1,restart_nmembers,1),nmembers,replace=FALSE)
      restart_x_previous <- ncvar_get(nc, "x_restart")
      x_previous <- restart_x_previous[sampled_nmembers,]
      if(INITIAL_CONDITION_UNCERTAINITY == FALSE & hist_days == 0){
        x_previous_1 <- colMeans(x_previous)
        for(m in 1:nmembers){
          x_previous[m,] <- x_previous_1
        }
      }
    }else if(restart_nmembers < nmembers){
      sampled_nmembers <- sample(seq(1,restart_nmembers,1),nmembers,replace=TRUE)
      restart_x_previous <- ncvar_get(nc, "x_restart")
      x_previous <- restart_x_previous[sampled_nmembers,]
      if(INITIAL_CONDITION_UNCERTAINITY == FALSE & hist_days == 0){
        x_previous_1 <- colMeans(x_previous)
        for(m in 1:nmembers){
          x_previous[m,] <- x_previous_1
        }
      }
    }else{
      x_previous <- ncvar_get(nc, "x_restart")   
      if(INITIAL_CONDITION_UNCERTAINITY == FALSE & hist_days == 0){
        x_previous_1 <- colMeans(x_previous) 
        for(m in 1:nmembers){
          x_previous[m,] <- x_previous_1
        }
      }
    }
    #Qt <- ncvar_get(nc, "Qt_restart")
    nc_close(nc)
  }else{
    x_previous <- read.csv(paste0(workingGLM,'/','restart_',year(full_time[1]),'_',month(full_time[1]),'_',day(full_time[1]),'_cold.csv'))
  }
  
  #Set initial conditions
  x[1,,] <- as.matrix(x_previous)
  x_prior[1,,] <- as.matrix(x_previous)
  

  
  #Matrix to store essemble specific surface height
  surface_height <- array(NA,dim=c(nsteps,nmembers))
  surface_height[1,] <- lake_depth_init
  
  #Create a copy of the NML to record starting parameters
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
    pars_corr <-  array(NA, dim = c(nmembers,npars))
    #Matrix to store calculated ensemble specific deviations and innovations
    dit <- array(NA,dim=c(nmembers,nstates))
    dit_pars<- array(NA,dim=c(nmembers,npars))
    
    for(m in 1:nmembers){
      
      tmp <- update_temps(curr_temps = round(x[i-1,m,temp_start:temp_end],3),the_depths_init,workingGLM)
      update_var(surface_height[i-1,m],'lake_depth',workingGLM)
      if(npars > 0){
        if(i > (hist_days+1)){
          new_pars <- x[i-1,m,(nstates+1):(nstates+npars)]
          if(PARAMETER_UNCERTAINITY == FALSE){
            new_pars <- mean(x[i-1,,(nstates+1):(nstates+npars)])
          }
        }else{
          new_pars <- rmvnorm(n=1, mean = c(x[i-1,m,(nstates+1):(nstates+npars)]),sigma=as.matrix(Qt_pars))
        }
        
        update_var(c(round(new_pars[1],3),round(new_pars[2],3)),'sed_temp_mean',workingGLM)
        update_var(round(new_pars[3],3),'sw_factor',workingGLM)
        update_var(round(new_pars[3],3),'lw_factor',workingGLM)
        pars_corr[m,] <- new_pars
      }
      
      if(include_wq){
        wq_init_vals <- round(c(x[i-1,m,wq_start[1]:wq_end[num_wq_vars]]),3)
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
      
      #Use GLM NML files to run GLM for a day
      # Only allow simulations without NaN values in the output to proceed.  Necessary due to random
      # Nan in AED output
      pass <- FALSE
      num_reruns <- 0
      
      while(!pass){
        unlink(paste0(workingGLM,'/output.nc')) 
        system(paste0(workingGLM,'/','glm'))
        
        if(file.exists(paste0(workingGLM,'/output.nc')) & !has_error(nc<-nc_open('output.nc'))){
          if(length(ncvar_get(nc,'time')) > 1){
            nc_close(nc)
            if(include_wq){
              GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = the_depths_init,vars = glm_output_vars)
              x_star[m,1:(nstates-npars)] <- c(GLM_temp_wq_out$output)
            }else{
              GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = the_depths_init,vars = 'temp')
              x_star[m,temp_start:temp_end] <- c(GLM_temp_wq_out$output)
            }
            
            surface_height[i,m] <- GLM_temp_wq_out$surface_height 
            if(length(which(is.na(x_star[m,])))==0){
              pass = TRUE
            }else{
              num_reruns <- num_reruns + 1
            }
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
      if(met_index > nMETmembers | WEATHER_UNCERTAINITY == FALSE){
        met_index <- 1
      }
      
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
    NQt <- rmvnorm(n=nmembers, sigma=as.matrix(Qt))
    x_corr <- x_star + NQt
    
    if(include_wq){
      for(m in 1:nmembers){
        for(wq in 1:num_wq_vars){
          index <- which(x_corr[m,] < 0.0)
          index <- index[which(index > wq_start[1])]
          x_corr[m,index] <- 0.0
        }
      }
    }
    
    x_prior[i,,] <- cbind(x_corr,pars_corr)
    
    if(i >= spin_up_days+1){
      #Obs for time step
      z_index <- which(!is.na(z[i,]))
      
      #if no observations at a time step then just propogate model uncertainity
      if(length(z_index) == 0 | i > (hist_days+1)){
        x[i,,] <- cbind(x_corr,pars_corr) 
        if(PROCESS_UNCERTAINITY == FALSE & i > (hist_days+1)){
          x[i,,] <- cbind(x_star,pars_corr)
        }
        if(i == (hist_days+1) & INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            x[i,m,] <- colMeans(cbind(x_star,pars_corr)) 
          }
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
        if(npars > 0){
          par_mean <- apply(pars_corr, 2, mean)
        }
        
        N_psi = t(rmvnorm(n=1,mean = zt, sigma=as.matrix(psi_t)))
        D_mat <- t(matrix(rep(N_psi,each=nmembers), nrow = nmembers, ncol=length(N_psi)))
        
        #Loop through ensemble members
        for(m in 1:nmembers){  
          
          #  #Ensemble specific deviation
          dit[m,] <- x_corr[m,]-ens_mean
          dit_pars[m, ] <- pars_corr[m,] - par_mean
          
          if(m == 1){
            Pit <- dit[m,] %*% t(dit[m,]) 
            Pit_pars <- dit_pars[m, ] %*% t(dit[m,])
          }else{
            Pit <- dit[m,] %*% t(dit[m,]) +  Pit 
            Pit_pars <- dit_pars[m, ] %*% t(dit[m,]) + Pit_pars
          }
        }
        
        #estimate covariance
        Pt <- Pit/(nmembers-1)
        Pt_pars <- Pit_pars/(nmembers-1)
        #Kalman gain
        Kt <- Pt %*% t(H) %*% solve(H%*%Pt%*%t(H)+psi_t)
        Kt_pars <- Pt_pars %*% t(H) %*% solve(H%*%Pt%*%t(H)+psi_t)
        
        #Update states array (transposes are necessary to convert between the dims here and the dims in the EnKF formulations)
        x[i,,1:nstates] <- t(t(x_corr) + Kt%*%(D_mat - H%*%t(x_corr)))
        for(pp in 1:npars){
          x[i,,(nstates+pp)] <- alpha[pp]*x[i-1,,(nstates+pp)] + (1-alpha[pp])*t(t(pars_corr) + Kt_pars%*%(D_mat - H%*%t(x_corr)))[,pp]
        }
        
        if(include_wq){
          for(m in 1:nmembers){
            for(wq in 1:num_wq_vars){
              index <- which(x[i,m,] < 0.0)
              index <- index[which(index > wq_start[1])]
              x[i,m,index] <- 0.0
            }
          }
        }
        
        #IF NO INITIAL CONDITION UNCERTAINITY THEN SET EACH ENSEMBLE MEMBER TO THE MEAN
        #AT THE INITIATION OF THE FUTURE FORECAST
        if(i == (hist_days+1) & INITIAL_CONDITION_UNCERTAINITY == FALSE){
          for(m in 1:nmembers){
            if(PARAMETER_UNCERTAINITY == FALSE){
            x[i,m,] <- colMeans(cbind(x_star,pars_corr)) 
            }else{
              x[i,m,] <- cbind(colMeans(x_star),x[i,m,(nstates+1):(nstates+npars)])
            }
          }
        }
        
        if(length(which(is.na(x[i,,]))) > 0){dies = i}
      }
    }else{
      x[i,,] <- cbind(x_star,pars_corr)
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
                        x_prior,
                        include_wq,
                        wq_start,
                        wq_end,
                        par1,
                        par2,
                        par3,
                        z,
                        nstates,
                        npars)
  
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
