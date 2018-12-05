GLM_EnKF <- function(
  nsteps,
  full_time,
  workingGLM,
  nmembers,
  nstates,
  npars,
  x,
  modeled_depths,
  surface_height,
  Qt_pars,
  wq_start,
  wq_end,
  met_file_names,
  include_wq,
  nMETmembers,
  Qt,
  spin_up_days,
  z,
  z_states,
  psi,
  alpha,
  WEATHER_UNCERTAINITY,
  PROCESS_UNCERTAINITY,
  INITIAL_CONDITION_UNCERTAINITY,
  PARAMETER_UNCERTAINITY){
  
  
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
      
      tmp <- update_temps(curr_temps = round(x[i-1,m,temp_start:temp_end],3),modeled_depths,workingGLM)
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
      if(i > (hist_days+1)){
        update_var(met_file_names[1+met_index],'meteo_fl',workingGLM)
        update_var(paste0('FCR_inflow.csv'),'inflow_fl',workingGLM)
        update_var(paste0('FCR_spillway_outflow.csv'),'outflow_fl',workingGLM)
      }else{
        update_var(met_file_names[1],'meteo_fl',workingGLM)
        update_var(paste0('FCR_inflow.csv'),'inflow_fl',workingGLM)
        update_var(paste0('FCR_spillway_outflow.csv'),'outflow_fl',workingGLM)
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
              GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = modeled_depths,vars = glm_output_vars)
              x_star[m,1:(nstates-npars)] <- c(GLM_temp_wq_out$output)
            }else{
              GLM_temp_wq_out <- get_glm_nc_var_all_wq(ncFile = 'output.nc',z_out = modeled_depths,vars = 'temp')
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
  
  return(list(x = x, x_restart = x_restart, Qt_restart = Qt_restart, x_prior = x_prior))
}