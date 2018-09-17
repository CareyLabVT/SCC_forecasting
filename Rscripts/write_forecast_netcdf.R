write_forecast_netcdf <- function(x,full_time,Qt,the_depths_init,save_file_name,x_restart,Qt_restart,time_of_forecast,hist_days){
  
  ncfname <- paste0(save_file_name,'.nc')
  #Set dimensions
  ens <- seq(1,dim(x)[2],1)
  depth <- the_depths_init
  t <- as.numeric(as.POSIXct(full_time),tz="EST5EDT",origin = '1970-01-01 00:00.00 UTC')
  states <- seq(1,dim(x)[3],1)
  
  #Set variable that states whether value is forecasted
  forecasted <- rep(1,length(t))
  forecasted[1:(hist_days+1)] <- 0
  
  #Create summary output
  mean_temp <- array(NA,dim=c(length(t),length(depth)))
  upper95_temp <- array(NA,dim=c(length(t),length(depth)))
  lower95_temp <- array(NA,dim=c(length(t),length(depth)))
  for(i in 1:length(depth)){
    for(j in 1:length(t)){
      mean_temp[j,i] <- mean(x[j,,i])
      lower95_temp[j,i] <- quantile(x[j,,i],0.025)
      upper95_temp[j,i] <- quantile(x[j,,i],0.975)
    }
  }
  
  #Define dims
  ensdim <- ncdim_def("ens",units = "",vals = ens, longname = 'ensemble member') 
  depthdim <- ncdim_def("z",units = "meters",vals = as.double(depth), longname = 'Depth from surface') 
  timedim <- ncdim_def("time",units = 'seconds', longname = 'seconds since 1970-01-01 00:00.00 UTC',vals = t)
  statedim <- ncdim_def("states",units = '', vals = states)
  
  #Define variables
  fillvalue <- 1e32
  dlname <- 'temperature'
  tmp_def <- ncvar_def("temp","deg_C",list(timedim,ensdim,depthdim),fillvalue,dlname,prec="single")
  dlname <- 'temperature_mean'
  tmp_mean_def <- ncvar_def("temp_mean","deg_C",list(timedim,depthdim),fillvalue,dlname,prec="single")
  dlname <- 'temperature_upperCI'
  tmp_upper_def <- ncvar_def("temp_upperCI","deg_C",list(timedim,depthdim),fillvalue,dlname,prec="single")
  dlname <- 'temperature_lowerCI'
  tmp_lower_def <- ncvar_def("temp_lowerCI","deg_C",list(timedim,depthdim),fillvalue,dlname,prec="single")
  dlname <- 'zone 1 temperature'
  par1_def <- ncvar_def("zone1temp","deg_C",list(timedim,ensdim),fillvalue,dlname,prec="single")
  dlname <- 'zone 2 temperature'
  par2_def <- ncvar_def("zone2temp","deg_C",list(timedim,ensdim),fillvalue,dlname,prec="single")
  dlname <- 'Kw'
  par3_def <- ncvar_def("Kw","unitless",list(timedim,ensdim),fillvalue,dlname,prec="single")
  dlname <- 'restart covariance matrix'
  Qt_restart_def <- ncvar_def("Qt_restart","-",list(statedim,statedim),fillvalue,dlname,prec="float")
  dlname <- 'matrix for restarting EnKF'
  x_def <- ncvar_def("x_restart","-",list(ensdim,statedim),fillvalue,dlname,prec="float")
  
  fillvalue <- -99
  dlname <- '0 = historical; 1 = forecasted'
  forecast_def <- ncvar_def("forecasted","-",list(timedim),fillvalue,longname = dlname,prec="integer")
  
  # create netCDF file and put arrays
  ncout <- nc_create(ncfname,list(tmp_def,forecast_def,tmp_mean_def,tmp_upper_def,tmp_lower_def,x_def,Qt_restart_def,par1_def,par2_def,par3_def),force_v4=T)
  ncvar_put(ncout,tmp_mean_def,mean_temp)

  ncvar_put(ncout,tmp_upper_def,upper95_temp)

  ncvar_put(ncout,tmp_lower_def,lower95_temp)

  ncvar_put(ncout,tmp_def,x[,,1:29])

  ncvar_put(ncout,par1_def,array(x[,,30]))
  
  ncvar_put(ncout,par2_def,array(x[,,31]))
  
  ncvar_put(ncout,par3_def,array(x[,,32]))

  ncvar_put(ncout,forecast_def,as.array(forecasted))

  ncvar_put(ncout,x_def,as.matrix(x_restart))

  ncvar_put(ncout,Qt_restart_def,as.matrix(Qt_restart))

  
  #Global file metadata
  ncatt_put(ncout,0,"title",'Falling Creek Reservoir forecast')
  ncatt_put(ncout,0,"institution",'Virginia Tech')
  ncatt_put(ncout,0,"source",'GLMV3')
  #ncatt_put(ncout,0,"references",references$value)
  history <- paste('Run date:',time_of_forecast, sep=", ")
  ncatt_put(ncout,0,"history",history)
  #ncatt_put(ncout,0,"Conventions",)
  
  nc_close(ncout)
}

