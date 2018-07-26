#FUNCTIONS Shared between the MCMC and EnKF

#Set GLM Initial conditions for temperature at each layer from first observations
update_temps <- function(curr_temps,curr_depths,workingGLM){
  origNML = read_nml(file.path(workingGLM,'glm3.nml'))
  index1 = NA; index2 = NA;index3 = NA; index4 = NA
  for (g in 1:length(origNML)) {
    for (q in 1:length(origNML[[g]])) {
      if (names(origNML[[g]][q]) == "the_temps") {
        temps = as.numeric(as.character(unlist(origNML[[g]][q])))
        index1 = g; index2 = q; 
      }
      if (names(origNML[[g]][q]) == "the_depths") {
        depths = as.numeric(as.character(unlist(origNML[[g]][q])))
        index3 = g; index4 = q
      }
    }
  }
  temp_inter = approxfun(curr_depths,curr_temps,rule=2)
  init_temps = temp_inter(depths)
  char_temps = paste(init_temps, collapse = ', ')
  holder2 = unlist(origNML[[index1]][index2])
  holder2[1:length(holder2)] = init_temps
  holder2 = list(holder2)
  holder3 = unlist(origNML[[index3]][index4])
  holder3[1:length(holder3)] = depths
  holder3 = list(holder3)
  origNML[[index1]][index2] = holder2
  origNML[[index3]][index4] = holder3
  write_nml(origNML, file.path(workingGLM,'glm3.nml'))
  return(list(depths,init_temps))
}


update_var <- function(var_value,var_name,workingGLM){
  origNML = read_nml(file.path(workingGLM,'glm3.nml'))
  index1 = NA; index2 = NA
  for (g in 1:length(origNML)) {
    for (q in 1:length(origNML[[g]])) {
      if (names(origNML[[g]][q]) == var_name) {
        index1 = g; index2 = q; 
      }
    }
  }
  holder2 = unlist(origNML[[index1]][index2])
  holder2[1:length(var_value)] = var_value
  holder2 = list(holder2[1:length(var_value)])
  origNML[[index1]][index2] = holder2
  write_nml(origNML, file.path(workingGLM,'glm3.nml'))
}

update_time <- function(start_value,stop_value,workingGLM){
  origNML = read_nml(file.path(workingGLM,'glm3.nml'))
  index1 = NA; index2 = NA; index3 = NA; index4 = NA
  for (g in 1:length(origNML)) {
    for (q in 1:length(origNML[[g]])) {
      if (names(origNML[[g]][q]) == 'start') {
        index1 = g; index2 = q; 
      }
      if (names(origNML[[g]][q]) == 'stop') {
        index3 = g; index4 = q; 
      }
    }
  }
  origNML[[index1]][index2] = start_value
  origNML[[index3]][index4] = stop_value
  write_nml(origNML, file.path(workingGLM,'glm3.nml'))
}

get_glm_nc_var <- function(ncFile,z_out,var = 'temp'){
  glm_nc <- nc_open(ncFile)
  tallest_layer <- ncvar_get(glm_nc, "NS")
  elev <- ncvar_get(glm_nc, "z")
  temp <- ncvar_get(glm_nc, var)
  nc_close(glm_nc)
  
  elev_surf = get_surface_height(ncFile)
  
  max_i <- tallest_layer[length(tallest_layer)]
  
  elev <- elev[1:max_i, length(tallest_layer)]
  temp <- temp[1:max_i, length(tallest_layer)]
  
  num_step <- length(tallest_layer)
  num_dep <- length(z_out)
  temp_out <- rep(NA,num_dep)
  tme = num_step
  elevs_out <- elev_surf[tme, 2] - z_out
  
  elevs = elev
  temps = temp
  
  num_z <- max_i
  layer_mids <- c(elevs[1]/2, elevs[1:num_z-1] + diff(elevs)/2)
  temps_re <- c(temps[1], temps, tail(temps,1))
  elevs_re <- c(0, layer_mids, tail(elevs, 1))
  temps <- approx(x = elevs_re, y = temps_re, xout = elevs_out)$y
  return(temps)
}

update_phyto <- function(p_initial,nml_name = 'aed2_phyto_pars.nml'){
  
  if(length(p_initial)<6){
    print('number of phyto group does not equal 6')
  }   
  nml_file <- file.path(workingGLM,nml_name)
  c <- file(nml_file, "r")
  fileLines <- readLines(c)
  close(c)
  lineStart <- substr(fileLines, 1, 1)
  ignoreLn <- lineStart == "!" | fileLines == ""
  lineStart <- lineStart[!ignoreLn]
  fileLines <- fileLines[!ignoreLn]
  l <- strsplit(fileLines, ",")
  
  l[[2]][2] = p_initial[1]
  l[[3]][2] = p_initial[2]
  l[[4]][2] = p_initial[3]
  l[[5]][2] = p_initial[4]
  l[[5]][2] = p_initial[5]
  
  zz <- file(nml_file, open = "wt")
  sink(zz)
  cat(noquote(paste(paste(l[[1]],collapse = ''),'\n')))
  cat(noquote(paste(paste(l[[2]],collapse = ','),'\n')))
  cat(noquote(paste(paste(l[[3]],collapse = ','),'\n')))
  cat(noquote(paste(paste(l[[4]],collapse = ','),'\n')))
  cat(noquote(paste(paste(l[[5]],collapse = ','),'\n')))
  cat(noquote(paste(paste(l[[6]],collapse = ','),'\n')))
  cat(noquote(paste(paste(l[[7]],collapse = ','),'\n')))
  sink()
}
