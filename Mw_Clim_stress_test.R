

#------------------------------------------------------------------------------#
# Title: Climate stress test of the alternative Mwache system designs
# By:    M.Umit Taner
# Date:  August, 2015
#------------------------------------------------------------------------------#
    
  # Hydrologic simulations are based on the abcd model
  # Water system operations are based on WEAP
  
  source("Mw_Input_pars.R")

  #Read-in previously generated climate forcings
  st_forcings <- lapply(1:nrow(clim_index), function(x)
    read.csv(paste0("Inputs/ST_forcings_022015/clim_",x,".csv"), header=T, skip=2))
    clim_num <- length(st_forcings)
     
  #DATA structure
   loopArray <- expand.grid(c=1:clim_num, d=1:nrow(designs))   
     
     
  #SIMULATE HYDROLOGY
  results_df <- array(NA, c(nrow(loopArray), 5))

  sim_end <- nrow(loopArray)

    
  #Progress bar & timing
  pb <- txtProgressBar(min = 1, max = sim_end, style = 3)
    
  ptm <- proc.time() 
  for (s in 1:sim_end) {
      
    setTxtProgressBar(pb, s)
    
    #Set climate (c) and design (d) counters
    c <- loopArray[s,1]
    d <- loopArray[s,2]
    
    #Calculate PET using HARGREAVES
    tdif <- st_forcings[[c]]$Tmax- st_forcings[[c]]$Tmin
    PE <- HARGREAVES(Dates = sim_dates,
      tavg=st_forcings[[c]]$Tavg, tdif = tdif, BasinLat = 3.5)
    
    #Simulate streamflow (mm)
    Q_sim_mm <- ABCD_QEST_NO_SNOW(
      par_calib, 
      P = st_forcings[[c]]$Prec, 
      PE = PE, 
      S_ini = 1, 
      G_ini = 1) 
    
    #Simulate streamflow (MCM/mo)
    Q_sim <- Q_sim_mm * (catchment_area * 1e-3)

    #Find yield for the given reliability and inflow sequence
    result <- YIELD_ESTIMATE(
      dates         = sim_dates,   
      Q             = Q_sim,                     
      capacity      = designs$Volume[d],   
      top_conserv   = designs$Volume[d],   
      top_buffer    = 35,   
      top_inactive  = 5,   
      buffer_coef   = 1,   
      hVol_curve    = e_approx,  
      Volh_curve    = s_approx,   
      evap_m        = netevap,   
      S_initial     = 1,      
      reliability   = 0.95,   
      max_iter      = 100,    
      double_cycle  = TRUE)  
    
    #Save results
    results_df[s, 1] <- c
    results_df[s, 2] <- d
    results_df[s, 3:5] <- result$Val
  
    close(pb)
    
  }
  proc.time() - ptm
  
  head(results_df)
  
  results_out <- data.frame(results_df) %>%
    as_data_frame() %>%
    select(c = X1, d = X2, iter = X3, rel = X4, yield = X5)
  
  write.csv(x=results_out, file = "safeyield_Sep25.csv")
  
  
  


