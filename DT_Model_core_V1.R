
#------------------------------------------------------------------------------#
#   Title: Decision Tree Project - Mwache Systems Model
#   By   : Mehmet Umit Taner 
#   Date : September 25, 2015                                                     
#------------------------------------------------------------------------------# 

# TASK - LIST 
# 1. Simulate streamflow from rainfall (one-node) 
# 2. Route through Mwache reservoir 
# 3. Demand nodes (domestic + irrigation)

  source("/Users/umit/GDrive/Research/Scripts/HYDROLOGY.R")
  source("/Users/umit/GDrive/Research/Scripts/SAFE_YIELD.R")
  source("Mw_input_pars.R")


# Hydrology simulation  --------------------------------------------------------

  #Set forcings 
  dates <- hist_dates
  forcings <- hist_forcings

  #Calculate PET using HARGREAVES
  Tdif <- forcings$Tmax - forcings$Tmin
  PE <- HARGREAVES(dates, tavg = forcings$Tavg, tdif = Tdif, BasinLat = 3.5)

  #Simulate streamflow, in mm_month
  Q_sim_mm  <- ABCD_QEST_NO_SNOW(parm = par_calib, P = forcings$Prec, PE = PE) 
  Q_sim_mcm <- Q_sim_mm * (catchment_area * 1e-3)
  
# Water Resources system simulation --------------------------------------------

  #1) Yield-reliablity analysis ************************************************
  
  #Simulate reservoir performance over the given period 
  yield <- c(seq(50, 75, 5), seq(80, 140, 10))
  
  rel <- vector("numeric", length = length(yield_range))
  out <- list(length(yield_range))
  
  for (k in 1:length(yield_range)) {
    
    out <- RESERVOIR_SIM(
      beg_year   = 1950,     
      beg_month  = 1,    
      capacity   = 120,    
      s_ini      = 0.8,        
      top_cons   = 120,      
      top_buff   = 100,      
      top_dead   = 20,     
      e_v_curve  = e_approx,      
      v_e_curve  = s_approx,  
      demand_yr  = yield[k], 
      irrig_mon  = irrig$demand,
      Inflow     = Q_sim_mcm,       
      nevap      = netevap,      
      cycle      = TRUE,
      env_flow   = 2)

    rel[k] <- with(out, length(which((Release - Demand) == 0))/ nrow(out))
  
  }
  

  
  