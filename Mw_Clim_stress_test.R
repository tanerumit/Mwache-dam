

#------------------------------------------------------------------------------#
# Title: Climate stress test of the alternative Mwache system designs
# By:    M.Umit Taner
# Date:  August, 2015
#------------------------------------------------------------------------------#

# Hydrologic simulations are based on the abcd model
# Water system operations are based on WEAP

  source("/Users/umit/GDrive/Research/Scripts/HYDROLOGY.R")
  source("/Users/umit/GDrive/Research/Scripts/YIELD3.R")
  
# Climate Stress Test Design ---------------------------------------------------

  # This section defines the parameters used for the climate stress test
  # Parameters: Annual DeltaT, Annual DeltaP, and climate variability

  # Historical climate means
  # Which years??? What data???
  hist_clim <- list(Temp=25.21428,  Prec=845.4606)

  #Delta factors for climatic changes (mid-values)
  DeltaT_mid <- seq(0,2.5,0.5)
  DeltaP_mid <- seq(70,140,10)
  T_mid <- hist_clim$Temp + DeltaT_mid
  P_mid <- hist_clim$Prec * DeltaP_mid/100

  #Climate change matrix
  cc_matrix <- expand.grid(Tavg = DeltaT_mid, Prec = DeltaP_mid)
  cc_index  <- expand.grid(Tavg=1:length(DeltaT_mid), Prec=1:length(DeltaP_mid))

  #Climate matrix (climate change + climate variability)
  nvar_num <- 10
  clim_matrix <- expand.grid(Tavg = DeltaT_mid, Prec = DeltaP_mid, NVar=1:nvar_num) %>%
    mutate(Clim = 1: nrow(.), CC = rep(1:nrow(cc_index),nvar_num))
  clim_index <- expand.grid(Tavg = 1:length(DeltaT_mid), Prec = 1:length(DeltaP_mid), NVar=1:nvar_num) %>%
    mutate(Clim = 1: nrow(.), CC = rep(1:nrow(cc_index),nvar_num))
  iv_index <- sapply(unique(clim_matrix$NVar), function(x) which(clim_matrix$NVar==x))

  #Read-in previously generated climate forcings
  st_forcings <- lapply(1:nrow(clim_index), function(x)
    read.csv(paste0("Inputs/ST_forcings_022015/clim_",x,".csv"), header=T, skip=2))
  clim_num <- length(st_forcings)

  #DATA structure
  loopArray <- expand.grid(c=1:clim_num, d=1:nrow(designs))

################################################################################

# SYSTEM SIMULATION ------------------------------------------------------------

  #SIMULATE HYDROLOGY
  results_df <- array(NA, c(nrow(loopArray), 5))
  loopArray <- loopArray[1:50,]
 
  ptm <- proc.time() 
  for (s in seq_len(nrow(loopArray))) {
  
    #Set climate (c) and design (d) counters
    c <- loopArray[s,1]
    d <- loopArray[s,2]
    
    #Calculate PET using HARGREAVES
    tdif <- st_forcings[[c]]$Tmax- st_forcings[[c]]$Tmin
    PE <- HARGREAVES(Dates = dates_sim,
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
      dates         = dates_sim,   
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
  
  }
  proc.time() - ptm
  
################################################################################  
  
  #Set target demand (million m3 per month):
  TDemand= 81 * monthly_factors

  #Prepare input/output data tables
  result <- as.list(rep(NA, nrow(loopArray)))

  #Progress bar & timing
  pb <- txtProgressBar(min = 1, max = nrow(loopArray), style = 3)

  #SIMULATE HYDROLOGY & RESERVOIR SYSTEM OPERATIONS
  ptm <- proc.time()
  for (s in seq_len(nrow(loopArray))) {

    #Set progress bar
    setTxtProgressBar(pb, s)

    #Set climate (c) and design (d) counters
    c <-  loopArray[s,1]
    d <-  loopArray[s,2]

    #Calculate PET using HARGREAVES
    tdif <- st_forcings[[c]]$Tmax- st_forcings[[c]]$Tmin
    PE <- HARGREAVES(Dates = dates_sim,
      tavg=st_forcings[[c]]$Tavg, tdif = tdif, BasinLat = 3.5)

    #Simulate Hydrology & water Res operations
    result[[s]] <- WEAPSIM(
      sim_years       = year(dates_sim),
      sim_months      = month(dates_sim),
      catchment_area  = catchment_area,
      hydro_par       = par_calib,
      prec            = st_forcings[[c]]$Prec,
      PET             = PE,
      storage_cap     = designs$Volume[d],
      storage_init    = designs$Volume[d],
      top_conserv     = designs$Volume[d],
      top_buffer      = 35,
      top_inactive    = 20,
      buffer_coef     = 1.0,
      elev_vol_curve  = e_approx,
      vol_elev_curve  = s_approx,
      demand_monthly  = TDemand)

    close(pb)

  }

  #CALCULATE RESERVOIR PERFORMANCE CRITERIA
  metrics <- sapply(seq_len(length(result)),
    function(x) RES_PERFORMANCE(result[[x]])$Value)


# ------------------------------------------------------------------------------

