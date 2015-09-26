
################################################################################  

# CONVERSION FACTORS -----------------------------------------------------------

  cmd_to_mcmm <- 365/1e6
  mcmm_to_cmd <- 1e6/365
  

################################################################################  
  
# CLIMATE PARAMETERS -----------------------------------------------------------

  #Historical climate
  hist_dates <- seq.Date(as.Date("1950-1-1"),as.Date("1999-12-1"), by="month")
  hist_forcings <- read.csv(
    file = "./inputs/GCM_forcings/Mwache_1_Obs.csv", header = FALSE)
  colnames(hist_forcings) <- c("Year", "Month", "Prec", "Tavg","Tmin","Tmax")

################################################################################  
  
# HYDROLOGY PARAMETERS ---------------------------------------------------------

  #ABCD MODEL PARAMETERS
  par_calib  <- c(0.54988, 293.87242, 0.9, 0.000014) #calibrated parameters
  
  #WATER RESOURCES MODEL PARAMETERS
  Days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  monthly_factors <- Days/sum(Days)
  
  #Net evaporation (m/month)
  catchment_area  <- 2250 #square kilometers
  netevap <- c(132, 168, 135, 44, -23, 51, 45, 74, 104, 77, 45, 88) * 1e-3
  

################################################################################  
  
# SYSTEM MODEL PARAMETERS ------------------------------------------------------

  #RESERVOIR DESIGN ALTERNATIVES (Labels, Volume (MCM), PV_COST (MUSD)
  designs <- data_frame(
    Labels = paste0(seq(40,140,20), "_MCM"),
    Volume = c(40, 60, 80, 100, 120, 140),
    PVCost = c(45.7778, 59.9176, 74.4364, 88.1147, 100.1927, 109.7938))

  #RATING CURVE (From Tahal, 2013) 
  rating_curve <- data.frame(
    volume    = c(0,1,5,15,30,50,81.5,110,160),
    elevation = c(14, 24, 34, 44, 54, 64, 75.6, 81.1, 84))
  
  #ESTIMATIONS BY APPROXFUN PACKAGE
  e_approx <- approxfun(rating_curve[,"volume"], rating_curve[,"elevation"])
  s_approx <- approxfun(rating_curve[,"elevation"], rating_curve[,"volume"])
  
  #DEMAND SECTOR
  demand_annual <- data_frame(
    Year = c(2015, 2035),
    Total_cmd = c(364000, 887000),
    Total_mcmm = Total_cmd * cmd_to_mcmm)

  irrig <- data_frame(
    months = month.abb,
    demand = c(3.95, 3.17, 2.69, 2.86, 1.93, 3.71, 4.79, 1.94, 1.93, 1.42, 0.82, 4.54))
  

