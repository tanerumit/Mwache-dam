
#-----------------------------------------------------------------------------#
# Title: Hydrology calibration for the Mwache river catchment
# By:    M.Umit Taner
# Date:  September, 2015
#-----------------------------------------------------------------------------#

  require(lubridate)
  source("/Users/umit/GDrive/Research/Scripts/HYDROLOGY.R")
  source("/Users/umit/GDrive/Research/Scripts/SCE_OPTIM.R")

# ABCD Model Calibration -------------------------------------------------------

  #Simulation dates
  clim_date  <- seq.Date(as.Date("1950-01-1"), as.Date("1999-12-1"), by="month")
  qobs_date  <- seq.Date(as.Date("1976-06-1"), as.Date("1990-05-1"), by="month")
  calib_date <- seq.Date(as.Date("1978-01-1"), as.Date("1989-12-1"), by="month")

  #Hydroclimate data
  area <- 2250
  climate <- read.csv("inputs/forcings_historic.csv") %>% tbl_df()
  q_obs   <- read.csv("inputs/observed_flow_cms.csv", skip=1) %>% tbl_df() %>%
    mutate(Flow_mm = Flow_obs * 86400 * days_in_month(qobs_date) * 10^3 / (area * 10^6))

  #Calibration period
  index_clim <- which(clim_date %in% calib_date)
  index_qobs <- which(qobs_date %in% calib_date)

  #Met data for calibratoin
  qobs_cal <- q_obs[index_qobs, "Flow_mm"] %>% unlist(use.names = FALSE)
  prec_cal <- climate[index_clim, "Prec"] %>% unlist(use.names = FALSE)
  tavg_cal <- climate[index_clim, "Tavg"] %>% unlist(use.names = FALSE)
  tmax_cal <- climate[index_clim, "Tmax"] %>% unlist(use.names = FALSE)
  tmin_cal <- climate[index_clim, "Tmin"] %>% unlist(use.names = FALSE)
  tdif_cal <- tmax_cal - tmin_cal

  #ET (mm/month)
  PE <- HARGREAVES(Dates = calib_date, tavg = tavg_cal, tdif = tdif_cal, BasinLat = 3.5) %>%
    unlist(use.names = FALSE)

  #abcd model parameter bounds
  lowerb    <- c(0.001  ,50   ,0.0 ,0.001)
  upperb    <- c(0.9999 ,5000 ,1.0 ,1)
  par.start <- c(0.2    ,200  ,0.7 ,0.002)

  #Model calibration (metric=Kling Gupta Efficiency (KGE))
  result <- SCEoptim(ABCD_SCE_CALIBRATE_NO_SNOW,
    par = par.start,
    P = prec_cal,
    PE = PE,
    Q.obs = qobs_cal,
    metric = 'KGE',
    lower = lowerb, upper = upperb,
    S_ini = 1, G_ini = 1)

  #Calibration (September 2015)
  calib_par <- result$par
  #calib_par  <- c(1.202155e-03, 2.517880e+02, 9.170601e-01, 1.089352e-04) #(09.2014)
  #result$par <- c(0.54988, 293.87242, 0.9, 0.000014)

  #Simulate Model and evaluate results
  qsim_cal <- ABCD_QEST_NO_SNOW(calib_par, P = prec_cal, PE = PE, S_ini =1, G_ini =1)

  round(calib_par, 3)
# Simulate calibration ---------------------------------------------------------

  plots <- HydroFitPlots(
    date = calib_date[-(1:12)],
    qobs = qobs_cal[-(1:12)],
    qsim = qsim_cal[-(1:12)],
    units = "cm day")

  stats <- HydroFitStats(
    qobs = qobs_cal[-(1:12)],
    qsim = qsim_cal[-(1:12)])


