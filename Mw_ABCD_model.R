
#-----------------------------------------------------------------------------#
# Title: Hydrology calibration for the Mwache river catchment
# By:    M.Umit Taner
# Date:  September, 2015
#-----------------------------------------------------------------------------#

  require(lubridate)
  scripts_folder <- "/Users/Umit/Dropbox/Research/Scripts/"
  source(paste0(scripts_folder, "HYDROLOGY.R"))
  source(paste0(scripts_folder, "SCE_OPTIM.R"))

# Calibration data -------------------------------------------------------------

  #Dates vector for the monthly climate forcings
  clim_date  <- seq.Date(as.Date("1950-01-1"), as.Date("1999-12-1"), by="month")
  #Dates vector for observed monthly streamflow data
  qobs_date  <- seq.Date(as.Date("1976-06-1"), as.Date("1990-05-1"), by="month")
  #Dates vector for the selected calibration period
  calib_date <- seq.Date(as.Date("1978-01-1"), as.Date("1989-12-1"), by="month")
  #Calibration period
  index_clim <- which(clim_date %in% calib_date)
  index_qobs <- which(qobs_date %in% calib_date)
  
  #Catchment area (in km2)
  area <- 2250
  #Read-in monthly climate forcings (Prec, Tavg, Tmin, Tmax)
  climate <- read.csv("inputs/historical_met_forcings.csv") %>% tbl_df()
   #Read-in monthly streamflow time-series (cms) & convert to mm
  q_obs <- read.csv("inputs/Obs_monthly_flow_cms.csv", skip=1) %>%  tbl_df() %>%
    mutate(Flow_mm = Flow * 86400 * days_in_month(qobs_date) * 10^3 / (area * 10^6))

  #Hydro-met data for calibration period
  qobs_cal <- q_obs[index_qobs, "Flow_mm"] %>% unlist(use.names = FALSE)
  prec_cal <- climate[index_clim,  "Prec"] %>% unlist(use.names = FALSE)
  tavg_cal <- climate[index_clim,  "Tavg"] %>% unlist(use.names = FALSE)
  tmax_cal <- climate[index_clim,  "Tmax"] %>% unlist(use.names = FALSE)
  tmin_cal <- climate[index_clim,  "Tmin"] %>% unlist(use.names = FALSE)
  tdif_cal <- tmax_cal - tmin_cal

  #Calculate PET with the Hargreaves equation (mm/month)
  PE_cal <- HARGREAVES(calib_date, tavg_cal, tdif_cal, BasinLat = 3.5) %>%
    unlist(use.names = FALSE)


# ABCD calibration (Princeton) -------------------------------------------------

  #Parameter bounds for the abcd model 
  lowerb  <- c(0.0000 ,50   ,0.0 ,0.000)
  upperb  <- c(0.9999 ,5000 ,1.0 ,1.000)
  par.ini <- c(0.2000 ,200  ,0.7 ,0.002)
 
  #Model calibration (metric=Kling Gupta Efficiency (KGE))
  result <- SCEoptim(
    ABCD_SCE_CALIBRATE_NO_SNOW,
    par    = par.ini,
    P      = prec_cal,
    PE     = PE_cal,
    Q.obs  = qobs_cal,
    metric = "KGE",
    lower  = lowerb, 
    upper  = upperb,
    S_ini  = 1, 
    G_ini  = 1)

  #Calibration (September 2015)
  calib_par <- result$par
  #calib_par  <- c(1.202155e-03, 2.517880e+02, 9.170601e-01, 1.089352e-04) #(09.2014)
  #result$par <- c(0.54988, 293.87242, 0.9, 0.000014)

  #Simulate Model and evaluate results
  qsim_cal <- ABCD_QEST_NO_SNOW(
    calib_par, 
    P = prec_cal, 
    PE = PE_cal, 
    S_ini = 400, G_ini = 200)

  gof <- Model_GOF(
    date = calib_date[-(1:12)],
    obs  = qobs_cal[-(1:12)],
    sim  = qsim_cal[-(1:12)],
    unit = "mm per day")

 gof[[1]]$Time_Series
 gof[[1]]$Flow_Duration_Curve
 gof[[1]]$QQ_plot
 gof[[1]]$Monthly_Means
 
 gof[[2]]$NSE
 gof[[2]]$RMSE
 gof[[2]]$KGE

################################################################################  
  
  #Calibration with different data sources
  prec_dat <- read.csv("./inputs/grid_prec_long_mm.csv") %>% 
    as_data_frame() %>%
    select(-X) %>%
    spread(key = Source, value = value) %>%
    mutate(Date = hist_dates) %>%
    filter(Date %in% calib_date)

 result <- list()
 cal_pars <- array(0, dim = c(4, ncol(prec_dat)-1))
 qsim_lst <- array(0, dim = c(length(prec_cal), ncol(prec_dat)-1))
 KGE_lst  <- list()
 
 #write.csv(x = cal_pars, file = "met_cal_pars.csv") 
 
  for (i in 1:3) {
   
    #Model calibration 
    result[[i]] <- SCEoptim(
      ABCD_SCE_CALIBRATE_NO_SNOW,
      par    = par.ini,
      P      = unlist(prec_dat[, i+1], use.names = FALSE),
      PE     = PE_cal,
      Q.obs  = qobs_cal,
      metric = "KGE",
      lower  = lowerb, 
      upper  = upperb,
      S_ini  = 1, 
      G_ini  = 1) 
  
    cal_pars[,i] <- result[[i]]$par
    
    #Simulated streamflow
    qsim_lst[, i] <- ABCD_QEST_NO_SNOW(
      parm  = cal_pars[, i], 
      P     = unlist(prec_dat[, i+1], use.names = FALSE),
      PE    = PE_cal, 
      S_ini = 1, 
      G_ini = 1)

    #GOF measures
    KGE_lst[[i]] <- Model_GOF(
      date = calib_date[-(1:12)],
      obs  = qobs_cal[-(1:12)],
      sim  = qsim_lst[, i][-(1:12)],
      unit = "mm per day")
    
  }

 ################################################################################ 

 
 cal_pars <- read.csv(file = "./inputs/met_calibration_pars.csv") 
 grid_sim_flow <- data_frame(Date = calib_date, CRU = NA, GPCC = NA, Princeton = NA)
 
 for (i in 1:3) {
 
 grid_sim_flow[,i+1] <- ABCD_QEST_NO_SNOW(
   parm  = cal_pars[, i], 
   P     = unlist(prec_dat[, i+1], use.names = FALSE),
   PE    = PE_cal, 
   S_ini = 1, 
   G_ini = 1)
 }
 
 grid_sim_flow$Observed <- qobs_cal
 
 plots <- pl_TimeSeries(date = as.Date(grid_sim_flow$Date), 
                        text_size = 12,
                        data = grid_sim_flow[,-1], unit = "mm per month", 
                        agg_type = "sum", titles = FALSE, Guides = TRUE,
                        color_set = c("#d7191c", "#a6d96a", "#2c7bb6", "#000000"))

ggsave(plot = plots$fdc, file = "grid_fdc.png", height = 4, width = 5)
ggsave(plot = plots$boxplot,  file = "grid_boxplot.png",  height = 5, width = 10) 
ggsave(plot = plots$Seasonal, file = "grid_seasonal.png", height = 5, width = 10)
