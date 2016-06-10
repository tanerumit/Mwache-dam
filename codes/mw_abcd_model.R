
#-----------------------------------------------------------------------------#
# Title: Hydrology model calibration for the Mwache river catchment
#-----------------------------------------------------------------------------#

  source("Mw_Input.R")
  scripts_dir <- "/Users/Umit/Dropbox/Research/Scripts/"
  source(paste0(scripts_dir, "PLOTTING.R"))
  source(paste0(scripts_dir, "WATER_RESOURCES.R"))
  source(paste0(scripts_dir, "SCE_OPTIM.R"))

# ------------------------------------------------------------------------------

# CLIMATE / GAGE DATA ----------------------------------------------------------

  #Precipitaton datasets available
  datasets <- c("CRU", "GPCC", "Princeton")

  #Read-in climate forcings, streamflow, and grid precip forcings
  clim_forcing_data <- read_csv("./inputs/forcing_gcm.csv", skip = 1) %>%
    filter(Index == 1)
  gage_data <- read_csv("./Inputs/flow_gage_cms.csv", skip=1)
  prec_grid_data <- read_csv("./inputs/precip_grid_mm.csv")

  #Climate forcings
  clim <- clim_forcing_data %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep="-")),
           Tdif = Tmax - Tmin,
           PE = HARGREAVES(Date, Tavg, Tdif, BasinLat)) %>%
    select(Date, Prec:Tmin, Tdif, PE)

 #write_csv(x = clim, path = "./inputs/forcing_princeton.csv")

  #Read-in monthly streamflow time-series (cms) & convert to mm
  gage <- gage_data %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep = "-")),
      Q_obs_mm = Flow_obs * (86.4 * days_in_month(month(Date))) /area,
      Q_int_mm = Flow_int * (86.4 * days_in_month(month(Date))) /area) %>%
    select(Date, Q_obs_mm, Q_int_mm)

# ------------------------------------------------------------------------------

# MODEL CALIBRATION  -----------------------------------------------------------

  #Set abcd model parameter bounds
  lowerb  <- c(0.0,      15, 0.0, 0.0)
  upperb  <- c(1.0,     500, 1.0, 0.00001)
  par.ini <- c(0.0010,  211, 0.9, 0.00005)
  #par  <- c(1.202155e-03, 2.517880e+02, 9.170601e-01, 1.089352e-04) #(09.2014)
  #par <- c(0.54988, 293.87242, 0.9, 0.000014)

  #Set the calibration period
  date_cal_beg <- as.Date("1978-01-1")
  date_cal_end <- as.Date("1989-12-1")
  date_cal <- seq.Date(date_cal_beg, date_cal_end, by = "month")

  #Calibration parameters
  calib_pars <- matrix(NA, nrow = 4, ncol =length(datasets))
  colnames(calib_pars) <- datasets

  #Loop through each grid-dataset
  for (i in 1:length(datasets)) {

    #Select associated data
    prec_i <- filter(prec_grid_data, Source == datasets[i]) %>%
      select(value) %>% rowSums()

    #Append prec data to the climate data table
    clim <- mutate(clim, Prec = prec_i)

    #Adjustments for the calibration period
    clim_cal <- filter(clim, Date %in% date_cal)
    gage_cal <- filter(gage, Date %in% date_cal)

    #Simulate streamflow (spin-up period only!)
    spin_up <- 24
    out <- ABCD_QEST(par.ini, print.all = TRUE, S_ini = 100, G_ini = 2,
      P = clim_cal$Prec[1:spin_up], PE = clim_cal$PE[1:spin_up])

    #Model calibration (metric=Kling Gupta Efficiency (KGE))
    result <- SCEoptim(ABCD_CALIBRATE, par = par.ini, na.rm = FALSE,
      P = clim_cal$Prec[-c(1:spin_up)], PE = clim_cal$PE[-c(1:spin_up)],
      Q.obs = gage_cal$Q_int_mm[-c(1:spin_up)],
      metric = "KGE",lower = lowerb, upper = upperb,
      S_ini  = out$S[spin_up], G_ini  = out$G[spin_up])

    #Calibration parameters
    calib_pars[, datasets[i]] <- result$par

  }

  #***********************************************************************
  #write_csv(x = as.data.frame(calib_pars), "./inputs/calib_pars_11_2015.csv")
  #***********************************************************************


# ------------------------------------------------------------------------------

# CALIBRATION /MET-DATASET PERFORMANCE -----------------------------------------

  #Load calibration parameters
  q_estimate <- list()

  #Choose dataset
  for (i in 1:length(datasets)) {

    #Select associated data
    prec_i <- filter(prec_grid_data, Source == datasets[i]) %>%
      select(value) %>% rowSums()

    #Append prec data to the climate data table
    clim_sim <- mutate(clim, Prec = prec_i) %>% filter(!is.na(Prec))

    #Parameters
    pars <- hydro_pars[[datasets[i]]]

    #Adjustments for the calibration period
    q_estimate[[datasets[i]]] <- clim_sim %>%
      mutate(Q_sim = ABCD_QEST(pars, S_ini  = 100, G_ini  = 2,
        P = clim_sim$Prec, PE = clim_sim$PE))
  }

  ### 1. CHECK CALIBRATION PERFORMANCES
  #calibration period
  date_cal_beg <- as.Date("1980-01-1")
  date_cal_end <- as.Date("1989-12-1")
  date_cal <- seq.Date(date_cal_beg, date_cal_end, by = "month")

  dataset_i <- "CRU"
  Q_obs <- gage %>% filter(Date %in% date_cal) %>% select(Q_int_mm) %>% rowSums()
  Q_sim <- q_estimate[[dataset_i]] %>% filter(Date %in% date_cal) %>%
    select(Q_sim) %>% rowSums()

  GOF <- Model_GOF(Date = date_cal, Observed = Q_obs, Simulated = Q_sim)
  GOF[[1]][[3]]
  GOF[[2]]

  #Long-term simulation to check "trends"
  pars <- hydro_pars[[dataset_i]]
  q_est <- ABCD_QEST(pars, print.all = F, S_ini = 100, G_ini = 2,
    P = clim_sim$Prec, PE = clim_sim$PE)
  df <- mutate(clim_sim, q_est = q_est)
  p  <- ggplot(df, aes(x = Date, y = q_est)) + geom_line(); p

  grid_qest <- bind_rows(q_estimate, .id = "Dataset") %>%
    filter(Date %in% date_cal) %>%
    select(Dataset, Date, Q_sim) %>%
    spread(key = Dataset, value = Q_sim)

  df <- grid_qest %>% mutate(Observed = Q_obs) %>%
    gather(key = Data, value = value, - Date) %>%
    mutate(Data = factor(Data, levels = c("Observed", "Princeton", "GPCC", "CRU")))

  #Time-series of flow
  p <- ggplot(df, aes(x = Date, y = value, color = Data))  +
    theme_bw(base_size = 11) +
    geom_line(size = 0.6) +
    scale_x_date(breaks = as.Date("1980-01-01") + years(seq(0,10,2)), date_labels = "%Y") +
    scale_color_manual(values = c(
      "Observed" = "black", "Princeton" = "#d7191c",
      "GPCC" = "#a6d96a", "CRU" = "#2c7bb6")) +
    labs(x = "", y = "flow (mm)")

  ggsave("flow_timeseries.png", p, width = 7, height = 4)

  #Monthly mean flow
  df2 <- df %>%
    mutate(month = month(Date)) %>%
    group_by(month, Data) %>%
    summarize(value = mean(value))

  p <- ggplot(df2, aes(x = as.factor(month), group = Data,
                       y = value, color = Data))  +
    theme_bw(base_size = 11) +
    geom_line(size = 0.6) +
    scale_x_discrete(labels = month.abb) +
    scale_color_manual(values = c(
      "Observed" = "black", "Princeton" = "#d7191c",
      "GPCC" = "#a6d96a", "CRU" = "#2c7bb6")) +
    labs(x = "", y = "flow (mm)")

  ggsave("flow_means.png", p, width = 7, height = 4)


  #ggsave(plot = pl$fdc, file = "grid_fdc.png", height = 4, width = 5)
  #ggsave(plot = pl$boxplot, file = "grid_boxplot.png",  height = 5, width = 10)
  #ggsave(plot = pl$Seasonal, file = "grid_seasonal.png", height = 5, width = 10)


# ------------------------------------------------------------------------------

# STAT ANALYSIS ----------------------------------------------------------------

 clim_monthly_means <- clim %>%
    mutate(Month = month(Date)) %>%
    select(-Date) %>%
    group_by(Month) %>%
    summarize_each(funs(mean))

