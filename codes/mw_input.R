
################################################################################
# MWACHE ANALYSIS - INPUT PARAMETERS
################################################################################

#-------------------------------------------------------------------------------
# REQUIRED DATA ----------------------------------------------------------------

  options(stringsAsFactors = FALSE)

  library(sirad, quietly = T, warn.conflicts = F)
  library(lubridate, quietly = T, warn.conflicts = F)
  library(LaplacesDemon, quietly = T, warn.conflicts = F)
  library(truncnorm, quietly = T, warn.conflicts = F)

  script_dir <- "/Users/umit/Dropbox/Research/Scripts/"
  source(paste0(script_dir, "PRIMARY.R"))
  source(paste0(script_dir, "WATER_RESOURCES.R"))
  source(paste0(script_dir, "PLOTTING.R"))
  source(paste0(script_dir, "reservoir_sim.R"))

  # climate datasets (x3)
  dataset_lst <- c("Princeton","GPCC","CRU")

  #Climate forcings (all)
  forcing <- read_csv("./data/clim_trajectories.csv") %>%
    select(dataset = Dataset, nclim = nClim, Year:Tmax)

  #Crop efficiency and Kc values for each crop type
  crop_data <- read_csv("./data/crop_data_irrig.csv", skip = 1)

  # Monthly parameters for Penmann - Montheith calculation (fixed)
  crop_pars_et0 <- read_csv("./data/crop_data_et0.csv", skip = 1)

  #gcm_meta <- read_csv("Inputs/gcm_metadata.csv", col_names = TRUE)

#-------------------------------------------------------------------------------
# STRESS TEST DESIGN -----------------------------------------------------------

  # This section defines the parameters used for the climate stress test
  # Parameters: Annual DeltaT, Annual DeltaP, and climate variability

  #Delta factors for climatic changes (full period changes)
  nvar_num <- 10
  #DelP <- seq(0.4, 1.8, 0.1)
  #DelT <- seq(0,4,1)
  DelP <- seq(0.7, 1.5, 0.1)
  DelT <- seq(0,5,1)

  # Look-up table for climate changes
  cc_tbl  <- expand.grid(
      ntemp = 1:length(DelT), nprec = 1:length(DelP)) %>%
    mutate(temp = DelT[ntemp], prec = DelP[nprec], ncc = 1:n()) %>%
    select(ncc, ntemp, nprec, temp, prec)

  # Look-up table for climate (variability + change)
  clim_tbl <- cc_tbl %>%
    expand_grid_df(., data.frame(nvar = 1:nvar_num)) %>%
    mutate(nclim = 1:n()) %>%
    select(nclim, ncc, ntemp, nprec, nvar, temp, prec)


  # ######## GENERATE CLIMATE TRAJECTORIES BASED ON CLIMATE GRID
  # future_period <- as.Date("2020-01-1") + c(0:599) * months(1)
  # future_years <- unique(year(future_period))
  # data <- read_csv("./inputs/WARM_realizations.csv")
  #
  # outmat <- expand_grid(
  #   dataset = dataset_lst, nclim = 1:nrow(clim_tbl)) %>%
  #   left_join(clim_tbl, by ="nclim") %>%
  #   mutate(dataset = as.character(dataset)) %>%
  #   arrange(dataset, nclim)
  #
  # traject <- vector(mode = 'list', length = nrow(outmat))
  #
  # #Loop through each row of the outmat and modify the climate data
  # end <- nrow(outmat)
  # for (i in 1:end) {
  #
  #   #Extract chance factors & variability index in each row
  #   dataset <- outmat$dataset[i]
  #   nclim <- outmat$nclim[i]
  #   temp <- outmat$temp[i]
  #   prec <- outmat$prec[i]
  #   nvar <- outmat$nVar[i]
  #
  #   #Vector of transient climate change factors
  #   t_trend <- rep(seq(0, temp, length.out=length(future_years)),each=12)
  #   p_trend <- rep(seq(1, prec, length.out=length(future_years)),each=12)
  #
  #   #Monthly climate table
  #   cli <- data %>% filter((Dataset == dataset) & (nVar == nvar)) %>%
  #     select(-nVar)
  #
  #   traject[[i]] <- cli %>%
  #     mutate(nclim = nclim) %>%
  #     mutate(prec = prec * p_trend, temp = temp + t_trend,
  #            Tmin = Tmin + t_trend, Tmax = Tmax + t_trend)
  # }
  #
  # #Bind all trajectories
  # traject_all <- bind_rows(traject)
  # write_csv(x = traject_all, path = "./inputs/WARM_trajectories.csv")

#------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------
# ABCD MODEL PARAMETERS --------------------------------------------------------

  BasinLat <- -3.6

  #Misc. conversions
  cmd_to_mcmm <- 365/1e6
  mcmm_to_cmd <- 1e6/365
  month_coef <- as.numeric(lubridate::days_in_month(1:12)/365)

  #Calibration parameters
  par_calib  <- c(0.54988, 293.87242, 0.9, 0.000014)
  area <- 2250 #in square kilometers
  Latitude = -4.20

  hydro_pars <- list(
    Princeton = c(0.000445674, 315.8826494,	0.905240881, 1.00E-06),
    GPCC      = c(1.19E-05,	   268.4657258,	0.897705472, 1.00E-06),
    CRU	      = c(7.01E-05,	   283.515018,	0.931386386, 1.00E-06))


#-------------------------------------------------------------------------------
# RESERVOIR MODEL PARAMETERS ---------------------------------------------------

  #Net evaporation (m/month) (NET ET = ET - precipitation)
  #netevap <- c(132, 168, 135, 44, -23, 51, 45, 74, 104, 77, 45, 88) * 1e-3

  #RESERVOIR DESIGN ALTERNATIVES (Labels, Volume (MCM), PV_COST (MUSD)
  designs <- data_frame(
    Labels = paste0(c(80,100,120,130,140), "_MCM"),
    Volume = c(80,100,120,130,140),
    PVCost = c(74.4364, 88.1147, 100.1927, 105.3601, 109.7938))

  #designs <- data_frame(
  #  Labels = paste0(seq(40,140,20), "_MCM"),
  #  Volume = c(40, 60, 80, 100, 120, 140),
  #  PVCost = c(45.7778, 59.9176, 74.4364, 88.1147, 100.1927, 109.7938))

  # Area Elevation Capacity Curve at Mwache Dam site (Final Design report, 2014)
  # elevation (m), volume (MCM), area (km^2)
  ratingc <- data_frame(
    elev = c(14,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,85.7,90,95),
    area = c(0,0.01,0.1,0.28,0.42,0.65,0.9,1.16,1.55,1.91,2.32,2.74,
      3.24,3.8,4.76,6.03,6.33,8.18,12.04),
    volume = c(0,0,0.24,1.14,2.85,5.49,9.35,14.47,21.21,29.84,40.39,
      53.02,67.95,85.52,106.86,133.77,138.1,169.29,219.8))

  #elevation = f(volume)
  f_elev <- approxfun(x = ratingc$volume, y = ratingc$elev)
  #volume = f(elevation)
  f_vol  <- approxfun(x = ratingc$elev, y = ratingc$volume)
  #area = f(volume)
  f_sarea  <- approxfun(x = ratingc$volume, y = ratingc$area)
  #p <- ggplot(rating_curve, aes(x=volume, y = elev)) + geom_line(); p

  #Critical reliability calculation
  critical_level <- 0.5

#-------------------------------------------------------------------------------
# DEMAND ESTIMATES--------------------------------------------------------------

  #80.3 MCM = 220,000 m3/d
  #67.89 MCM = 186,000 m3/d

  year_mf <- c(31,28,31,30,31,30,31,31,30,31,30,31) / 365

  #Monthly irrigation demand estimates by CES, 2014
  demand_ces_mcm <- data_frame(Month = as.factor(1:12),
    irr = c(3.95,3.17,2.69,2.86,1.93,3.71,4.79,1.94,1.93,1.42,0.82,4.54),
    env = c(0.83, 0, 0, 1.69, 1.49, 0.82, 0.76, 0.55, 0, 0, 0.4, 1.26),
    dom = 68.3 * year_mf)

  mf <- demand_ces_mcm %>%
    mutate(irr = irr/sum(irr), eco = env / sum(env), dom = year_mf)

  df <- demand_ces_mcm %>%
    select(Month, Irrigation = irr, Domestic = dom, Environmental = env) %>%
    gather(key = Uses, value = value, -Month)

  p <- ggplot(df, aes(x= Month, y = value, fill = Uses)) +
    theme_bw() +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(labels = month.abb) +
    scale_y_continuous(limits = c(0,6), breaks = seq(0,6,2)) +
    labs(x="", y="MCM")

  ggsave(plot = p, filename = "demand_ces.tiff", height = 3.5, width = 7)

  # DOMESTIC DEMAND SCENARIOS +++++++++++++++++

  # demand levels for stress test
  #demand_levels %<>% filter(Year %in% 2020:2069) %>%
  #  expand_grid_df(data.frame(Month = 1:12)) %>% as.tbl() %>%
  #  select(level, Year, Month, value_yr = value) %>%
  #  arrange(level, Year, Month) %>%
  #  mutate(value = value_yr * month_coef[Month]) %>%
  #  select(level, Year, Month, value)


  #IRRIGATION CALCULATIONS ++++++++++++++++++++++
  #Crop area & efficiency
  crop_area <- crop_data %>%
    select(Season, Crop = c, Area_ha, Eff = Efficiency)

  #Crop Kc coefficient calculation
  crop_Kc <- crop_data %>%select(Season, Crop = c, Jun:May) %>%
    gather(key = Month, value = Kc, -Season, -Crop) %>%
    rowwise() %>%
    mutate(Month = which(Month == month.abb)) %>%
    mutate(Month = as.integer(Month))

  # p1 <- ggplot(demand_ces_mcm, aes(x = Month, y = irr)) +
  #    theme_bw() +
  #    geom_bar(stat="identity", color = "darkred", fill = "darkred") +
  #    labs(x = "Months", y = "MCM per month") +
  #    scale_y_continuous(limits= c(0,8), breaks = seq(0,8, 2)) ; p1
  #
  # p2 <- ggplot(demand_ces_mcm, aes(x = Month, y = dom)) +
  #    theme_bw() +
  #    geom_bar(stat="identity", color = "darkblue", fill = "darkblue") +
  #    labs(x = "Months", y = "MCM per month") +
  #    scale_y_continuous(limits= c(0,8), breaks = seq(0,8, 2)) ; p2

#-------------------------------------------------------------------------------
# SOCIO-ECONOMIC PARAMETERS ----------------------------------------------------

  #From Kenya Gazette (2012), upper limit for m3/houshold, price (USD/m3)
  water_tariff <- data_frame(id = 1:6,
    upper = c(6,20, 50, 100, 300, 1000),
    value = c(0.68,0.78,0.90,1.3,1.6,1.85))


#-------------------------------------------------------------------------------

# PLOTTING PARAMETERS ----------------------------------------------------------

  # Parameters required for climate response surfaces
  alabs <- list(
    x= expression("Temperature change (" * degree * C *")"),
    y= "Precipitation change (%)")
  alim  <- list(temp = c(-0.5,5.5), prec = c(0.65,1.55))
  atick <- list(temp = DelT, prec = DelP, prec2 = (DelP - 1) * 100)

