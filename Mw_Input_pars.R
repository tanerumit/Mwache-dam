
################################################################################
# VARIOUS INPUT PARAMETERS FOR THE MWACHE MODEL
################################################################################


# SCRIPTS & PACKAGES NEEDED ----------------------------------------------------

  options(stringsAsFactors = FALSE)

  #detach("package:dplyr", unload=TRUE)
  #require(Hmisc)
  #require(readr)
  #require(devtools)
  #require(LaplacesDemon)
  #require(mvtnorm)
  #require(truncnorm)
  #require(lazyeval)
  #require(dplyr)

  script_dir <- "/Users/umit/Dropbox/Research/Scripts/"
  source(paste0(script_dir, "HYDROLOGY.R"))
  source(paste0(script_dir, "RESERVOIR.R"))
  source(paste0(script_dir, "STRESS_TEST.R"))
  source(paste0(script_dir, "GGPLOT.R"))

################################################################################
################################################################################

# CLIMATE DATA  ----------------------------------------------------------------

  #Obtained from Brent Boehlert, Industrial Economics

  #145 Model projections: 600 months (50 years) x 4 variables in the order of:
  #1 Observed (based on the Princeton dataset)
  #23 CMIP5 hindcasts
  #56 CMIP3 GCM runs
  #43 CMIP5 GCM runs
  #22 CMIP5 GCM runs processed by Bruce Hewitson's group at UCT
  #(using pressure patterns to drive precip rather
  #than precip directly from the GCMs)

  #Read-in meta-data for the climate projections
  gcm_meta <- read.csv("Inputs/gcm_projections_metadata.csv", header=TRUE)
  hist_mean <- list(Temp = gcm_meta$Mean.Temp[1], Prec = gcm_meta$Mean.Prec[1])

  #Date-Time parameters
  date_gcm  <- seq.Date(as.Date("2000-1-1"),as.Date("2050-12-1"), by="month")
  date_hist <- seq.Date(as.Date("1950-1-1"),as.Date("1999-12-1"), by="month")

  gcm_data <- gcm_meta %>%
    select(Index, Name, Ensemble, Scenario, Temp = Mean.Temp, Prec = Mean.Prec) %>%
    mutate(
      Delta_T  = round(Temp - hist_mean$Temp,3),
      Delta_P  = round(100 * (Prec - hist_mean$Prec)/hist_mean$Prec,3),
      Ensemble = factor(Ensemble,
        levels = c("CMIP3", "CMIP5", "hindcast", "observed")),
      Scenario = factor(Scenario,
        levels =  c("historic", "a1b", "a2", "b1", "rcp45", "rcp85")))

  #Read-in each GCM time-series
  forcings_gcm_list <- lapply(1:145, function(x) read.csv(
      paste0("Inputs/forcings/GCM_runs/", gcm_meta$FileName[x]), header=FALSE)) %>%
    lapply(setNames, nm = c("Year", "Month", "Prec", "Tavg", "Tmin", "Tmax"))

  #GCM times-series (monthly & annual)
  forcings_gcm_monthly <- bind_rows(forcings_gcm_list, .id = "ID")
  forcings_gcm_annual  <- forcings_gcm_monthly %>%
    group_by(ID, Year) %>%
    summarise(Prec=sum(Prec), Tavg=mean(Tavg), Tmin=mean(Tmin), Tmax=mean(Tmax)) %>%
    ungroup() %>%
    arrange(ID, Year)

  #subset 20 GCM runs (select ones with historic emissions)
  #select_index  <- c(2:11,13,15:20,22,23,24,81:100, 102,103,105:116,118:123)

  #Ensemble stats
  gcm_means <- gcm_data %>%
    #gcm_data[select_index,]
    select(Index, Ensemble, Temp, Prec) %>%
    gather(Variable, Value, -Index, -Ensemble) %>%
    group_by(Ensemble,Variable) %>%
    summarise(Mean=mean(Value), Min=min(Value), Max=max(Value)) %>%
    arrange(Variable, Ensemble)

  #Climate Envelopes
  #ensemble_names <- c('CMIP5','FULL')
  #ensembles <- setNames(vector('list', length(ensemble_names)), ensemble_names)
  #ensembles[['CMIP5']] <- filter(gcms_df, Ensemble =='CMIP5') %>%
  #  select(Index, Temp, Prec)
  #ensembles[['FULL']]  <- gcm_data %>%
  #  select(Index, Temp, Prec)

  #Envelope <- lapply(ensembles,
  #    function(data) data[chull(data$Temp, data$Prec),-1]) %>%
  #  lapply(., function(x) rbind(x, x[1,])) %>%
  #  mapply(cbind, ., "Ensemble"=names(ensembles), SIMPLIFY=F) %>%
  #  do.call('rbind', .) %>%
  #  as_data_frame() %>%
  #  mutate(Ensemble = ordered(Ensemble, levels = names(ensembles)))

################################################################################
################################################################################

# STRESS TEST DESIGN -----------------------------------------------------------

  # This section defines the parameters used for the climate stress test
  # Parameters: Annual DeltaT, Annual DeltaP, and climate variability

  #Delta factors for climatic changes (mid-values)
  DeltaT_mid <- seq(0,2.5,0.5)
  DeltaP_mid <- seq(70,140,10)
  T_mid <- hist_mean$Temp + DeltaT_mid
  P_mid <- hist_mean$Prec * DeltaP_mid/100

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

################################################################################
################################################################################

# HYDROLOGY MODEL PARAMETERS -------------------------------------------------

  cmd_to_mcmm <- 365/1e6
  mcmm_to_cmd <- 1e6/365

  #ABCD MODEL PARAMETERS
  par_calib  <- c(0.54988, 293.87242, 0.9, 0.000014) #calibrated parameters

  #WATER RESOURCES MODEL PARAMETERS
  Days <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  monthly_factors <- Days/sum(Days)

  #Drainage area in sq km
  catchment_area  <- 2250

  #Net evaporation (m/month) (NET ET = ET - Precipitation)
  netevap <- c(132, 168, 135, 44, -23, 51, 45, 74, 104, 77, 45, 88) * 1e-3

  Latitude = -4.20


################################################################################
################################################################################

# WR SYSTEMS MODEL PARAMETERS ----------------------------------------------

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
    Total_mcmm = Total_cmd * cmd_to_mcmm,
    Coastal_def = Total_mcmm * 0.59)

  #Monthly irrigation demand
  irrig_monthly <- data_frame(
    months = month.abb,
    demand = c(3.95, 3.17, 2.69, 2.86, 1.93, 3.71, 4.79, 1.94, 1.93, 1.42, 0.82, 4.54))

  #Based on Dam prefeasibilty report (1-129)
  envir_monthly <- data_frame(
    months = as.ordered(factor(month.abb, levels = month.abb)),
    demand = c(0.83, 0, 0, 1.69, 1.49, 0.82, 0.76, 0.55, 0, 0, 0.4, 1.26))

  #p <- ggplot(envir_monthly, aes(x = months, y = demand)) +
  #  theme_bw() +
  #  geom_bar(stat="identity", color = "gray50", fill = "gray50") +
  #  labs(x = "", y = "MCM per month") +
  #  scale_y_continuous(limits= c(0,2), breaks = seq(0,2, 0.5)); p

################################################################################
################################################################################

# GRAPHICAL PARAMETERS ---------------------------------------------------------

#   #Plotting variables
#   alabels <- c(expression(paste("Mean temp. increase (", ~degree~C, ")")),
#                paste("Mean precip. (%)"))
#   aticks <- list(x = seq(0,2.5,0.5), y = seq(70,140,10))
#
#   #Color & shape schemes
#   color_2 <- c("#0072B2", "#D55E00")
#   color_4 <- c('#000000', '#666666', "#D55E00","#0072B2")
#   ColorQual6 <- brewer.pal(6, "Spectral")
#
#   shape_4 <- c(1,0,3,8)
#   shape_2 <- c(1,3,8)
#   shape_2 <- c(3,16)
#
#   #Template
#   pl_surface <- ggplot(mapping=aes(x,y)) +
#     scale_x_continuous(expand = c(0,0), limits = range(aticks$x), breaks = aticks$x) +
#     scale_y_continuous(expand = c(0,0), limits = range(aticks$y), breaks = aticks$y) +
#     labs(x = alabels[1], y = alabels[2]) +
#     theme_bw(base_size = 12) +
#     theme(
#       panel.grid.minor  = element_blank(),
#       panel.grid.major  = element_blank(),
#       plot.title        = element_text(face="bold"),
#       legend.background = element_rect(fill=NA),
#       panel.border      = element_rect(color="black", fill=NA),
#       strip.background  = element_rect(fill="white", color='black'),
#       panel.margin      = unit(1, "lines"))
#

# ################################################################################
