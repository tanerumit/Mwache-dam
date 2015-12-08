
# ------------------------------------------------------------------------------

# PACKAGES & SOURCES -----------------------------------------------------------

  options(stringsAsFactors = FALSE)
  #require(Hmisc)
  #require(readr)
  #require(devtools)
  #require(LaplacesDemon)
  #require(mvtnorm)
  #require(truncnorm)
  #require(lazyeval)
  #require(dplyr)

  script_dir <- "/Users/umit/Dropbox/Research/Scripts/"

  source(paste0(script_dir, "PRIMARY.R"))
  source(paste0(script_dir, "WATER_RESOURCES.R"))
  source(paste0(script_dir, "PLOTTING.R"))

# ------------------------------------------------------------------------------

# CLIMATE DATA  ----------------------------------------------------------------

  #Obtained from Brent Boehlert, Industrial Economics

  #145 Model projections: 600 months (50 years) x 4 variables in the order of:
  #1  Observed (based on the Princeton dataset)
  #23 CMIP5 hindcasts
  #56 CMIP3 GCM runs
  #43 CMIP5 GCM runs
  #22 CMIP5 GCM runs processed by Bruce Hewitson's group at UCT
  #(using pressure patterns to drive precip rather
  #than precip directly from the GCMs)

  #Read-in meta-data for the climate projections
  gcm_meta <- read_csv("Inputs/gcm_metadata.csv", col_names = TRUE)
  historic <- list(Temp = gcm_meta$Temp[1], Prec = gcm_meta$Prec[1])

  #GCM genology
  gcm_groups <- list(
    Grp1  = c("ACCESS1-0", "ACCESS1-3", "HadGEM2-AO","HadGEM2-CC", "HadGEM2-ES"),
    Grp2  = c("GISS-E2-H", "GISS-E2-H-CC", "GISS-E2-R","GISS-E2-R-CC", "giss_aom",
              "giss_model_e_h", "giss_model_e_h", "giss_model_e_r"),
    Grp3  = c("cnrm_cm3", "CNRM-CM5"),
    Grp4  = c("EC-EARTH"),
    Grp5  = c("MIROC5", "MIROC5-ESM", "MIROC5-ESM-CHEM", "MIROC-ESM-CHEM", "MIROC-ESM",
              "MIROC-ESM-CHEM", "miroc3_2_hires", "miroc3_2_medres"),
    Grp6  = c("MRI-CGCM3", "mri_cgcm2_3_2a"),
    Grp7  = c("inmcm3_0"),
    Grp8  = c("FGOALS-g2","FGOALS-s2", "csiro_mk3_5"),
    Grp9  = c("CanESM2"),
    Grp10 = c("CSIRO-Mk3-6-0", "csiro_mk3_0"),
    Grp11 = c("bcc-csm1-1", "bcc-csm1-1-m", "CCSM4","CESM1-BGC", "CESM1-CAM5",
              "FIO-ESM", "NorESM-M", "NorESM-ME", "NorESM1-ME","NorESM1-ME","NorESM1-M"),
    Grp12 = c("MPI-ESM-LR", "MPI-ESM-MR", "CMCC-CM", "CMCC-CESM","CMCC-CMS", "mpi_echam5"),
    Grp13 = c("GFDL-CM3", "GFDL-ESM2G", "GFDL-ESM2M", "gfdl_cm2_0", "gfdl_cm2_1" ),
    Grp14 = c("IPSL-CM5A-LR", "IPSL-CM5A-MR", "IPSL-CM5B-LR", "ipsl_cm4"),
    Grp15 = c("BNU-ESM"))

  gcm_groups_all <- unlist(gcm_groups, use.names = F)
  gcm_meta$Model[!gcm_meta$Model %in% gcm_groups_all]

  CMIP5_meta  <- gcm_meta %>% filter(Scenario %in% c("rcp45", "rcp85"))
  group_index <- sapply(1:nrow(CMIP5_meta), function(y)
    which(sapply(1:15, function(x) CMIP5_meta$Model[y] %in% gcm_groups[[x]])))

  CMIP5_meta <- CMIP5_meta %>% select(Model, Ensemble:Temp) %>%
    mutate(Group = group_index)

  #Date-Time parameters
  date_gcm  <- seq.Date(as.Date("2000-1-1"),as.Date("2049-12-1"), by="month")
  date_hist <- seq.Date(as.Date("1950-1-1"),as.Date("1999-12-1"), by="month")

  #GCM summary table
  ensemble_levels <- c("CMIP3", "CMIP5", "hindcast", "observed")
  scenario_levels <- c("historic", "a1b", "a2", "b1", "rcp45", "rcp85")
  gcm_data <- gcm_meta %>%
    select(Index, Name, Ensemble:Temp) %>%
    mutate(Delta_T  = round(Temp - historic$Temp,3),
           Delta_P  = round(100 * (Prec - historic$Prec)/historic$Prec,3),
           Ensemble = factor(Ensemble, levels = ensemble_levels),
           Scenario = factor(Scenario, levels = scenario_levels))

  #Read-in GCM time-series:
  forcings_gcm <- read_csv("./Inputs/forcings_gcm.csv", skip = 1)

  #GCM times-series (monthly & annual)
  forcings_gcm_monthly <- forcings_gcm
  forcings_gcm_annual  <- forcings_gcm_monthly %>%
    group_by(ID, Year) %>%
    summarise(Prec=sum(Prec), Tavg=mean(Tavg), Tmin=mean(Tmin), Tmax=mean(Tmax)) %>%
    ungroup() %>%
    arrange(ID, Year)

  #subset 20 GCM runs (select ones with historic emissions)
  #select_index <- c(2:11,13,15:20,22,23,24,81:100, 102,103,105:116,118:123)

  #Ensemble stats
  gcm_ensemble_means <- gcm_data %>%
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

  #Monthly coefficients for conversion
  coeff.m <- as.numeric(days_in_month(1:12) / 365)

# ------------------------------------------------------------------------------

# STRESS TEST DESIGN -----------------------------------------------------------

  # This section defines the parameters used for the climate stress test
  # Parameters: Annual DeltaT, Annual DeltaP, and climate variability

  #Delta factors for climatic changes (mid-values)
  nvar_num <- 10
  DelT_mid <- seq(0,2.5,0.5)
  DelP_mid <- seq(0.7,1.4,0.1)
  T_mid <- historic$Temp + DelT_mid
  P_mid <- historic$Prec * DelP_mid

  #Look-up table for climate changes
  cc_tbl  <- expand.grid(
      nTavg = 1:length(DelT_mid), nPrec = 1:length(DelP_mid)) %>%
    mutate(Tavg = DelT_mid[nTavg], Prec = DelP_mid[nPrec], nCC = 1:n()) %>%
    select(nCC, nTavg, nPrec, Tavg, Prec)

  #Look-up table for climate (variability + change)
  clim_tbl <- cc_tbl %>%
    expand.grid.df(., nVar = 1:nvar_num) %>%
    mutate(nClim = 1:n()) %>%
    select(nClim, nCC, nTavg, nPrec, nVar = y, Tavg, Prec)

  iv_index <- sapply(1:nvar_num, function(x) which(clim_tbl$nVar == x))

################################################################################
################################################################################

# ------------------------------------------------------------------------------

# ABCD MODEL PARAMETERS --------------------------------------------------------

  #Misc. conversions
  cmd_to_mcmm <- 365/1e6
  mcmm_to_cmd <- 1e6/365
  month_coef <- as.numeric(days_in_month(1:12)/365)

  #Calibration parameters
  par_calib  <- c(0.54988, 293.87242, 0.9, 0.000014)
  area <- 2250 #in square kilometers
  Latitude = -4.20

# ------------------------------------------------------------------------------

# WR SYSTEMS MODEL PARAMETERS ----------------------------------------------

  BasinLat <- -3.6

  #Net evaporation (m/month) (NET ET = ET - Precipitation)
  netevap <- c(132, 168, 135, 44, -23, 51, 45, 74, 104, 77, 45, 88) * 1e-3

  #RESERVOIR DESIGN ALTERNATIVES (Labels, Volume (MCM), PV_COST (MUSD)
  designs <- data_frame(
    Labels = paste0(seq(40,140,20), "_MCM"),
    Volume = c(40, 60, 80, 100, 120, 140),
    PVCost = c(45.7778, 59.9176, 74.4364, 88.1147, 100.1927, 109.7938))

  #RATING CURVE (From Tahal, 2013)
  rating_curve <- read_csv("./Inputs/reservoir_vol_elev_curve.csv", skip = 1) %>%
    select(elev = Elevation_m, volume = Volume_mcm, area = Area_skm)
  #p <- ggplot(rating_curve, aes(x=volume, y = elev)) + geom_line(); p

  #Domestic demand estimates
  demand_annual_ces <- data_frame(
    Year = c(2015, 2035),
    Total_cmd = c(364000, 887000),
    Total_mcmm = Total_cmd * cmd_to_mcmm,
    Coastal_def = Total_mcmm * 0.59)

  #Monthly irrigation demand estimates by CES, 2014
  demand_ces_mcm <- data_frame(Month = as.factor(1:12),
    irr = c(3.95,3.17,2.69,2.86,1.93,3.71,4.79,1.94,1.93,1.42,0.82,4.54),
    env = c(0.83, 0, 0, 1.69, 1.49, 0.82, 0.76, 0.55, 0, 0, 0.4, 1.26),
    dom = 80.3/12)

   p <- ggplot(demand_ces_mcm, aes(x = Month, y = irr)) +
     theme_bw() +
     geom_bar(stat="identity", color = "darkred", fill = "darkred") +
     labs(x = "Months", y = "MCM per month") +
     scale_y_continuous(limits= c(0,8), breaks = seq(0,8, 2)) ; p


  p <- ggplot(demand_ces_mcm, aes(x = Month, y = dom)) +
     theme_bw() +
     geom_bar(stat="identity", color = "darkblue", fill = "darkblue") +
     labs(x = "Months", y = "MCM per month") +
     scale_y_continuous(limits= c(0,8), breaks = seq(0,8, 2)) ; p




################################################################################
################################################################################

# ------------------------------------------------------------------------------

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

  # alabels <- c(expression(paste("Mean Temp (", ~degree~C, ")")),
  #                              paste("Mean precip. (mm)"))
  #
  # var_col <- colorRampPalette(brewer.pal(10,"Spectral"))(15)

  # ggplot(CMIP5_meta, aes(x = Temp, y = Prec)) +
  #   geom_point(aes(color = as.factor(Group)), size =2) +
  #   scale_y_continuous(expand = c(0,0), limits = range(400, 1400), breaks = seq(400,1400,200)) +
  #   scale_x_continuous(expand = c(0,0), limits = range(25, 28.5), breaks = seq(25,28,1)) +
  #   labs(x = alabels[1], y = alabels[2]) +
  #   theme_bw(base_size = 12) +
  #   scale_color_manual(values = var_col, name="GCM groups") +
  #   stat_ellipse(aes(fill= as.factor(Group)), type = "norm", geom = "polygon", alpha = 0.2)+
  #   geom_text(aes(color = as.factor(Group), label=Model),
  #             position = position_jitter(width=1, height=0),
  #             hjust=0.5, vjust=0.5, size =2.5) +
  #   guides(fill = FALSE) +
  #   facet_grid(~ Scenario, labeller = label_both)

# ------------------------------------------------------------------------------
