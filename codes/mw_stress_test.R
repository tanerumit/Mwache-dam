

#------------------------------------------------------------------------------#
# Title: Climate stress test of the alternative Mwache system designs
# By:    M.Umit Taner
# Date:  August, 2015
# Note:  Implementation of stress test (future conditions)
#        Hydrologic simulations are based on the abcd model
#        Water system operations are based on WEAP
#------------------------------------------------------------------------------#

  library(lubridate)
  library(foreach)
  library(doParallel)

  # required parameters for the analysis
  source("./Codes/Mw_input.R")

# ------------------------------------------------------------------------------
# 1) SAFE YIELD ANALYSIS -------------------------------------------------------

  # Identify safeyield of the system (99% reliability) across all climates
  # No distinction between uses
  # TOTAL RUNS: 8,100 = 3 (datasets) * 54 (cc) * 10 (nvar) * 5 (designs)
  # TOTAL RUNS: 11,340 = 3 (datasets) * 54 (cc) * 10 (nvar) * 7 (designs)

  spin_dataset <- dataset_lst
  spin_K      <- c(80,90,100,110,120,130,140)
  spin_nclim  <- 1:nrow(clim_tbl)

  results_tbl <- expand_grid(dataset=spin_dataset, nClim=spin_nclim, K=spin_K) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, nClim, K)

  #Loop through climate futures
  end <- nrow(results_tbl)

  cl <- makeCluster(6)
  registerDoParallel(cl)
  OUT <- foreach(s=1:end, .combine = rbind, .packages = c("lubridate")) %dopar% {

    #Set parameters based on current iteration
    dataset  <- results_tbl$dataset[s]
    dataset_ind <- which(dataset_lst == dataset)

    c <- results_tbl$nClim[s]
    K <- results_tbl$K[s]

    #abcd calibration parameters
    par_calib <- hydro_pars[[dataset]]

    #climate forcings
    clim <- forcing %>% filter((nclim == c) & (dataset == dataset))

    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

    #Find safe-yield at given target reliability
    OUT <- BINARY_SEARCH(
      FUN = RESERVOIR_RELIABILITY, par = "Tar_dom",
      sign = "+", lower=10, upper=150, target = 95, tol=0.1, max_iter=500,
      beg.y =2020,
      beg.m = 1,
      S_fr=0.8, buffer=1,
      K = K,
      Q = clim$Q_sim, K_d = 20, K_cons = K, K_buff = K *.9,
      evap.m = (clim$PE - clim$Prec) * 1e-3,
      f_elev = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      Tar_irr  = rep(0,12), Tar_eco  = rep(0,12), cycle = FALSE)

    c(dataset_ind, c, K, OUT$objective, OUT$iter, OUT$input_par)

  }

  stopCluster(cl)

  df <- OUT %>% as.data.frame() %>%
    as_data_frame() %>%
    mutate(Dataset = dataset_lst[V1]) %>%
    select(Dataset, nClim = V2, K = V3, REL = V4, iter = V5, safeyield = V6)

  write_csv(x = df, path = "stresstest_safeyield.csv")

# ------------------------------------------------------------------------------
# 2) RELIABILTY ANALYSIS  ------------------------------------------------------

  # Calculte system reliability at different demand increase scenarios
  # TOTAL RUNS: 68040 = 3 * 540 * 7 * 6

  spin_dataset <- c("Princeton","GPCC","CRU")
  spin_K <- c(80,90,100,110,120,130,140)
  spin_nclim <- 1:nrow(clim_tbl)
  spin_deml = c(57,68,80,91,103,114)

  #DATA structure
  results_tbl <- expand_grid(dataset = spin_dataset, nclim = spin_nclim,
      size = spin_K, demand = spin_deml) %>%
    mutate(rel = NA, crel = NA) %>%
    left_join(clim_tbl, by = "nclim") %>%
    select(dataset, size, demand, nclim, nvar, temp, prec, rel, crel) %>%
    arrange(dataset, size, demand, nclim)

  end <- nrow(results_tbl)

  #Loop through future conditions
  cl <- makeCluster(6)
  registerDoParallel(cl)
  OUT <- foreach(s=1:end, .combine = rbind,
                 .packages = c("lubridate")) %dopar% {

    # Set dataset, climate (c) and design (d), demand target
    dataset <- results_tbl$dataset[s] %>% as.character()

    c <- results_tbl[s,"nclim"] %>% as.numeric()
    K <- results_tbl[s,"size"] %>% as.numeric()
    demand_dom <- results_tbl[s,"demand"] %>% as.numeric()

    #abcd calibration parameters
    par_calib <- hydro_pars[[dataset]]

    #climate forcings
    clim <- forcing %>% filter((nclim == c) & (dataset == dataset))

    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3
    clim$Q_sim <- as.numeric(clim$Q_sim)

    # Simulate reservoir operations
    sim <- RESERVOIR_SIM(
      beg.y    = 2020,
      beg.m    = 1,
      Q        = clim$Q_sim,
      K_d      = 20,
      K        = K,
      K_cons   = K,
      K_buff   = 70,
      buffer   = 1,
      Tar_dom  = demand_dom,
      Tar_irr  = 0,
      Tar_eco  = 0,
      evap.m = (clim$PE - clim$Prec) * 1e-3,
      f_elev = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      cycle = FALSE)


    #Calculate volumetric reliablity
    demand <- sim$Tar_dom
    supply <- sim$R_dom

    nofail <- length(which(demand-supply == 0))
    REL  <- 100 * nofail / length(demand)

    nofail <- length(which(supply > critical_level * demand))
    cREL <- 100 * nofail / length(demand)

    c(s, REL, cREL)

  }

  stopCluster(cl)
  df <- results_tbl %>%
    mutate(rel = OUT[,2], crel = OUT[,3])

  write_csv(x = df, './Inputs/stresstest_reliability.csv')


# ------------------------------------------------------------------------------

  ### MERGE STRESS TEST ANALYSIS RESULTS
  dat1 <- read_csv("./data/stresstest_safeyield_old.csv") %>%
    rename(size = K, dataset = Dataset, nclim = nClim, rel = REL)

  dat2 <- read_csv("./data/stresstest_reliability_old.csv")

  dat <- dat2 %>%
    left_join(select(dat1, -rel, -iter), by = c("dataset", "size", "nclim"))

  write_csv(dat, "./data/stresstest_ffd2.csv")

# ------------------------------------------------------------------------------

################################################################################
################################################################################
################################################################################

  #LOAD-DATA FOR PLOTTING
  stest_data <- read_csv("./data/stresstest_ffd.csv")
  CMIP5_data <- read_csv("./data/cmip5_deltaclim.csv")
  CMIP5_gen  <- read_csv("./data/cmip5_genealogy.csv")

# RESPONSE SURFACES ------------------------------------------------------------

  # Parameters required for climate response surfaces
  alabs <- list(
    x= expression("Temperature change (" * degree * C *")"),
    y= "Precipitation change (%)")
  alim  <- list(temp = c(-0.5,5.5), prec = c(0.65,1.55))
  atick <- list(temp = DelT, prec = DelP)

  theme_set(theme_bw(base_size = 11))
  theme_update(
    plot.title        = element_text(face = "bold"),
    panel.grid.minor  = element_blank(),
    panel.grid.major  = element_blank(),
    panel.border      = element_rect(color = "black", fill=NA),
    panel.margin      = unit(1.5, "lines"),
    strip.background  = element_rect(fill = "white", color='black'),
    legend.background = element_rect(fill = NA),
    legend.text       = element_text(size=9),
    legend.title      = element_text(size=9)
  )

 # GCM PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 p1 <- ggplot(CMIP5_df, aes(x=temp, y=prec, fill = GCM, shape = scenario)) +
   theme_bw(base_size = 12) +
   geom_point(size = 2.5, stroke = 1) +
   labs(x = alabs$x, y = alabs$y, shape = "Scenarios") +
   scale_fill_manual(limits = CMIP5_groups$GCM, values = CMIP5_groups$color) +
   geom_vline(xintercept = 0, linetype = "dashed") +
   geom_hline(yintercept = 1, linetype = "dashed") +
   scale_shape_manual(name = "Scenario", values = c(21,21,22,24)) +
   guides(fill = guide_legend(
     override.aes = list(shape = 21,linetype = "dashed"), order = 1, ncol = 2, keyheight = 0.7),
     shape  = guide_legend(order = 2, ncol = 1, keyheight = 0.7)) +
   labs(x = alabs$x, y = alabs$y, shape = "Scenarios") +
   scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
   scale_y_continuous(expand=c(0,0), limits = alim$prec, breaks = atick$prec,
                      labels = seq(-30,50,10))

  ggsave(filename = 'gcm_scatter_all.tiff', p1, height = 5, width = 8)

  CMIP5_df_means <- CMIP5_df %>%
    na.omit() %>%
    group_by(Family, scenario) %>%
    summarize(temp = mean(temp), prec = mean(prec)) %>%
    ungroup() %>%
    mutate(Family = as.factor(Family))


  p2 <- ggplot(CMIP5_df_means, aes(x=temp, y=prec, fill = Family, shape = scenario)) +
    theme_bw(base_size = 12) +
    geom_point(size = 2.5, stroke = 1) +
    labs(x = alabs$x, y = alabs$y) +
    scale_shape_manual(name = "Scenario", values = c(21,21,22,24)) +
    scale_fill_manual(limits = CMIP5_groups$Family, values = CMIP5_groups$color) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = alim$prec, breaks = atick$prec,
                       labels = seq(-30,50,10)) +
    guides(fill = FALSE) +
    stat_ellipse(aes(group = 1), level = 0.99)

  ggsave(filename = 'gcm_scatter_means.tiff', p2, height = 6, width = 7)

  # RELIABILITY PLOTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #Define bin & color scheme
  bin1 <- seq(45,95,10)
  bin2 <- seq(96,100,1)
  col1 <- colorRampPalette(c("red", "white"))(length(bin1))
  col2 <- colorRampPalette(c("lightsteelblue1", "blue"))(length(bin2))
  var_col <- c(col1,col2)
  var_bin <- c(bin1, bin2)

  #Parameters
  dataset <- "Princeton"

  #~~~~~~~ Plot a single surface
  data_st <- stest_data %>%
    filter(dataset == dataset, size == 120) %>%
    group_by(demand, temp, prec) %>%
    #rename(value = rel) %>%
    rename(value = crel) %>%
    summarize(value = mean(value)) %>%
    mutate(Bins = cut(value, breaks = var_bin, dig.lab = 5, include.lowest = T)) %>%
    ungroup()

  p <- ggplot(filter(data_st, demand == 68), aes(temp,prec)) +
    geom_tile(aes(fill = Bins), color = "gray80") +
    scale_fill_manual(values=var_col, name="Reliability (%)", drop = F)  +
    guides(fill = guide_legend(reverse = TRUE , keyheight = 1)) +
    labs(x = alabs$x, y = alabs$y) +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = alim$prec, breaks = atick$prec,
      labels = seq(-30,50,10))

  #ggsave(paste0("surfacepl_crel_",dataset,".tiff"), plot=p, height=5, width=6)

  #~~~~~~~ Plot multiple surfaces (facet)
  data_st2 <- data_st %>%
    mutate(demand = factor(demand, levels = c(57,68,91,114),
                           labels = c("50% (57 MCM)",
                                      "80% (68 MCM)",
                                      "140% (91 MCM)",
                                      "200% (114 MCM)")))


  p <- ggplot(mapping = aes(temp,prec)) +
    theme_bw(base_size = 10) +
    geom_tile(aes(fill = Bins), data = data_st2, color = "gray80") +
    facet_wrap(~ demand, switch = NULL) +
    scale_fill_manual(values=var_col, name="Critical\nreliability (%)", drop = F) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(reverse = TRUE , keyheight = 1)) +
    labs(x = alabs$x, y = alabs$y) +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = alim$prec, breaks = atick$prec,
       labels = seq(-30,50,10))

  ggsave(paste0("surfacepl_crel_perDem_",dataset,".tiff"), plot=p, height=7, width=6)

  # SAFEYIELD SURFACE PLOTS ++++++++++++++++++++++++++++++++++++++++++++++++++++

  #Define color scheme
  bin1 <- seq(60,70,5)
  bin2 <- seq(80,120,10)
  col1 <- colorRampPalette(c("red", "white"))(length(bin1))
  col2 <- colorRampPalette(c("lightsteelblue1", "blue"))(length(bin2))
  var_col <- c(col1,col2)
  var_bin <- c(bin1, bin2)

  dataset <- "Princeton"

  data_st <- stest_data %>%
    filter(dataset == dataset, size == 120) %>%
    group_by(demand, temp, prec) %>%
    rename(value = safeyield) %>%
    summarize(value = mean(value)) %>%
    mutate(Bins = cut(value, breaks = var_bin, dig.lab = 5, include.lowest = T))

  p <- ggplot(data_st, aes(temp,prec)) +
    geom_tile(aes(fill = Bins), color = "gray20") +
    scale_fill_manual(values=var_col, name="Safe yield\n(MCM)", drop = F)  +
    guides(fill = guide_legend(reverse = TRUE , keyheight = 1)) +
    labs(x = alabs$x, y = alabs$y) +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = alim$Prec, breaks = atick$prec,
      labels = seq(-30,50,10))

  ggsave(height = 5, width = 6,
      filename = paste0("surface_syield_",dataset,".tiff"))

# ------------------------------------------------------------------------------
# PLOT EMPRICAL CDFs -----------------------------------------------------------

  df2 <- data_frame(
    scenario = c("low", "medium", "high"),
    value = c(68, 91 ,114)) %>%
    mutate(scenario = factor(scenario,
      levels = c("low", "medium", "high")))

  safeyield_dat <- stest_data %>% filter(demand == 91)

  df <- safeyield_dat %>%
    filter(size == 120) %>%
    rename(value = safeyield) %>%
    arrange(desc(value)) %>%
    mutate(ExceedP = (1:n())/(n()+1))

  p <- ggplot(cdf_df, aes(x= ExceedP, y = value)) +
    theme_bw(base_size = 11) +
    #geom_line(size =1, color = "black") +
    geom_smooth(size = 1, se = FALSE, color = "black") +
    geom_hline(yintercept = 80.3, linetype = "dashed") +
    labs(x = "Probability of Exceedance", y = "Safeyield (MCM)",
      color = "Demand levels") +
    scale_y_continuous(limits = c(50,130), breaks = seq(50,130,20)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))

  ggsave(height = 4, width = 7, plot = p, filename = paste0("CDF_syield.png"))

  df %<>% mutate(
    Demand = factor(Demand, levels= c(1,1.5,2, 2.5),
    labels = c("+0%", "+50%", "+100%", "150%")))

  p1 <- ggplot(df, aes(x= ExceedP, y = value, color = Demand)) +
    theme_bw(base_size = 11) +
    geom_line(size = 1) +
    labs(x = "Probability of Exceedance", y = "Reliability (%)") +
    geom_hline(yintercept = 99, linetype = "dashed") +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1))

  #### CRITICAL RELIABILITY METRIC
  REL_dat <- stest_data

  df <- REL_dat %>% filter(size == 120) %>%
    rename(value = rel) %>%
    mutate(demand = as.factor(demand))

  cdf_df <- df %>%
    group_by(demand) %>%
    arrange(desc(value)) %>%
    mutate(ExceedP = (1:n())/(n()+1)) %>%
    ungroup()

  cdf_df %<>%
    mutate(demand = factor(demand, levels = c(57,68,91,114),
                           labels = c("50% (57 MCM)",
                                      "80% (68 MCM)",
                                      "140% (91 MCM)",
                                      "200% (114 MCM)")))


  p2 <- ggplot(cdf_df, aes(x= ExceedP, y = value, color = demand)) +
    theme_bw(base_size = 11) +
    geom_line(size = 1) +
    labs(x = "Probability of Exceedance", y = "Reliability (%)") +
    geom_hline(yintercept = 95, linetype = "dashed") +
    scale_y_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
    labs(color = "Demand increases")

  ggsave("exceedP_rel.tiff", p2, width = 7, height = 4)


  ### CDF of safeyield
  df <- syield_dat  %>%
    rename(value = safeyield) %>%
    mutate(Capacity = as.factor(K)) %>%
    select(Capacity, value)

  cdf_df <- df %>%
    group_by(Capacity) %>%
    arrange(desc(value)) %>%
    mutate(ExceedP = (1:n())/(n()+1)) %>%
    ungroup() %>%
    mutate(Capacity = factor(Capacity, levels= c(80, 100, 120, 130, 140)))


  p1 <- ggplot(cdf_df, aes(x= ExceedP, y = value, color = Capacity)) +
    theme_bw(base_size = 11) +
    geom_smooth(size = 1) +
    labs(x = "Probability of Exceedance", y = "safeyield (MCM)") +
    geom_hline(yintercept = 68, linetype = "dashed") +
    scale_y_continuous(limits = c(20,140), breaks = seq(20,140,20), expand = c(0,0)) +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.1), expand = c(0,0)) +
    scale_color_brewer(name = "Capacity\n(MCM)", palette = "Set1")

  #ggsave(plot = p1, filename = "yield_designs.png", height = 4, width = 7)

  # Median yield vs cost
  df <- cdf_df %>%
    group_by(Capacity) %>%
    summarize(
      q50 = median(value),
      q10 = quantile(value, 0.1),
      q25 = quantile(value, 0.25),
      q75 = quantile(value, 0.75),
      q90 = quantile(value, 0.9)) %>%
    mutate(Volume = as.numeric(as.character(Capacity))) %>%
    left_join(designs, by = "Volume") %>%
    select(-Volume, -Labels) %>%
    gather(key = Conditions, value = value, -Capacity, -PVCost) %>%
    mutate(Conditions = factor(Conditions,
      levels = c('q90', 'q75', 'q50', 'q25', 'q10'),
      labels = c('Very wet', 'Wet', 'Median', 'Dry', 'Very dry')))

  p <- ggplot(df, aes(x = PVCost, y = value, group = Conditions)) +
    theme_bw(base_size = 11) +
    geom_point(aes(color = Capacity), size = 3, shape = 15) +
    geom_line(aes(linetype = Conditions),color = "black", size = 1) +
    labs(x = "Project capital cost (M.USD)", y = "safeyield (MCM)") +
    scale_y_continuous(limits = c(20,140), breaks = seq(20,140,20), expand = c(0,0)) +
    scale_x_reverse() +
    scale_color_brewer(name = "Capacity\n(MCM)", palette = "Set1")

  ggsave(plot = p, filename = "yield_pvcost.png", height = 5, width = 6)


# ------------------------------------------------------------------------------




