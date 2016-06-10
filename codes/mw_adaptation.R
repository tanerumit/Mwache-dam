
#------------------------------------------------------------------------------#
# Title: Exploring alternative adaptation options for Mwache
# By:    M.Umit Taner
# Date:  April, 2016
#------------------------------------------------------------------------------#

  library(foreach)
  library(doParallel)
  library(lubridate)

# required parameters for the analysis
  source("./Codes/Mw_input.R")

# GCM-based likelihoods
  nclim_select <- clim_tbl %>%
    filter(temp <= 3 && prec >= 0.6 && prec <= 1.5) %$% nclim

# 1) ALLOCATION ANALYSIS -------------------------------------------------------


# SIMULATION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  # TOTAL RUNS: 540 cli x 1 design x 1 demand levels x 6 rules =

  # DEMAND LEVEL (MEDIUM INCREASE 140% - 91 MCM in 2035)
  demand_yr <- data_frame(year = c(2015,2035), value = c(38, 91)) %>%
    do(Demand = as.data.frame(approx(x=.$year, .$value, xout=2015:2035))) %>%
    unnest(Demand) %>%
    select(year = x, value = y) %>%
    bind_rows(data_frame(year = 2036:2069, value = 91)) %>%
    filter(year %in% 2020:2069)

  demand_mon <- demand_yr %>%
    expand_grid_df(data.frame(month = 1:12)) %>% as.tbl() %>%
    select(year, month, value_yr = value) %>%
    arrange(year, month) %>%
    mutate(value = value_yr * month_coef[month]) %>%
    select(year, month, value)

  span_dataset <- "Princeton"
  span_K <- 120
  span_clim = 1:nrow(clim_tbl)
  span_prio = 1:3
  span_alpha = c(0.5,0.9)

  prio_lst <- list(
    c("T_eco", "T_dom", "T_irr"),   #EDI
    c("T_dom", "T_eco", "T_irr"),   #DEI
    c("T_dom", "T_irr", "T_eco"))   #DIE


  results_tbl <- expand_grid(
    dataset = span_dataset, K = span_K,
    nclim = span_clim, prio = span_prio, alpha = span_alpha) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, K, nclim, prio, alpha) %>%
    mutate(Id = 1:n())

  #Loop through climate futures
  end <- nrow(results_tbl)

  ptm <- proc.time()
  cl <- makeCluster(6)
  registerDoParallel(cl)

  OUT <- foreach(s=1:end, .combine = rbind,
    .packages = c("lubridate", "sirad")) %dopar% {

    #Set parameters based on current iteration
    dataset  <- "Princeton"
    par_calib <- hydro_pars[[dataset]]
    c <- results_tbl$nclim[s]
    K <- results_tbl$K[s]
    prio <- prio_lst[[results_tbl$prio[s]]]
    alpha = results_tbl$alpha[[s]]

    #climate forcings
    clim <- forcing %>% filter((nclim == c) & (dataset == "Princeton"))
    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3


    #### CALCULATE IRRIGATION DEMAND ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate ET0 (Penmann - Monteith equation (FAO))
    irrig_cal <- clim %>%
      left_join(crop_pars_et0, by = "Month") %>%
      mutate(
        EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
        Wind  = Wind * 0.01157407,
        v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
        vap_p = Humid/100 * v_sat,
        ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
          tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
          days  = Date, lat = Latitude) * as.numeric(days_in_month(Month))) %>%
      select(Year, Month, Prec, Tavg, Tmin, Tdif, PE, EPrec, ET0)

    #Calculate CRW (mm/month)
    irr_demand <- irrig_cal %>%
      left_join(crop_Kc, by = "Month") %>%
      left_join(crop_area, by = c("Season","Crop")) %>%
      mutate(
        Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
      mutate(
        ETc = ET0 * Kc,
        NIR_mm = pmax(ETc - EPrec, 0),
        GIR_mm = NIR_mm / (Eff/100),
        GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
      mutate(Month = as.integer(Month)) %>%
      group_by(Year, Month) %>%
      summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm

    # Simulate reservoir operations
    sim <- reservoir_sim(
      Q        = clim$Q_sim,
      K        = 120,
      K_dead   = 20,
      T_dom    = demand_mon$value,
      T_irr    = irr_demand,
      T_eco    = demand_ces_mcm$env,
      evap.m   = (clim$PE - clim$Prec) * 1e-3,
      f_elev   = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      C_alloc  = alpha,
      priority = prio)

    #Calculate reliabilities
    nofail <- length(which(sim$R_dom >= 1 * sim$T_dom))
    dom_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= 1 * sim$T_irr))
    irr_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= 1 * sim$T_eco))
    eco_rel  <- 100 * nofail / nrow(sim)

    #Calculate critical reliabilities
    dom_crel <- 100 * sum(sim$R_dom) / sum(sim$T_dom)
    eco_crel <- 100 * sum(sim$R_eco) / sum(sim$T_eco)
    irr_crel <- 100 * sum(sim$R_irr) / sum(sim$T_irr)

    #nofail <- length(which(sim$R_dom >= alpha * sim$T_dom))
    #dom_crel  <- 100 * nofail / nrow(sim)
    #nofail <- length(which(sim$R_irr >= alpha * sim$T_irr))
    #irr_crel  <- 100 * nofail / nrow(sim)
    #nofail <- length(which(sim$R_eco >= alpha * sim$T_eco))
    #eco_crel  <- 100 * nofail / nrow(sim)

    c(s, dom_rel, irr_rel, eco_rel, dom_crel, irr_crel, eco_crel)
  }
  stopCluster(cl)

  df <- OUT %>% as.data.frame() %>% as_data_frame() %>%
    select(Id = V1, rel_dom = V2, rel_irr = V3, rel_eco = V4,
           crel_dom = V5, crel_irr = V6, crel_eco = V7) %>%
    left_join(results_tbl, ., by = "Id") %>%
    select(Id, dataset:alpha, rel_dom:crel_eco)

  #write_csv(x = df, path = "./data/adaptation_prio.csv")


#EVALUATE RESULTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  prio_df <- read_csv('./data/adaptation_prio.csv') #%>%
    #filter(nClim %in% nClim_select)

  #checks
  temp <- prio_df %>%
    group_by(prio, alpha) %>%
    summarize(rel_dom = mean(rel_dom), rel_irr = mean(rel_irr),
      rel_eco = mean(crel_eco), crel_dom = mean(crel_dom),
      crel_irr = mean(crel_irr), crel_eco = mean(crel_eco))

  df <- prio_df %>%
   rowwise() %>%
    mutate(Rule = switch(as.character(prio),"1" ="EDI", "2"="DEI", "3"="DIE"),
      Rule = paste(Rule, alpha, sep = "-")) %>%
    select(Id, nclim, Rule, crel_dom, crel_irr, crel_eco) %>%
    gather(key = variable, value = value, -Id, -nclim, -Rule) %>%
    mutate(variable = factor(variable,
      levels = c("crel_eco", "crel_dom", "crel_irr"),
      labels = c("Environ.", "Domestic", "Irrigation")))


  #Box-plot of reliabilities
  p <- ggplot(df, aes(x = variable, y = value, color = variable)) +
    theme_bw(base_size = 10) +
    facet_wrap( ~ Rule, labeller = label_both, ncol = 2, scales = "free") +
    geom_boxplot() +
    stat_summary(fun.y=mean, geom="point", size=2, shape = 16, stroke = 1) +
    labs(x = "", y = "Reliability (%)") +
    guides(color = FALSE) +
    scale_y_continuous(limits = c(50,100), breaks = seq(50,100,25))

  #ggsave(plot = p, filename = 'alloc_boxplot.tiff', height = 5, width = 7)

  df <- prio_df %>% filter(nclim == 244) %>%
    rowwise() %>%
    mutate(Rule = switch(as.character(prio),"1" ="EDI", "2"="DEI", "3"="DIE")) %>%
    ungroup() %>%
    select(Id, nclim, Rule, crel_dom, crel_irr, crel_eco, alpha) %>%
    gather(key = variable, value = value, crel_dom:crel_eco) %>%
    mutate(variable = factor(variable,
      levels = c("crel_eco", "crel_dom", "crel_irr"),
      labels = c("Environ.", "Domestic", "Irrigation")))

  p <- ggplot(df, aes(x=variable, y = value, group = Rule, color = Rule)) +
    theme_bw() +
    facet_wrap( ~ alpha, labeller = label_both, ncol = 2) +
    geom_point(size = 1) +
    geom_line() +
    labs(x = "", y = "Reliability (%)") +
    guides(color = FALSE) +
    scale_y_continuous(limits = c(80,100), breaks = seq(80,100,10))


  ggsave(plot = p, filename = 'alloc_median.tiff', height = 2.5, width = 7)


  df <- data %>%
    select(Id:demL, Domestic = cRel_dom, Irrig = cRel_irr, Envir = cRel_eco) %>%
    gather(key = variable, value = value, Domestic:Envir) %>%
    mutate(
      variable = factor(variable, levels = c("Envir", "Domestic", "Irrig")),
      K = factor(K, levels = c(120, 130),labels = c("120 MCM", "130 MCM")),
      demL = factor(demL, levels=1:3, labels = c("+60%", "+120%", "+180%"))) %>%
    select(Capacity = K, nClim, Demand = demL, variable, value, Rule = Prio)

  p <- ggplot(df, aes(x = variable, y = value, color = variable)) +
    facet_wrap( ~ Rule, labeller = label_both, nrow = 1) +
    theme_bw(base_size = 10) +
    geom_boxplot() +
    #stat_summary(fun.y=mean, geom="point", size=2, shape = 21, stroke = 1) +
    labs(x = "", y = "Critical Reliability (%)") +
    scale_y_continuous(limits = c(20, 100), breaks = seq(20,100,20)) +
    guides(color = FALSE) +
    scale_color_manual(values = var_col)

    #ggsave(plot = p, filename = 'crel_prio.png', width = 7, height = 4)

  nClim_q <- c(Q10, Q25, Q50, Q75, Q90) #525 583 637 316 385

  p <- list()

  for (i in 1:3) {

    df2 <- df %>%
      filter(nClim %in% nClim_q) %>%
      spread(key = variable, value = value) %>%
      mutate(ID = 1:n()) %>%
      filter(Capacity == "120 MCM") %>%
      filter(Demand == "+120%") %>%
      gather(key = variable, value = value, Envir:Irrig) %>%
      filter(Rule == i) %>%
      mutate(variable = factor(variable, levels = prio_labels[[i]]))

    p[[i]] <- ggplot(df2, aes(x=variable, y = value, group = ID)) +
      facet_grid( ~ Rule) +
      theme_bw() +
      geom_line(aes(linetype = as.factor(nClim)))

  }

  df2 <- df %>%
    filter(nClim %in% nClim_q) %>%
    spread(key = variable, value = value) %>%
    #mutate(ID = 1:n()) %>%
    #filter(Capacity == "120 MCM") %>%
    #filter(Demand == "+120%") %>%
    gather(key = variable, value = value, Envir:Irrig) %>%
    spread(key = variable, value = value)

  ggplot(data = df2, aes(x = Domestic, y = Envir)) +
    theme_bw() +
    geom_point()


  nClim_q

  # QUANTILES OF STREAMFLOW ++++++++++++++++++++++++++++++++++

  forcing_i <- forcing %>% filter(dataset == "Princeton")
  Qmean <- vector(mode = "numeric", length = 540)


  for (i in 1:540) {

    clim <- filter(forcing_i, nclim == i)

    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

    Qmean[i] <- mean(clim$Q_sim)
  }

  Qmean_df <- data_frame(nclim = 1:540, value = Qmean) %>%
    arrange(desc(value)) %>%
    mutate(ExceedP = (1:n())/(n()+1))

  ExP <- Qmean_df$ExceedP

  Q10 <- Qmean_df[['nclim']][which.min(abs(0.10 - ExP))]
  Q25 <- Qmean_df[['nclim']][which.min(abs(0.25 - ExP))]
  Q50 <- Qmean_df[['nclim']][which.min(abs(0.50 - ExP))]
  Q75 <- Qmean_df[['nclim']][which.min(abs(0.75 - ExP))]
  Q90 <- Qmean_df[['nclim']][which.min(abs(0.90 - ExP))]

  nClim_q <- c(Q10, Q25, Q50, Q75, Q90)

  nClim_q <- c(373, 150, 244, 338,  11)

  filter(Qmean_df, nclim %in% nClim_q)

  out[c(1,4,3,2,5), c(1,2,4,3,5)]

# 2) BUFFER COEFFICIENT ANALYSIS -----------------------------------------------

  span_dataset <- "Princeton"
  span_K <- c(120)
  span_Clim = nClim_q
  span_demL = 1:3
  span_buffer <- c(0.3,0.5,0.7,0.8,0.9,1)

  results_tbl <- expand_grid(dataset = span_dataset, K = span_K,
      nClim = span_Clim, Buffer = span_buffer, demL = span_demL) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, K, nClim, demL, Buffer)

  end <- nrow(results_tbl)

  ptm <- proc.time()
  cl <- makeCluster(1)
  registerDoParallel(cl)

  OUT <- foreach(
    s=1:end, .combine = rbind, .packages = c("lubridate", "sirad")) %dopar% {

    #Set parameters based on current iteration
    dataset  <- results_tbl$dataset[s]
    par_calib <- hydro_pars[,dataset]

    c <- results_tbl$nClim[s]
    K <- results_tbl$K[s]
    bufferC <- results_tbl$Buffer[s]
    dom_demand <- filter(demand_levels, level == results_tbl$demL[s]) %$% value

    #climate forcings
    clim <- forcing %>% filter((nClim == c) & (Dataset == dataset))
    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

    #### CALCULATE IRRIGATION DEMAND ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate ET0 (Penmann - Monteith equation (FAO))
    irrig_cal <- clim %>%
      left_join(crop_pars_et0, by = "Month") %>%
      mutate(
        EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
        Wind  = Wind * 0.01157407,
        v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
        vap_p = Humid/100 * v_sat,
        ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
        tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
          days  = Date, lat = Latitude) * as.numeric(days_in_month(Month))) %>%
          select(Year, Month, Prec, Tavg, Tmin, Tdif, PE, EPrec, ET0)

    #Calculate CRW (mm/month)
    irr_demand <- irrig_cal %>%
      left_join(crop_Kc, by = "Month") %>%
      left_join(crop_area, by = c("Season","Crop")) %>%
    mutate(
      Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
    mutate(
      ETc = ET0 * Kc,
      NIR_mm = pmax(ETc - EPrec, 0),
      GIR_mm = NIR_mm / (Eff/100),
      GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
      mutate(Month = as.integer(Month)) %>%
      group_by(Year, Month) %>%
      summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm

    # Simulate reservoir operations
    sim <- RESERVOIR_SIM(
      beg.y    = 2020,
      beg.m    = 1,
      Q        = clim$Q_sim,
      K_d      = 20,
      K        = K,
      K_cons   = K,
      K_buff   = 70,
      buffer   = bufferC,
      T_dom  = dom_demand,
      T_irr  = irr_demand,
      T_eco  = demand_ces_mcm$env,
      evap.m = (clim$PE - clim$Prec) * 1e-3,
      f_elev = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      cycle = FALSE)

    #Calculate reliabilities
    nofail <- length(which(sim$R_dom >= 1 * sim$T_dom))
    dom_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= 1 * sim$T_irr))
    irr_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= 1 * sim$T_eco))
    eco_rel  <- 100 * nofail / nrow(sim)

    #Calculate critical reliabilities
    nofail <- length(which(sim$R_dom >= critical_level * sim$T_dom))
    dom_crel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= critical_level * sim$T_irr))
    irr_crel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= critical_level * sim$T_eco))
    eco_crel  <- 100 * nofail / nrow(sim)

    c(s, dom_rel, irr_rel, eco_rel, dom_crel, irr_crel, eco_crel)
  }
  stopCluster(cl)

  results_tbl %<>% mutate(Id = 1:n())
  df <- OUT %>% as.data.frame() %>%
    as_data_frame() %>%
    select(Id = V1, Rel_dom = V2, Rel_irr = V3, Rel_eco = V4,
           cRel_dom = V5, cRel_irr = V6, cRel_eco = V7) %>%
    left_join(results_tbl, ., by = "Id") %>%
    select(Id, dataset:demL, Rel_dom:cRel_eco)
  #write_csv(x = df, path = "./inputs/adaptation_buffer.csv")

  #### PLOTS FOR BUFFER COEFFICIENT EXPLORING

  buffer_dat <- df
  #buffer_dat <- read_csv("./inputs/adaptation_buffer.csv")

  data <- buffer_dat %>%
    filter((demL == 1)) %>%
    filter(nClim %in% nClim_q)

  p <- ggplot(data, aes(x = as.factor(Buffer), color = as.factor(nClim))) +
    theme_bw(base_size = 11) +
    scale_y_continuous(limits = c(40,100), breaks = seq(40, 100, 20)) +
    scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3',
      '#ff7f00'), name="Condition", breaks=c("525", "583", "637", "316", "385"),
      labels=c("Very wet", "Wet", "Median", "Dry", "Very Dry"))

  p1 <-  p +  geom_point(aes(y = Rel_dom), size = 1, stroke = 1, shape = 21) +
    geom_line(aes(y = Rel_dom, group = as.factor(nClim))) +
    labs(x = "Buffer Coefficient", y = "Reliability (%)")

  p2 <-  p +  geom_point(aes(y = cRel_dom), size = 1, stroke = 1, shape = 21) +
    geom_line(aes(y = cRel_dom, group = as.factor(nClim))) +
    labs(x = "Buffer Coefficient", y = "Critical Reliability (%)")

  ggsave(filename = "buffer_rel.png", p1, height = 4, width = 7)
  ggsave(filename = "buffer_crel.png", p2, height = 4, width = 7)

  #Q10, Q25, Q50, Q75, Q90
  #525 583 637 316 385






# 3) DEMAND REDUCTION ANALYSIS -------------------------------------------------


  span_dataset <- "Princeton"
  span_K <- 120
  span_Clim = c(373, 150, 244, 338, 111)
  span_reduc <- seq(0.4,1,0.1)

  results_tbl <- expand_grid(dataset = span_dataset, K = span_K,
      nclim = span_Clim, reduction = span_reduc) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, K, nclim, reduction)

  end <- nrow(results_tbl)
  out <- data_frame(s = 1:end, nclim = NA, reduction = NA,
    Qmean = NA, dom_rel = NA, eco_rel = NA, irr_rel = NA)

  for (s in 1:end) {

    #Set parameters based on current iteration
    dataset  <- "Princeton"
    par_calib <- hydro_pars[[dataset]]
    c <- results_tbl$nclim[s]
    K <- results_tbl$K[s]

    dem_reduce <- results_tbl$reduction[s]

    #climate forcings
    clim <- forcing %>% filter((nclim == c) & (dataset == "Princeton"))
    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

    #### CALCULATE IRRIGATION DEMAND ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate ET0 (Penmann - Monteith equation (FAO))
    irrig_cal <- clim %>%
      left_join(crop_pars_et0, by = "Month") %>%
      mutate(
        EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
        Wind  = Wind * 0.01157407,
        v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
        vap_p = Humid/100 * v_sat,
        ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
          tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
          days  = Date, lat = Latitude) * as.numeric(days_in_month(Month))) %>%
      select(Year, Month, Prec, Tavg, Tmin, Tdif, PE, EPrec, ET0)

    #Calculate CRW (mm/month)
    irr_demand <- irrig_cal %>%
      left_join(crop_Kc, by = "Month") %>%
      left_join(crop_area, by = c("Season","Crop")) %>%
      mutate(
        Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
      mutate(
        ETc = ET0 * Kc,
        NIR_mm = pmax(ETc - EPrec, 0),
        GIR_mm = NIR_mm / (Eff/100),
        GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
      mutate(Month = as.integer(Month)) %>%
      group_by(Year, Month) %>%
      summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm

    # Simulate reservoir operations
    sim <- RESERVOIR_SIM(
      beg.y = 2020,
      beg.m = 1,
      Q        = clim$Q_sim,
      K        = 120,
      K_dead   = 20,
      T_dom    = demand_mon$value * dem_reduce,
      T_irr    = irr_demand * dem_reduce,
      T_eco    = demand_ces_mcm$env,
      evap.m   = (clim$PE - clim$Prec) * 1e-3,
      f_elev   = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      C_alloc  = 0.8)

    #Calculate reliabilities
    nofail <- length(which(sim$R_dom >= 1 * sim$T_dom))
    dom_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= 1 * sim$T_irr))
    irr_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= 1 * sim$T_eco))
    eco_rel  <- 100 * nofail / nrow(sim)

    out[s,'Qmean'] <- mean(clim$Q_sim)
    out[s,'dom_rel'] <- dom_rel
    out[s,'eco_rel'] <- eco_rel
    out[s,'irr_rel'] <- irr_rel
    out[s,'nclim'] <- c
    out[s,'reduction'] <- dem_reduce
  }

  filter(out, nclim == 150)


  df <- OUT %>% as.data.frame() %>% as_data_frame() %>%
    select(Id = V1, rel_dom = V2, rel_irr = V3, rel_eco = V4,
           crel_dom = V5, crel_irr = V6, crel_eco = V7) %>%
    left_join(results_tbl, ., by = "Id") %>%
    select(Id, dataset:alpha, rel_dom:crel_eco)


  #EXPERIMENTAL DESIGN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  span_dataset <- "Princeton"
  span_K <- 120
  span_Clim = nClim_select
  span_demL = 2
  span_irrig <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
  span_demL  <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

  results_tbl <- expand_grid(dataset = span_dataset, K = span_K,
      nClim = span_Clim, IrrMult = span_irrig, DemMult = span_demL) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, K, nClim, IrrMult, DemMult)

  #Loop through climate futures
  end <- nrow(results_tbl)

  ptm <- proc.time()
  cl  <- makeCluster(6)
  registerDoParallel(cl)

  OUT <- foreach(s=1:end, .combine = rbind, .packages = c("lubridate", "sirad"))      %dopar% {

    #Set parameters based on current iteration
    dataset  <- results_tbl$dataset[s]
    par_calib <- hydro_pars[,dataset]

    c <- results_tbl$nClim[s]
    K <- results_tbl$K[s]
    demM <- results_tbl$DemMult[s]
    irrM <- results_tbl$IrrMult[s]
    dom_demand <- filter(demand_levels, level == 2) %>%
      mutate(value = value * demM) %$% value

    #climate forcings
    clim <- forcing %>% filter((nClim == c) & (Dataset == dataset))
    Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
    clim$Tdif <- clim$Tmax - clim$Tmin
    clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
    clim$Q_sim <- ABCD_QEST(
      par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

    #### CALCULATE IRRIGATION DEMAND ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #Calculate ET0 (Penmann - Monteith equation (FAO))
    irrig_cal <- clim %>%
      left_join(crop_pars_et0, by = "Month") %>%
      mutate(
        EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
        Wind  = Wind * 0.01157407,
        v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
        vap_p = Humid/100 * v_sat,
        ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
          tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
          days  = Date, lat = Latitude) * as.numeric(days_in_month(Month))) %>%
      select(Year, Month, Prec, Tavg, Tmin, Tdif, PE, EPrec, ET0)

    #Calculate CRW (mm/month)
    irr_demand <- irrig_cal %>%
      left_join(crop_Kc, by = "Month") %>%
      left_join(crop_area, by = c("Season","Crop")) %>%
      mutate(
        Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
      mutate(
        ETc = ET0 * Kc,
        NIR_mm = pmax(ETc - EPrec, 0),
        GIR_mm = NIR_mm / (Eff/100),
        GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
      mutate(Month = as.integer(Month)) %>%
      group_by(Year, Month) %>%
      summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm

    irr_demand <- irr_demand * irrM

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
      T_dom  = dom_demand,
      T_irr  = irr_demand,
      T_eco  = demand_ces_mcm$env,
      evap.m = (clim$PE - clim$Prec) * 1e-3,
      f_elev = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      cycle = FALSE)

    #Calculate reliabilities
    nofail <- length(which(sim$R_dom >= 1 * sim$T_dom))
    dom_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= 1 * sim$T_irr))
    irr_rel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= 1 * sim$T_eco))
    eco_rel  <- 100 * nofail / nrow(sim)

    #Calculate critical reliabilities
    nofail <- length(which(sim$R_dom >= critical_level * sim$T_dom))
    dom_crel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_irr >= critical_level * sim$T_irr))
    irr_crel  <- 100 * nofail / nrow(sim)
    nofail <- length(which(sim$R_eco >= critical_level * sim$T_eco))
    eco_crel  <- 100 * nofail / nrow(sim)

    c(s, dom_rel, irr_rel, eco_rel, dom_crel, irr_crel, eco_crel)
  }
  stopCluster(cl)

  results_tbl %<>% mutate(Id = 1:n())

  df <- OUT %>% as.data.frame() %>%
    as_data_frame() %>%
    select(Id = V1, Rel_dom = V2, Rel_irr = V3, Rel_eco = V4,
           cRel_dom = V5, cRel_irr = V6, cRel_eco = V7) %>%
    left_join(results_tbl, ., by = "Id") %>%
    select(Id, dataset:DemMult, Rel_dom:cRel_eco)
  #write_csv(x = df, path = "./inputs/adaptation_dreduction.csv")


  # RESULTS ANALYSIS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  reduction_dat <- read_csv("./inputs/adaptation_dreduction.csv")

  data <- reduction_dat %>%
    select(nClim:cRel_eco) %>%
    group_by(IrrMult, DemMult) %>%
    summarize_each(funs(median)) %>%
    select(-nClim, -cRel_eco, -Rel_eco) %>%
    ungroup() %>%
    mutate(Id = 1:n(),
           IrrMult  = factor(IrrMult, levels = span_irrig),
           DemMult  = factor(DemMult, levels = span_demL),
           Dom_acc = ifelse(Rel_dom >= 99, 1, 0),
           Irr_acc = ifelse(Rel_irr >= 75, 1, 0),
           Combo_acc = ifelse((Dom_acc == 1) & (Irr_acc == 1), 1, 0))

  scaling <- seq(0.3,1,0.1) * 100

  p <- ggplot(data, aes(x=DemMult, y = IrrMult)) +
    theme_bw(base_size = 8) +
    scale_x_discrete(expand=c(0,0), labels = scaling) +
    scale_y_discrete(expand=c(0,0), labels = scaling) +
    scale_fill_manual(values = c("darkgreen", "red"), limits = c('1','0'),
      labels = c("Acceptable", "Not acceptable"), name = "Condition") +
    labs(x = "Domestic reduction (%)", y = "Irrigation reduction (%)") +
    guides(fill = FALSE)

  p1 <- p + geom_tile(aes(fill = as.factor(Dom_acc)), color = "white")
  p2 <- p + geom_tile(aes(fill = as.factor(Irr_acc)), color = "white")
  p3 <- p + geom_tile(aes(fill = as.factor(Combo_acc)), color = "white")

  library(cowplot)
  px <- plot_grid(p1, p2, labels = c("a)","b)","c)"), nrow = 2)
  ggsave("robust_demand_dom.png", p1, height = 3, width = 3)
  ggsave("robust_demand_irr.png", p2, height = 3, width = 3)


    filter((demL == 1)) %>%
    filter(nClim %in% nClim_q)

  p <- ggplot(data, aes(x = as.factor(Buffer), color = as.factor(nClim))) +
    theme_bw(base_size = 11) +
    scale_y_continuous(limits = c(40,100), breaks = seq(40, 100, 20)) +
    scale_color_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3',
      '#ff7f00'), name="Condition", breaks=c("525", "583", "637", "316", "385"),
      labels=c("Very wet", "Wet", "Median", "Dry", "Very Dry"))


# 4) DROUGHT ANALYSIS ----------------------------------------------------------

# SIMULATION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  #Packages required
  library(zoo)

  #Baseline data
  forcing_dat <- forcing %>% filter(dataset == "Princeton")

  forcing_yr <- forcing_dat %>%
      group_by(nclim, Year) %>%
      summarize(Prec = sum(Prec))

  clim_select_index <- clim_tbl %>% filter(temp == 3) %$% nclim

  forcing_rm <- forcing_yr %>%
    filter(nclim %in% clim_select_index) %>%
    group_by(nclim) %>%
    mutate(value = rollmean(x=Prec, 3, align="right", fill = NA)) %>%
    ungroup() %>%
    arrange(value) %>%
    mutate(value = round(value,2))

  df <- forcing_rm %>% na.omit() %>%
    mutate(dummy = "3-yr rolling mean")

  p <- ggplot(df, aes(x = dummy, y = value)) +
    theme_bw() +
    geom_violin() +
    labs(x = "", y = "MCM") +
    geom_jitter(width = 0.2, alpha = 0.5, size = 0.4, color = "gray") +
    geom_point(aes(y = 668), color = "red", size = 2); p

  #ggsave(plot = p, filename = "rollmean.tiff", width = 5, height = 5)

  #Find 10th quantile
  q10 <- quantile(forcing_rm$value, probs = 0.1, na.rm = TRUE, type = 3) %>%
    as.numeric()
  forcing_q10 <- filter(forcing_rm, value == q10)

  forcing_q10_ts <- forcing_dat %>%
    filter(nclim == forcing_q10$nclim) %>%
    filter(Year %in% {forcing_q10$Year + c(-4,-3,-2,-1,0)})

  #Find quantile range: 8 to 12th
  forcing_x <- forcing_rm %>%
    mutate(rank = percent_rank(value)) %>%
    filter(rank > 0.08 & rank < 0.12) %>%
    mutate(clim_index = 1:n())

  #Prepare grid for the simulations
  span_dataset <- "Princeton"
  span_K <- c(80,100,120,130,140)
  span_nclim <- forcing_x$clim_index
  span_demand <- 80
  span_buffer <- 1

  results_tbl <- expand_grid(
    dataset = span_dataset, K = span_K, clim_index = span_nclim) %>%
    mutate(dataset = as.character(dataset)) %>%
    arrange(dataset, K, clim_index) %>%
    left_join(forcing_x, by = "clim_index") %>%
    mutate(Id = {1:n() %>% as.character()})

  #Loop through climate futures
  end <- nrow(results_tbl)
  sim <- vector(mode = "list", length = end)

  for (s in 1:end) {

  #prepare climate data
  clim_num <- filter(forcing_x, clim_index==results_tbl$clim_index[s]) %$% nclim
  clim_year <- filter(forcing_x, clim_index==results_tbl$clim_index[s]) %$% Year

  clim_i <- forcing_dat %>%
    filter(nclim == clim_num) %>%
    filter(Year %in% {clim_year + c(-4,-3,-2,-1,0)})

  #Calculate hydrology
  clim <- clim_i
  Date <- as_date(paste(clim$Year, clim$Month, "1", sep="-"))
  clim$Tdif <- clim$Tmax - clim$Tmin
  clim$PE = HARGREAVES(Date, clim$Tavg, clim$Tdif, lat = 3.5)
  clim$Q_sim <- ABCD_QEST(
    par_calib, clim$Prec, clim$PE, S_ini=100, G_ini=2) * area * 1e-3

  #### CALCULATE IRRIGATION DEMAND ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Calculate ET0 (Penmann - Monteith equation (FAO))
  irrig_cal <- clim %>%
    left_join(crop_pars_et0, by = "Month") %>%
    mutate(
      EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
      Wind  = Wind * 0.01157407,
      v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
      vap_p = Humid/100 * v_sat,
      ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
        tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
        days = Date, lat = Latitude) * as.numeric(days_in_month(Month))) %>%
    select(Year, Month, Prec, Tavg, Tmin, Tdif, PE, EPrec, ET0)

  #Calculate CRW (mm/month)
  irr_demand <- irrig_cal %>%
    left_join(crop_Kc, by = "Month") %>%
    left_join(crop_area, by = c("Season","Crop")) %>%
    mutate(
      Area_ha = ifelse((Year > 2033) & Crop == "Mango", 0, Area_ha)) %>%
    mutate(
      ETc = ET0 * Kc,
      NIR_mm = pmax(ETc - EPrec, 0),
      GIR_mm = NIR_mm / (Eff/100),
      GIR_mcm = GIR_mm * Area_ha / 10^5) %>%
    mutate(Month = as.integer(Month)) %>%
    group_by(Year, Month) %>%
    summarize(Irr_mcm = sum(GIR_mcm)) %$% Irr_mcm

    #Set parameters based on current iteration
    K <- results_tbl$K[s]
    demand <- 80
    buffer <- 1

    # Simulate reservoir operation
    sim[[s]] <- reservoir_sim(
      beg.y = 2031,
      beg.m = 1,
      Q  = clim$Q_sim,
      K  = K,
      buffer = buffer,
      pool_buff = 0.5,
      pool_cons = 0.5,
      K_dead = 20,
      T_dom = demand,
      T_irr = 0,
      T_eco = 0,
      evap.m = (clim$PE - clim$Prec) * 1e-3,
      f_elev = f_elev, f_vol = f_vol, f_sarea = f_sarea,
      C_alloc = 0.8)
  }

  sim_df <- sim %>% bind_rows(.id = 'Id') %>%
    left_join(select(results_tbl, -Year), by = "Id")

  write_csv(x = sim_df, path = "./data/adaptation_drought.csv")

# RESULTS PLOTTING +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  date_lims <- c(as.Date("2031-01-01"), as.Date("2036-01-01"))
  date_breaks <- as.Date("2031-01-01") + years(0:5)

  df_data <- sim_df %>%
  filter(clim_index == 1)
  #df_data <- read_csv("./inputs/adaptation_drought.csv")

  df_temp <- df_data %>%
    mutate(date = as.Date(paste(Year, Month, "01",sep = "-")))

  df <- df_temp %>%
    mutate(Domestic = T_dom - R_dom) %>%
    select(date, Capacity = K, Storage = S, Domestic) %>%
    gather(key = variable, value = value, Storage:Domestic)
    #select(date, BufferC = buffer, Storage = S, Domestic) %>%
    #gather(key = variable, value = value, -BufferC, -date)


  p1 <- ggplot(filter(df, variable == "Storage"), aes(x = date, y = value,
      group = as.factor(Capacity), color = as.factor(Capacity))) +
    theme_bw(base_size = 10) +
    geom_line() +
    scale_x_date(limits = date_lims, breaks = date_breaks, date_labels = "%Y") +
    geom_hline(yintercept = 20, linetype = "dashed") +
    scale_y_continuous(limits = c(20,140), breaks = seq(20,140,30)) +
    labs(x = "", y = "Storage (MCM)", color = "Capacity\n(MCM)")

  p2 <- ggplot(filter(df, variable == "Domestic"), aes(x = date, y = value,
    group = as.factor(Capacity), color = as.factor(Capacity))) +
    theme_bw(base_size = 10) +
    geom_line() +
    scale_x_date(limits = date_lims, breaks = date_breaks, date_labels = "%Y") +
    labs(x = "", y = "Deficit (MCM)", color = "Capacity\n(MCM)")

  library(cowplot)
  grobs <- ggplotGrob(p1 + theme(legend.position="right"))$grobs
  legend_b <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

  p_i <- plot_grid(p1 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   align = "h", ncol = 1, labels = c("a", "b"))

  p <- plot_grid(p_i, legend_b, rel_widths = c(1,.2)); p
  ggsave(filename = "dry_year_bufferC.tiff", plot = p, width = 7, height = 5)

  #DRY-YEAR QUANTILES

  #Average deficit
  #Duration of failure

  #Year, Month, Inflow, Demand, Release,  Spill

  df <- sim_df %>%
    mutate(Id = as.integer(Id))

  out <- data_frame(Id = 1:max(df$Id), Resilience = NA, Deficit = NA)

  for (i in 1:max(df$Id)) {

    data <- filter(df, Id == i)
    data <- data[-(1:24),]

    sf_mon <- data$T_dom - data$R_dom

    if(length(which(sf_mon > 0)) == 0) {
      Res <- 1
    } else {
      fail_indices  <- which(sf_mon > 0)
      fail_recovery <- length(which(sf_mon[fail_indices + 1] == 0))
      Res <- fail_recovery / length(fail_indices)
    }

    mean_deficit <- ifelse(length(which(sf_mon > 0)) == 0, 0,
      sum(sf_mon/data$T_dom) / length(which(sf_mon > 0)))

    out[i,"Resilience"] <- Res
    out[i,"Deficit"] <- mean_deficit
  }

  dat <- out %>%
    left_join(mutate(results_tbl, Id = as.numeric(Id)), by = "Id") %>%
    group_by(K) %>%
    summarize(Resilience = mean(Resilience), Deficit = mean(Deficit))



  df <- df_data %>%
    select(Year, Month, Inflow = pool, Demand = T_dom, Release = R_dom, Spill = Spl)

  sf_mon <- df$Demand - df$Release



  df_data <- sim_df %>%
  filter(clim_index == 1)
  #df_data <- read_csv("./inputs/adaptation_drought.csv")




  #df_data <- read_csv("./inputs/adaptation_drought.csv")






