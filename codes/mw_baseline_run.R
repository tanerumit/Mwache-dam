
#------------------------------------------------------------------------------#
#   Title   :    Decision Tree Project - MWACHE HYDROLOGY/SYSTEMS MODEL
#   By      :    Mehmet Umit Taner
#   Date    :    September 25, 2015
#   Version :    1.0
#   Notes   :    baseline analysis (HISTORICAL CONDITIONS)
#
#------------------------------------------------------------------------------#

  source("Mw_Input_pars.R")

  #Monthly precip forcings & calibration parameters
  prec_grid <- read_csv("./inputs/grid_monthly_prec_mm.csv")
  parm_grid <- read_csv("./Inputs/calib_pars_11_2015.csv")
  name_grid <- names(parm_grid)

# ******************************************************************************


# HYDROLOGY SIMULATION ---------------------------------------------------------

  #Simulation period
  date_sim <- as.Date("1950-01-01") + c(0:599) * months(1)

  #List object to store hydroclimate table
  hydroclim_sim <- as.list(rep(NA,3))

  #Loop through the grid data-sets
  for (i in 1:length(name_grid)) {

    #Select met. dataset to use
    met_data <- name_grid[i]

    #Calibration paramters
    parm_sim <- parm_grid[[met_data]]

    #Append precip to the Princeton dataset
    hydroclim <- filter(forcings_gcm_monthly, ID == 1) %>%
      mutate(
        Date = as.Date(paste(Year, Month, "1", sep="-")),
        Tdif = Tmax - Tmin,
        Prec = (filter(prec_grid, Source == met_data) %>% .[["value"]]),
        PE = HARGREAVES(Date, Tavg, Tdif, lat = BasinLat)) %>%
      filter(!is.na(Prec)) %>%
      mutate(key = met_data)

    #Simulated streamflow
    hydroclim_sim[[i]] <- hydroclim %>%
      filter(Date %in% date_sim) %>%
      mutate(Q_mm =  ABCD_QEST(parm = parm_sim, P = Prec, PE = PE,
        S_ini = 1, G_ini = 1), Q_mcm = Q_mm * area * 1e-3) %>%
      select(key, Date, Prec:Tmax, Tdif, PE, Q_mcm)

    #p <- ggplot(hydroclim_sim, aes(x=Date, y = Prec)) + geom_line(); p
    #df <- data_frame(Date = date_per, Q = Q_sim_mcm)
    #write_csv(df, paste0("./Inputs/", "hydro_clim_", met_data,".csv"))

  }

  #df <- bind_rows(hydroclim_sim) %>%
  #  select(key, Date, Q = Q_mcm) %>%
  #  spread(key = key, value = Q) %>%
  #  write_csv("Historical_inflows_cms.csv")


# YIELD-RELIABILITY ANALYSIS ---------------------------------------------------

  #Yield & environmental demand levels
  yieldT  <- c(seq(50, 130, 10))
  envir.m <- demand_ces_mcm$env

  # Create the input-output matrix([# d x # yield_levels ~ met_data]
  out_tbl <- expand.grid(Design = designs$Volume, Yield = yieldT) %>%
    as_data_frame() %>%
    mutate(iter = 1:n(), Rel = NA, df = rep(list(NA),  n())) %>%
    select(iter, Design, Yield, Rel, df)

  out_lst <- as.list(rep(NA,3))

  #Progress bar & timing
  pb <- txtProgressBar(min = 1, max = nrow(out_tbl), style = 3)

  # Explore yield-reliability relationship
  for (s in 1:length(name_grid)) {

    #Select met. dataset & parameters
    met_data <- name_grid[s]
    Q <- hydroclim_sim[[s]]$Q_mcm

    #Loop through each scenario
    for (i in 1:nrow(out_tbl)) {

      setTxtProgressBar(pb, i)

      #Set current values for iteration
      K <- out_tbl$Design[[i]]
      yield.m <- out_tbl$Yield[[i]] %>%
        multiply_by(as.numeric(days_in_month(1:12)/365))

      out_tbl$df[[i]] <- RESERVOIR_SIM(beg.y  = 1950, beg.m = 1,
        K = K, Q = Q, K_d = 20, K_cons = K, K_buff = K *.9,
        ratingc = rating_curve, evap.m = netevap,
        Tar_dom  = yield.m, Tar_irr  = rep(0,12), Tar_eco  = envir.m)

      #Calculate reliablity
      Rel <- with(out_tbl$df[[i]], R_eco + R_irr + R_dom)
      Tar <- with(out_tbl$df[[i]], Tar_eco + Tar_irr + Tar_dom)
      out_tbl$Rel[[i]] <- (1 - (sum((Rel/Tar) < 1)/length(Tar))) * 100

    }

  close(pb)

  #Save results to the list object
  out_lst[[s]] <- out_tbl

  }

  #Plotting **********************************

  df <- bind_rows(out_lst, .id = "Key") %>%
    mutate(Design = as.factor(Design),
           Source = as.factor(name_grid[as.numeric(Key)]))

  #Plot yield-reliability curve for all designs
  p <- ggplot(df, aes(x = Yield, y = Rel, group = Design, color = Design)) +
    theme_bw(base_size = 10) +
    geom_line(size = 0.8) +
    facet_wrap( ~ Source, labeller = label_both, ncol = 1) +
    geom_vline(xintercept=80.3, show.legend = F, linetype="dashed", size=1) +
    scale_y_continuous(limits = c(60,100), breaks = seq(60,100,10)) +
    scale_x_continuous(limits = c(50,120), breaks = seq(50,120,10)) +
    labs(x = "Yield (MCM)",  y= "Reliability (%)") +
    scale_color_brewer(name="Design\nCapacity\n(MCM)", palette="Spectral")

  # Compare yields for the de, in which the likelihood estimates of the
  # plausible future outcomes and their impacts are incomplete, sign options
  Trel <- c(85, 90, 95, 99)
  df2 <- expand.grid(Trel=Trel, Design = as.factor(designs$Volume),
    Source = name_grid) %>% as_data_frame() %>%
    mutate(Yield = NA)

  for (i in 1:nrow(df2)) {
    df3 <- filter(df, (Source == df2$Source[i]) & (Design == df2$Design[i]))
    df2$Yield[i] <- approx(x = df3$Rel, y = df3$Yield, xout = df2$Trel[i])$y
  }

  p <- ggplot(df2, aes(x= Yield, y = Trel)) +
    theme_bw() +
    geom_point(aes(color = Source, size = Design)) +
    scale_y_continuous(limits = c(80,100), breaks=seq(80,100,5)) +
    scale_color_brewer(name="Data\nSource", palette="Dark2") +
    labs(x = "Yield (MCM)", y = "Reliability (%)"); p

 results <- filter(df2, (Trel == 95) & (Design == 120))


# ******************************************************************************


# ******************************************************************************


# ALLOCATION SIMULATION --------------------------------------------------------

  #Select the met-dataset to analyze
  met_data <- "Princeton"
  hydroclim_i <- hydroclim_sim[[which(met_data == name_grid)]]

  #DOMESTIC-USE PARAMETERS
  par <- list()
  par$dom_baseline <- 80.3  # Target domestic yield (MCM/year)
  par$dom_scale  <- 1.0     # Scaling factor for yield
  par$dom_recycle  <- 0.1   # fraction of recycled water use

  #IRRIGATION PARAMETERS
  par$irr_baseline <- 5425  # planned irrigation area (ha)
  par$irr_scale    <- 0.5   # Scale coef. for irrigation area
  par$irr_spr_eff  <- 0.65  # efficiency for the sprinkler system
  par$irr_dip_eff  <- 0.90  # efficiency for the dip system

  #ENVIRONMENTAL RELEASE PARAMETERS
  par$eco_baseline <- demand_ces_mcm$env #monthly release requirement (MCM/mo)
  par$eco_scale <- 1.0      # Scale environmental release targets

  #Reservoir parameters
  par$res_K <- 120          # Reservoir storage capacity (MCM)

  #Allocation scenarios
  vec_irr_area <- c(1000,2000,3000,4000,5000,5425, 6000)
  vec_irr_scale <- vec_irr_area/ par$irr_baseline

  #Define the output tbl to store inputs & outputs
  out <- data_frame(Iter = 1:length(vec_irr_scale))
  out %<>% mutate(par = rep(list(NA), n()), sim = par, result = par)

################################################################################

  #Loop through all iterations
  pb <- txtProgressBar(min = 1, max = nrow(out), style = 3)
  for (i in 1:nrow(out)) {

    setTxtProgressBar(pb, i)

    par$irr_scale  <- vec_irr_scale[i]
    source("Mw_prep_input.R")

    sim <- RESERVOIR_SIM(beg.y  = input$Year[1], beg.m = input$Month[1],
      K_d = 20, ratingc = rating_curve, evap.m = netevap,
      K = par$res_K, K_cons = par$res_K, K_buff = par$res_K*0.9,
      Q = input$Q_mcm, Tar_dom  = input$Dom_mcm,
      Tar_irr  = input$Irr_mcm, Tar_eco  = input$Eco_mcm, cycle = FALSE,
      priority = c("Tar_eco", "Tar_dom", "Tar_irr"))

    #RELIABILITY CALCULATIONS
    result <- list()
    result$Rel_eco <- with(sim,
      (1 -(sum(ifelse(Tar_eco == 0, 1, R_eco/Tar_eco) < 1)/ length(Tar_eco))))
    result$Rel_dom <-  with(sim,
      (1 -(sum(ifelse(Tar_dom == 0, 1, R_dom/Tar_dom) < 1)/ length(Tar_dom))))
    result$Rel_irr <- with(sim,
      (1 -(sum(ifelse(Tar_irr == 0, 1, R_irr/Tar_irr) < 1)/ length(Tar_irr))))

    out[[i,"par"]] <- par
    out[[i,"sim"]] <- sim
    out[[i,"result"]] <- result

  }
  close(pb)

################################################################################

# ******************************************************************************

################################################################################


# RESULTS VISUALIZATION --------------------------------------------------------


  #devtools::install_github("wilkelab/cowplot")
  library(cowplot)

  # VISUALIZE RELIABILITIES
  Rel_df <- data_frame(Iter = 1:nrow(out))
  Rel_df$Ecological <- sapply(1:nrow(out), function(x) out$result[x][[1]]$Rel_eco)
  Rel_df$Domestic   <- sapply(1:nrow(out), function(x) out$result[x][[1]]$Rel_dom)
  Rel_df$Irrigation <- sapply(1:nrow(out), function(x) out$result[x][[1]]$Rel_irr)

  df <- gather(Rel_df, key = Demand, value = value, -Iter) %>%
    mutate(Scenario = as.factor(Iter), Demand = as.factor(Demand), value = value * 100)

  p2 <- ggplot(df, aes(x = Scenario, y = value, color = Scenario)) +
    theme_set(theme_gray(base_size = 10)) +
    facet_wrap(~ Demand, ncol = 1, labeller = label_both) +
    scale_y_continuous(limits = c(70,100), breaks = seq(70,100,10)) +
    scale_color_discrete(name = "Irrig (ha)",
                         labels = vec_irr_scale * par$irr_baseline)  +
    geom_point(size = 3.0) +
    labs(x = "Scenario", y = "Reliability"); p2

  Alloc_eco <- sapply(1:nrow(out), function(x) out$sim[x][[1]]$Tar_eco) %>%
    as.data.frame() %>% as_data_frame() %>%
    mutate(Date = date_sim, Demand = "Ecological") %>%
    gather(key = Iter, value = value, -Date, -Demand)

  Alloc_dom <- sapply(1:nrow(out), function(x) out$sim[x][[1]]$Tar_dom) %>%
    as.data.frame() %>% as_data_frame() %>%
    mutate(Date = date_sim, Demand = "Domestic") %>%
    gather(key = Iter, value = value, -Date, -Demand)

  Alloc_irr <- sapply(1:nrow(out), function(x) out$sim[x][[1]]$Tar_irr) %>%
    as.data.frame() %>% as_data_frame() %>%
    mutate(Date = date_sim, Demand = "Irrigation") %>%
    gather(key = Iter, value = value, -Date, -Demand)

  Alloc_df <- bind_rows(Alloc_eco, Alloc_dom, Alloc_irr) %>%
    mutate(Year = year(Date), Month = month(Date)) %>%
    group_by(Year, Demand, Iter) %>%
    summarize(value = sum(value)) %>%
    ungroup()

  Res_df <- data_frame(Iter = 1:nrow(out))
  Res_df$Storage <- sapply(1:nrow(out), function(x) out$sim[x][[1]]$S)


  p1 <- ggplot(Alloc_df, aes(x = Year, y = value, color = Iter, group = Iter)) +
    theme_set(theme_gray(base_size = 10)) +
    facet_wrap(~ Demand, ncol = 1, labeller = label_both) +
    scale_y_continuous(limits = c(0,75), breaks = c(0,25,50,75)) +
    scale_color_discrete(name = "Allocation \nScenario", labels = 1:11) +
    guides(color = FALSE) +
    geom_line(size = 1) + labs(x="", y ="MCM"); p1

  plot_grid(p1, p2, nrow = 1, align = "h", labels = c("A","B"))

  ##############################################################################

  Storage_df <- sapply(1:nrow(out), function(x) out$sim[x][[1]]$S) %>%
    as.data.frame() %>% as_data_frame() %>%
    mutate(Date = date_sim) %>%
    gather(key = Iter, value = value, -Date)

  p3 <- ggplot(Storage_df, aes(x = Date, y = value, color = Iter, group = Iter)) +
    theme_set(theme_gray(base_size = 10)) +
    geom_line(size = 1, alpha = 0.5) + labs(x="", y ="MCM") +
    scale_color_discrete(name = "Allocation \nScenario", labels = 1:11);p3

# ******************************************************************************


# RESULTS VISUALIZATION (SINGLE RUN) -------------------------------------------

  #Results file
  df_monthly <- out %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep ="-")))

  res_zones <- select(df_monthly, Date) %>%
    mutate(Capacity = K, Buffer = K * 0.9, Inactive = K_d) %>%
    gather(key = Levels, value = S, -Date)

  #Reservoir storage
  p1 <- ggplot(mapping = aes(x = Date, y = S)) +
    scale_y_continuous(limits =c(0,K+10), breaks = seq(0,K+10, length.out =5)) +
    geom_line(data = df_monthly) +
    labs(x = "Date", y = "Storage (MCM)"); p1

  #With operational levels
  p1 <- p1 + geom_line(data = res_zones, aes(color=Levels)); p1

  df_annual <- df_monthly %>%
    group_by(Year) %>%
    summarize(S=mean(S), Q=sum(Q), R_eco = sum(R_eco), Evap = sum(L),
              Spill = sum(Spl),R_dom = sum(R_dom), R_irr = sum(R_irr)) %>%
    gather(key = Type, value = value, -Year)

  p2 <- ggplot(data = df_annual, aes(x=Year, y=value, color = Type))  +
    geom_line(); p2

  #Deficit time-series
  df_deficits <- df_monthly %>%
    mutate(Def_eco = Tar_eco - R_eco,
           Def_dom = Tar_dom - R_dom,
           Def_irr = Tar_irr - R_irr) %>%
    select(Date, Def_eco : Def_irr) %>%
    gather(key = Deficits, value = value, -Date) %>%
    separate(Date, c("Year","Month","Day")) %>%
    mutate(Year = as.numeric(Year))

  p_deficit <- ggplot(df_deficits, aes(x = Year)) +
    geom_bar(aes(y = value, fill = Deficits), stat = "identity") +
    labs(x = "Year", y = "MCM") +
    scale_y_continuous(limits=c(0,80))

################################################################################

# ******************************************************************************

