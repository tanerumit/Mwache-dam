

################################################################################

  source("Mw_Input_pars.R")

# ------------------------------------------------------------------------------

# CLIMATE DATA CLEANING --------------------------------------------------------


#OBSERVED CLIMATE DATA (Monthly, from 10 rain gauge stations)

  obs_prec_long <- read.csv("./Inputs/Obs_monthly_precip_mm.csv") %>%
    gather(key = Month, value = value, -Year, -Station_ID, -Station_Name) %>%
    mutate(Month = as.vector(sapply(as.character(Month),
      function(x) which(x == month.abb)))) %>%
    transform(Date = as.Date(paste(Year, Month, 1, sep="-"))) %>%
    select(Station_ID, Station_Name, Date, value) %>%
    arrange(Station_ID, Date) %>%
    as_data_frame()

#GRIDDED CLIMATE DATA
  GPCC_3.25 <- read.table("./Repository/Climate/GPCC/GPCC_-3.25_39.25") %>%
    as_data_frame() %>%
    select(Year = V1, Month = V2, Value = V3) %>%
    filter(Year %in% 1950:1999) %>%
    mutate(Source = "GPCC_3.25")

  GPCC_3.75 <- read.table("./Repository/Climate/GPCC/GPCC_-3.25_39.75") %>%
    as_data_frame() %>%
    select(Year = V1, Month = V2, Value = V3) %>%
    filter(Year %in% 1950:1999) %>%
    mutate(Source = "GPCC_3.75")

  CRU_3.25 <- read.table("./Repository/Climate/CRU/CRU_-3.25_39.25") %>%
    as_data_frame() %>%
    select(Year = V1, Month = V2, Value = V3) %>%
    filter(Year %in% 1961:1999) %>%
    mutate(Source = "CRU_3.25")

  CRU_3.75 <- read.table("./Repository/Climate/CRU/CRU_-3.25_39.75") %>%
    as_data_frame() %>%
    select(Year = V1, Month = V2, Value = V3) %>%
    filter(Year %in% 1961:1999) %>%
    mutate(Source = "CRU_3.75")

  grid_prec_long <- bind_rows(GPCC_3.25, GPCC_3.75, CRU_3.25, CRU_3.75) %>%
    spread(key = Source, value = Value) %>%
    mutate(Princeton = hist_forcings$Prec, Date = hist_dates) %>%
    transmute(Date, CRU = CRU_3.75, GPCC = GPCC_3.75, Princeton) %>%
    gather(key = Source, value = value, -Date)


  #write.csv(x = grid_prec_long,
  #  file = "./inputs/grid_prec_long_mm.csv", row.names=FALSE)
  #write.csv(x = obs_prec_long,
  #file = "./inputs/obs_prec_long_mm.csv", row.names=FALSE)

# ------------------------------------------------------------------------------

# STAKEHOLDER DATA -------------------------------------------------------------

  #Packages to load
  require(readxl)
  require(dplyr)
  require(tidyr)
  require(ggplot2)

  p <- list()

  belief_df <- read_excel("inputs/Workshop_survey_v2.xlsx", sheet = 1) %>%
    gather(key = Particip, value = value, x1:x28) %>%
    mutate(value = as.factor(.$value))

  for (i in 1:12) {

    plot_name <- unique(belief_df$Uncertainty)[i]
    plot_data <- filter(belief_df, Uncertainty == plot_name)

  p[[i]] <- ggplot(plot_data, aes(x = value, fill = Factor)) +
    ggtitle(plot_name) +
    geom_bar() +
    facet_grid(.~ Factor) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_discrete(expand=c(0,0)) +
    labs(x="Significance score", y = "count (out of 28)") +
    theme(plot.title = element_text(face="bold"),
          strip.text = element_text(size = 6)) +
    guides(fill = FALSE)

    ggsave(paste0("factor_",i,".png"), p[[i]], width = 8, height = 4, units = "in")

  }

# ------------------------------------------------------------------------------

# MET_DATA_ANALYSIS ------------------------------------------------------------

  library(corrplot)
  library(scales)

  #Read-in data
  grid_prec_long <- read_csv("./inputs/grid_monthly_prec_mm.csv")
  obs_prec_long  <- read_csv("./inputs/ground_monthly_prec_mm.csv") %>%
    select(Date, Source = Station_ID, value) %>%
    mutate(Source = paste0("G_", Source))

  #Read-in monthly streamflow time-series (cms) & convert to mm
  q_obs_df <- read_csv("inputs/flow_gage_cms.csv", skip=1) %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep="-")),
           Flow = Flow_obs * 86.4 * days_in_month(Date) / area)

  #Data in wide-format
  grid_prec <- grid_prec_long %>% spread(key = Source, value = value)
  obs_prec  <- obs_prec_long %>% spread(key = Source, value = value)

  #Correlation matrix of met-stations
  cor_matrix <- cor(obs_prec[,-1], use="pairwise.complete.obs")
  corrplot(corr = cor_matrix, type = "upper",   addCoef.col="black",
    cl.lim = c(0,1), cl.length = 5, tl.cex = 1, tl.pos = "lt",
    method = "shade", tl.col="black", tl.srt = 45, diag = TRUE)

  #Scatter plot of data by pairs
  #ids <- expand.grid(x = colnames(obs_prec)[-1], y = colnames(obs_prec)[-1])
  #df <-  obs_prec_wide
  #for (i in 1:100) {
  #
  #  p <- GG_SCATTER(df[as.character(ids[i,1])], df[as.character(ids[i,2])],
  #    save.plot = TRUE, height = 10, width = 10)
  #}

  #Diagnostic plots (observed precip)
  data <- gather(grid_prec, key = Source, value = value, -Date)

  p <- viz_timeSeries(data = data, x = "Date", y = "value", type = "Source"); p
  p <- viz_flowDuration(data = data, x = "Date", y = "value", type = "Source"); p
  p <- viz_boxplot(data = data, x = "Date", y = "value", type = "Source"); p




  plots <- pl_TimeSeries(date = grid_prec$Date, data = grid_prec[,-1],
    unit = "mm", agg_type = "sum", titles = FALSE, Guides = TRUE)
  #ggsave("timeseries_grid_prec.png", plots[[1]], height = 10, width = 20, unit = "cm")

  #Ground-data versus gridded-data
  df <- obs_prec %>% left_join(grid_prec, by = "Date") %>%
    left_join(q_obs_df, by ="Date")

  #Ground-data versus gridded-data
  viz_scatter(df, x ="G_9339048", y= "GPCC", equal.scale = FALSE)

  GG_SCATTER(x=df["G_9339048"], y=df["GPCC"], save.plot = TRUE, height = 10, width = 10)
  GG_SCATTER(x=df["G_9339048"], y=df["CRU"], save.plot = TRUE, height = 10, width = 10)
  GG_SCATTER(x=df["G_9339048"], y=df["Princeton"], save.plot = TRUE, height = 10, width = 10)

  #Gridded-data versus gridded-data
  GG_SCATTER(x=df["Princeton"], y=df["GPCC"], save.plot = T, height = 10, width = 10)
  GG_SCATTER(x=df["Princeton"], y=df["CRU"], save.plot = T, height = 10, width = 10)
  GG_SCATTER(x=df["GPCC"], y=df["CRU"], save.plot = T, height = 10, width = 10)

  #Gridded-data vs observed flow
  GG_SCATTER(x=df["Flow_obs"], y=df["GPCC"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["Princeton"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["CRU"], save.plot = T, height = 10, width = 10, equal.scale = F)

  #Ground data vs observed flow
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339000"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339002"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339006"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339023"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339048"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339049"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339051"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339055"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339057"], save.plot = T, height = 10, width = 10, equal.scale = F)
  GG_SCATTER(x=df["Flow_obs"], y=df["G_9339078"], save.plot = T, height = 10, width = 10, equal.scale = F)

# ------------------------------------------------------------------------------

# WAVELET ANALYSIS  ------------------------------------------------------------

  #Annual precipitation series, 50 years
  selection <- c(4553, 1189, 533, 1450, 311, 3880, 215, 893, 1659, 4350, 3555)
  warm_prec_annual <- data.frame(read.csv(
    paste0(getwd(),'/inputs/WARM_SIM_50_5000.csv'), skip=1, header=TRUE))
  warm_prec_annual <- t(warm_prec_annual[selection,])
  colnames(warm_prec_annual) <- seq_len(ncol(warm_prec_annual))
  stoc_labels <- paste0('Real_',1:ncol(warm_prec_annual))
  P_annual_stoc <- matrix(warm_prec_annual,
    c(nrow(warm_prec_annual), ncol(warm_prec_annual)))

  #save.image('Mb_Climate_Data.Rdata')

  #Annual precipitation time-series
  ANNUAL_PRCP <- forcings_gcm_annual %>% dplyr::filter(ID == 1) %>% .[["Prec"]]

  #Box-cox transformation of the original data
  lambda <- powerTransform(ANNUAL_PRCP)$lambda
  TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)

  #Check normality of original ts
  data <- ANNUAL_PRCP
  model <- data.frame(x=1:length(data), y=data)
  pl_1 <- ggplot(model, aes(x,y)) + geom_line() + labs(y='mm/year', x='') +
    ggtitle('time-series') + theme_bw()

  pl_2 <- ggplot(model, aes(x=y)) +
    geom_histogram(fill='gray90', color='blue') + theme_bw() +
    labs(x='mm/year', y='count') + ggtitle('histogram')
  pl_3 <- ggplot(model, aes(sample = y)) + stat_qq()  +  theme_bw() +
    geom_abline(intercept = mean(model$y), slope = sd(model$y)) +
    ggtitle('qqplot')

  p <- arrangeGrob(pl_1, arrangeGrob(pl_2, pl_3, ncol=2), ncol=1)

  #ACF ON HISTORIC TIME-SERIES
  lambda <- powerTransform(ANNUAL_PRCP)$lambda
  TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)
  pl_historic_ACF_org  <- PLOT_ACF(
    ANNUAL_PRCP,lag=10, name.arg='Historical')[['ACF']]
  pl_historic_ACF_trns <- PLOT_ACF(
    TRANSFORM_DATA,lag=10, name.arg='Historical')[['ACF']]

  print(p)

  #WAVELET ON HISTORIC PRECIPITATION +++++++++++++++++++++++++++++++++++++++++++
  CLIMATE_VARIABLE <- ANNUAL_PRCP
  source(paste0(scripts_dir, "/WAVELET_ANALYSIS.R"))
  gws_obs <- GWS
  gws_sig <- signif_GWS
  gws_period <- period

  var_title <- c('GCMs', 'Stochastic')
  var_list <- list(P_annual_gcms, P_annual_stoc)
  gws_list <- list(GCMs=NA, Stochastic=NA)

  pl_WARM_GWS  <- list(GCMs=NA, Stochastic=NA)
  pl_WARM_IND  <- list(GCMs=NA, Stochastic=NA)
  pl_WARM_COR  <- list(GCMs=NA, Stochastic=NA)

  #WAVELET ON STOCHASTIC REALIZATIONS & HISTORICAL GCM RUNS ++++++++++++++++++++
  for (r in 1:length(var_title)) {

    var_annual <- var_list[[r]]
    gws_sim <- array(NA, c(length(period), ncol(var_annual)))
    if (r==1) {
      colnames(gws_sim) <- gcms_data$Name[1:20]
    } else {
      colnames(gws_sim) <- paste0('Stochastic_', seq_len(ncol(var_annual)))
    }

    for (x in seq_len(ncol(var_annual))) {

      CLIMATE_VARIABLE <- as.numeric(unlist(var_annual[,x]))
      source(paste0(scripts_dir, "/WAVELET_ANALYSIS.R"))
      gws_sim[,x] <- GWS

    }

    gws_list[[r]] <- gws_sim
    gws_mean <- apply(gws_sim, 1, mean)
    gws_min <- apply(gws_sim, 1, function(x) quantile(x, 0.025))
    gws_max <- apply(gws_sim, 1, function(x) quantile(x, 0.975))
    gws_df <- data.frame(gws_period, gws_sig, gws_obs, gws_mean, gws_min, gws_max)

    #GLOBAL WAVELET SPECTRA (ALL Time-series)
    p <- ggplot(gws_df, aes(x=gws_period)) +
      theme_bw(base_size=12) +
      ggtitle(var_title[r]) +
      geom_ribbon(aes(ymin=gws_min, ymax=gws_max), alpha=0.15) +
      geom_line(aes(y=gws_mean), color='black', size=0.6) +
      geom_line(aes(y=gws_sig), color="red", linetype="dashed", size=0.6) +
      geom_line(aes(y=gws_obs), color="blue", size=0.6) +
      labs(x="Period (years)", y=expression(paste("Power (", mm^2,")"))) +
      scale_x_continuous(breaks=seq(0,50,5), expand=c(0,0)) +
      scale_y_log10(limits = c(5*10**2,5*10**5),labels = comma)

    pl_WARM_GWS[[r]]  <- p

    #GLOBAL WAVELET SPECTRA (Plot each time-series seperately)
    gws_df2 <- data.frame(Period=gws_period,gws_sim/10^3) %>%
       gather(Type, Simulated, -Period)

    p <- ggplot(gws_df2, aes(x=Period, y=Simulated)) +
      theme_bw(base_size=font_size) +
      geom_line(color='black') +
      geom_line(aes(x=gws_period, y=gws_obs/10^3), gws_df, color='blue') +
      geom_line(aes(x=gws_period, y=gws_sig/10^3), gws_df, linetype='dashed', color='red') +
      facet_wrap( ~ Type, scales='free') +
      labs(x = expression(paste("Historical Power (", mm^2,")")),
           y = expression(paste("Simulated Power (", mm^2,")"))) +
      scale_y_continuous(expand=c(0,0), label=comma) +
      scale_x_continuous(expand=c(0,0), limits=c(0,40), breaks=seq(0,40,10),
        label=comma) +
      theme(strip.text.x = element_text(face="bold"),
            strip.text.y = element_text(face="bold"),
            strip.background = element_rect(fill="white"))

    pl_WARM_IND[[r]] <- p

    #CORRELATION ANALYSIS
    models <- lapply(seq_len(ncol(gws_sim)), function(k)
      lm(y ~ x, data.frame(x=gws_obs, y=gws_sim[,k])))
    r_sqrs <- round(sapply(models, function(k)
      as.numeric(format(summary(k)$r.squared, digits = 3))),3)
    r_sqrs <- paste0(r_sqrs, collapse=' ') #paste0('r^2=', r_sqrs)

    gws_cor <- data.frame(x=gws_obs/10^3,gws_sim/10^3) %>%
      gather(Type, y, -x)

    p <- ggplot(gws_cor, aes(x,y)) +
      theme_bw(base_size=font_size) +
      geom_point(shape=1, size=0.8) +
      stat_smooth(method = "lm", se=TRUE, color="black", formula = y ~ x) +
      facet_wrap( ~ Type, scales='free') +
      #geom_text(aes(x, y, label = r_sqrs), vjust=2, size=font_size-4,
      #  data=data.frame(x=75, y=Inf, label=r_sqrs, Type=unique(gws_cor$Type)),) +
      labs(x = expression(paste("Historical Power (", mm^2,")")),
           y = expression(paste("Simulated Power (", mm^2,")"))) +
      scale_y_continuous(expand=c(0,0), label=comma) +
      scale_x_continuous(expand=c(0,0), limits=c(0,100),
                         breaks=seq(0,100,25), label=comma) +
      theme(plot.title = element_text(face="bold"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.text.x = element_text(face="bold"),
            strip.text.y = element_text(face="bold"),
            strip.background = element_rect(fill="white"))


    pl_WARM_COR[[r]] <- p

  }

  plots_GWS <- arrangeGrob(pl_WARM_GWS[[1]], pl_WARM_GWS[[2]], ncol=2)


  #MEAN, STANDARD DEVIATION, KURTOSIS, SKEW
  GCMs_stats <- P_annual_gcms %>%
    summarise_each(funs(mean, sd, skewness, kurtosis)) %>%
    gather(variable, value) %>%
    separate(variable, c("var", "stat"), sep = "\\_") %>%
    spread(stat, value) %>%
    mutate(Type = 'GCMs') %>%
    select(Type, var, mean, sd, kurtosis, skewness)


  stoc_stats <- data.frame(P_annual_stoc) %>%
    summarise_each(funs(mean, sd, skewness, kurtosis)) %>%
    gather(variable, value) %>%
    separate(variable, c("var", "stat"), sep = "\\_") %>%
    spread(stat, value) %>%
    mutate(Type = 'Stochastic') %>%
    select(Type, var, mean, sd, kurtosis, skewness)

  hist_stats <- data.frame(ANNUAL_PRCP) %>%
    summarise_each(funs(mean, sd, skewness, kurtosis)) %>%
    mutate(Type='Historic', var = 1) %>%
    select(Type, var, mean, sd, kurtosis, skewness)

  P_annual_stats <- rbind(GCMs_stats, stoc_stats, hist_stats)

  p1 <- ggplot(P_annual_stats, aes(color=Type, x=mean, y=sd)) +
    geom_point(size=4) + theme_bw() + guides(color=FALSE) +
    labs(x='mean', y='standard deviation')

  p2 <- ggplot(P_annual_stats, aes(color=Type, x=skewness, y=kurtosis)) +
    geom_point(size=4) + theme_bw() +
    labs(x='skewness', y='kurtosis')

  pl_stat <- arrangeGrob(p1, p2, ncol=2, widths=c(1,1.3))

 # 3) AUTOCORRELATION ANALYSIS ++++++++++++++++++++++++++++++++++++++++++++++++
  #ACF OF HISTORICAL TIME-SERIES
  pl_ACF <- PLOT_ACF(ANNUAL_PRCP,lag=10, name.arg='Historical')
  pl_ACF[['ACF']]

  #ON CMIP5 HISTORICAL RUNS
  gcms_acf <- lapply(seq_len(ncol(P_annual_gcms)),
    function(x) PLOT_ACF(P_annual_gcms[,x],lag=10, name.arg=gcms_labels[x]))
  gcms_acf <- lapply(gcms_acf, function(x) x[['ACF']])
  pl_gcms_acf <- do.call("arrangeGrob", c(gcms_acf, ncol=4))

  #ON STOCHASTIC REALIZATIONS
  stoc_acf <- lapply(seq_len(ncol(P_annual_stoc)),
    function(x) PLOT_ACF(P_annual_stoc[,x],lag=10, name.arg=stoc_labels[x]))
  stoc_acf <- lapply(stoc_acf, function(x) x[['ACF']])
  pl_stoc_acf <- do.call("arrangeGrob", c(stoc_acf, ncol=4))

# ------------------------------------------------------------------------------

# SEDIMENT MANAGEMENT ----------------------------------------------------------

  #Trapping efficiency
  #Baseline-conditions:
  Q_df <- read_csv("./Inputs/historic_flow_mcm_princeton.csv") %>%
    mutate(Year = year(Date), Month = month(Date)) %>%
    select(-Date) %>% group_by(Year) %>%
    summarize(Q = sum(Q))

  #Trapping efficiency (Vorosmorthy et al., 2003)
  TE_dam <- 1 - (0.05 * (designs$Volume/mean(Q_df$Q))^0.5)

  # Land use classes (sediment rate given m3/km2-year)
  land_use <- data_frame(
    Type  = c("Low_Forest", "Closed_Forest", "Woodland", "Closed_Grassland",
      "Open_Grassland", "Urban", "Cropland"),
    Area  = c(45, 434.25, 402.75, 1055.25, 2.25, 112.50, 225),
    Erosion  = c("Low", "Low", "Low", "Med", "Med", "High", "High"),
    Sed_rate = c(500, 500, 500, 1000, 1000, 1500, 1500))

  #Annual sedimentation load
  sed_flux <- land_use %>%
    mutate(Load = (Area * Sed_rate) / 10^6) %>%
    summarize(Load = sum(Load) * TE_dam[2])

# ------------------------------------------------------------------------------


#
# # MASS-BALANCE ANALYSIS --------------------------------------------------------
#
#   #Mass-balance accounting for the projected demand & supplies
#   year_set <- 2015:2065
#   file_path <- "./inputs/BWSS_balance.xlsx"
#
#   demand_dat <- read_excel(path = file_path, sheet = "Demand", skip = 1)
#
#   demand_tbl <- data_frame(Year = year_set) %>%
#     mutate(Pop = approxExtrap(
#       x=demand_dat$Year, demand_dat$Population, xout = year_set)$y) %>%
#     mutate(Demand =  approxExtrap(
#       x=demand_dat$Year, demand_dat$Demand_mcm, xout = year_set)$y)
#
#   supply_tbl <- read_excel(path = file_path, sheet = "Supply", skip = 1) %>%
#     mutate_each(funs(. * (365/1e6)), -Year) %>%
#     select(Year, Tiwi, Marere, Others, Mzima, Baricho, Msambweni, Mwache) %>%
#     gather(key = Source, value = Yield, -Year)
#
#   var_cols <- c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd",
#                 "#969696","#737373","#253494")
#
#   #Plot cummulative yield from alternative sources
#   p3 <- ggplot(supply_tbl, aes (x = Year)) +
#     theme_bw() +
#     geom_area(mapping=aes(y = Yield, fill = Source)) +
#     scale_fill_brewer(palette = "Blues") +
#     scale_x_continuous(expand = c(0,0),
#                        limits = c(2015, 2065), breaks = seq(2015,2065,10)) +
#     scale_y_continuous(expand = c(0,0),
#                        limits = c(0,350), breaks = seq(0, 350, 50)) +
#     geom_line(data = demand_tbl, mapping = aes(y=Demand),
#               linetype = "longdash", size = 1) +
#     labs( x = "Year", y = "Cummulative yield (MCM)")
#
#   #Surplus storage from Mwache
#   balance_tbl <- supply_tbl %>%
#     spread(key = Source, value = Yield) %>%
#     mutate(Other_supply = Tiwi + Marere + Mzima + Baricho + Others + Msambweni,
#            Demand       = demand_tbl$Demand,
#            Net          = Other_supply + Mwache - Demand,
#            Balance      = as.factor(ifelse(Net > 0, 1, 0))) %>%
#     select(Year, Demand, Other_supply, Mwache, Net, Balance)
#
#   net.rle = rle(balance_tbl$Net < 0)
#   balance_tbl$group = rep.int(1:length(net.rle$lengths), times=net.rle$lengths)
#
#   p4 <- ggplot(
#     data = balance_tbl,
#     mapping = aes(x = Year, y = Net, group = group, fill = Balance)) +
#     geom_area() +
#     scale_fill_manual(values  = c('red', 'blue'), labels = c("Def.", "Surp")) +
#     scale_y_continuous(limits = c(-150,50), breaks = seq(-150,50,50)) +
#     labs(x = "Year", y = "Surplus supply (MCM)") +
#     guides(fill = FALSE)
#
#
#
