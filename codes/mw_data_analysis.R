
# Required input data
source("./Codes/Mw_input.R")
require(readxl)
library(corrplot)
library(Kendall)

# STAKEHOLDER SURVEY DATA ------------------------------------------------------

  ### 2015 SURVEY

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

    ggsave(paste0("factor_",i,".png"), p[[i]], width=8, height=4, units="in")

  }

  ### 2016 SURVEY

  results <- read_excel("data/workshop_2016.xlsx",
    sheet = "Survey1") %>%
    gather(key = Particip, value = value, x1:x23) %>%
    mutate(value = as.factor(.$value))

  p <- list()
  for (i in 1:10) {

    plot_name <- unique(results$Uncertainty)[i]
    plot_data <- filter(results, Uncertainty == plot_name) %>%
      filter(!is.na(value))

    df <- plot_data %>% group_by(Factor, value) %>%
      summarize(Freq = n()) %>%
      group_by(Factor) %>%
      mutate(TFreq = (Freq / sum(Freq)) * 100)

    p[[i]] <- ggplot(df, aes(x = value, y = TFreq, fill = Factor)) +
      theme_bw() +
      ggtitle(paste0("Q",i,". ",plot_name)) +
      geom_bar(stat = "identity") +
      facet_wrap(~ Factor, nrow = 1) +
      scale_x_discrete(expand=c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
      labs(x="rating score", y = "(%)") +
      theme(plot.title = element_text(face="bold"),
            strip.text = element_text(size = 9)) +
      guides(fill = FALSE)

    ggsave(paste0("./Outputs/factor_",i,".png"),
      p[[i]], width = 7, height = 3, units = "in")

  }

  p_1 <- plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], align = "v", ncol = 1)
  ggsave(plot = p_1, file="all_figures.png")

  survey_dat <- read_excel("data/workshop_2016.xlsx", sheet = 1)
  df <- survey_dat %>%
    filter(Allocation != "NA" & Uses != "Other", Ranking != "NA") %>%
    mutate(Uses = factor(Uses,
      levels = c("Domestic", "Environmental", "Irrigation", "Hydropower")),
           Group = as.factor(Group),
           Allocation = as.numeric(Allocation),
           Ranking = as.numeric(Ranking))

  p1 <- ggplot(df, aes(x = Ranking, fill = Uses)) +
    facet_wrap( ~ Uses) +
    geom_bar(stat = "count") +
    guides(fill = FALSE) +
    scale_y_continuous(limits = c(0,8), breaks = seq(0,8,2)) +
    labs(x = "Preference ranking", y = "count")
  ggsave(filename = "survey_2016_1.png", plot = p1, height = 4, width = 7)

  p2 <- ggplot(df, aes(x = Uses, y = Allocation, fill = Uses)) +
    geom_boxplot() +
    guides(fill = FALSE) +
    labs(x = "Beneficial uses", y = "fraction of supply allocated (%)")
  ggsave(filename = "survey_2016_2.png", plot = p2, height = 4, width = 7)


# MASS-BALANCE ANALYSIS ------------------------------------------------------

  #Mass-balance accounting for the projected demand & supplies
  year_set <- 2015:2069
  file_path <- "./inputs/water_balance.xlsx"

  demand_dat <- data_frame(
    Year  = seq(2015, 2065, 5),
    value = c(84, 104, 137, 161, 187, 212, 237, 262, 287, 312, 338))

  demand_tbl <- demand_dat %>%
    do(Demand = as.data.frame(
      Hmisc::approxExtrap(x=.$Year, .$value, xout=year_set))) %>%
    unnest(Demand) %>%
    select(Year = x, value = y)

  supply_tbl <- read_excel(path = file_path, sheet = "Supply", skip = 1) %>%
    mutate_each(funs(. * (365/1e6)), -Year) %>%
    select(Year, Tiwi, Marere, Others, Mzima, Baricho, Msambweni, Mwache) %>%
    gather(key = Source, value = Yield, -Year)

  #Plot cummulative yield from alternative sources
  p3 <- ggplot(supply_tbl, aes (x = Year)) +
    theme_bw(base_size = 10) +
    geom_area(mapping=aes(y = Yield, fill = Source)) +
    scale_fill_brewer(palette = "Set3") +
    scale_x_continuous(limits = c(2020, 2070), breaks = seq(2020,2070,10)) +
    scale_y_continuous(limits = c(0,400), breaks = seq(0, 400, 100)) +
    geom_line(data = demand_tbl, mapping = aes(y=value),
              linetype = "longdash", size = 1) +
    labs( x = "Year", y = "MCM")

  #Surplus storage from Mwache
  balance_tbl <- supply_tbl %>%
    spread(key = Source, value = Yield) %>%
    mutate(Other_supply = Tiwi + Marere + Mzima + Baricho + Others + Msambweni,
           Demand       = demand_tbl$value,
           Net          = Other_supply + Mwache - Demand,
           Balance      = as.factor(ifelse(Net > 0, 1, 0))) %>%
    select(Year, Demand, Other_supply, Mwache, Net, Balance)

  net.rle = rle(balance_tbl$Net < 0)
  balance_tbl$group = rep.int(1:length(net.rle$lengths), times=net.rle$lengths)

  p4 <- ggplot(
    data = balance_tbl,
    mapping = aes(x = Year, y = Net, group = group, fill = Balance)) +
    theme_bw(base_size = 10) +
    geom_area() +
    scale_fill_manual(values  = c('red', 'blue'), labels = c("Def.", "Surp")) +
    scale_x_continuous(limits = c(2020, 2070), breaks = seq(2020,2070,10)) +
    scale_y_continuous(limits = c(-150,50), breaks = seq(-150,50,50)) +
    labs(x = "Year", y = "Surplus supply (MCM)") +
    guides(fill = FALSE)

  library(cowplot)
  p <- plot_grid(p3, p4, align = "v", ncol = 1)
  p



# IRRIGATION DATA  -------------------------------------------------------------

  require(sirad)

  #Crop efficiency and Kc values for each crop type
  crop_data <- read_csv("./inputs/crop_data_irrig.csv", skip = 1)

  #Monthly parameters for Penmann - Montheith calculation
  crop_pars_et0 <- read_csv("./inputs/crop_pars_et0.csv", skip = 1)

  #Crop coefficients
  crop_Kc <- crop_data %>%select(Season, Crop = c, Jun:May) %>%
    gather(key = Month, value = Kc, -Season, -Crop) %>%
    rowwise() %>%
    mutate(Month = which(Month == month.abb)) %>%
    mutate(Month = as.numeric(Month))

 #Crop area & efficiency
  crop_area <- crop_data %>%
    select(Season, Crop = c, Area_ha, Method) %>%
    mutate(Eff=ifelse(Method == "sprinkler", par$irr_spr_eff, par$irr_dip_eff)) %>%
    select(-Method)

  #Calculate ET0 (Penmann - Monteith equation (FAO))
  forcing <- hydroclim_i %>%
    mutate(Year = as.numeric(year(Date)), Month = as.numeric(month(Date))) %>%
    left_join(crop_pars_et0, by = "Month") %>%
    mutate(
      EPrec = ifelse(Prec > 75, Prec * 0.8 -25, pmax(0, Prec * 0.6 - 10)),
      Wind  = Wind * 0.01157407,
      v_sat = 0.6108 * 2.7183^(17.27*Tavg/(Tavg + 273.3)),
      vap_p = Humid/100 * v_sat,
      ET0 = et0(Tmax=Tmax, Tmin=Tmin, vap_pres=vap_p, sol_rad = Solar,
        tal = cst(Solar, Date, radians(Latitude)),z = 150, uz = Wind,
        days  = Date, lat = Latitude) * as.numeric(days_in_month(Month)))

  #Calculate CRW (mm/month)
  irr_per_crop <- forcing %>%
    left_join(crop_Kc, by = "Month") %>%
    left_join(crop_area, by = c("Season","Crop")) %>%
    mutate(
      ETc = ET0 * Kc,
      NIR_mm = pmax(ETc - EPrec, 0),
      GIR_mm = NIR_mm / Eff,
      GIR_mcm = GIR_mm * Area_ha * par$irr_scale / 10^5)

  # Monthly total irrigation demand
  demand.m <- irr_per_crop %>%
    group_by(Year, Month) %>%
    summarize(Irr_mcm = sum(GIR_mcm)) %>%
    ungroup() %>%
    mutate(Year = as.numeric(Year), Month = as.numeric(Month))

  #DOMESTIC
  dom_y <- par$dom_baseline * par$dom_scale * (1 - par$dom_recycle)
  demand.m %<>% mutate(Dom_mcm = dom_y * coeff.m[Month])

  #ENVIRONMENTAL
  demand.m %<>% mutate(Eco_mcm = par$eco_baseline[Month] * par$eco_scale)

  input <- forcing %>%
    left_join(demand.m, by=c("Year","Month")) %>%
    select(key, Year, Month, Prec:Tmax, PE, Q_mcm, ET0, Irr_mcm, Dom_mcm, Eco_mcm)

# ENVIRONMENTAL FLOWS ----------------------------------------------------------

  # With historical flows, in mcm/month
  historical_flow <- read_csv("inputs/flow_gage_cms.csv", skip=1) %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep="-")),
           value = Flow_int * 86400 * days_in_month(Date) / 10^6,
           Year = as.factor(Year), Month = as.factor(Month)) %>%
    select(Year, Month, value) %>%
    filter(value != "NA")

  # With simulated streamflow, in mcm/month
  forcing <- read_csv("./Inputs/forcings_.csv", skip = 1)
  clim <- filter(forcing, Index == 1) %>%
    mutate(Date = as.Date(paste(Year, Month, "1", sep="-")),
           Tdif = Tmax - Tmin,
           PE = HARGREAVES(Date, Tavg, Tdif, lat = 3.5))

  #Simulate streamflow (mm) and (MCM/mo)
  simulated_flow <- clim %>%
    mutate(Q_sim_mm = ABCD_QEST(
      par_calib, P = clim$Prec, PE, S_ini = 100, G_ini = 2),
      Q_mcm = Q_sim_mm * area * 1e-3) %>%
    select(Year, Month, value = Q_mcm) %>%
    mutate(Year=as.factor(Year), Month=as.factor(Month), value=as.numeric(value))

  #df <- historical_flow
  df <- simulated_flow

  df1 <- df %>%
    group_by(Month) %>%
    summarize(lower = quantile(value, 0.05), upper = quantile(value, 0.95)) %>%
    mutate(Type ="0.05-0.95")
  df2 <- df %>%
    group_by(Month) %>%
    summarize(lower = quantile(value, 0.01), upper = quantile(value, 0.99)) %>%
    mutate(Type ="0.01-0.99")
  df3 <- df %>%
    group_by(Month) %>%
    summarize(lower = quantile(value, 0.25), upper = quantile(value, 0.75)) %>%
    mutate(Type ="0.25-0.75")

  df_median <- df %>% group_by(Month) %>% summarize(value = median(value))
  df_ranges <- bind_rows(df1, df2, df3)

  ggplot(mapping = aes(Month)) +
    theme_bw() +
    geom_ribbon(aes(ymin = lower, ymax = upper, color = Type, group = Type),
                data = df_ranges, alpha = 0.1) +
    geom_line(aes(y = value, group = 1), df_median, color = "blue", size = 1) +
    scale_x_discrete(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    labs(fill = "Quantile range", color = "Quantile range") +
    scale_color_manual(values = c("0.25-0.75" = "steelblue",
                                  "0.05-0.95" = "green", "0.01-0.99" = "orange"))

  env_flows <- data_frame(Month = 1:12,
    Policy1 = df1$lower * 0.3, Policy2 = df1$lower * 0.5, Policy3 = df1$lower,
    CES = demand_ces_mcm$env)

  p <- env_flows %>%
    gather(key = variable, value = value, -Month) %>%
    {ggplot(., aes(x=Month, y=value, color =variable)) +
        geom_line()};p


# HYDROCLIMATE DATA ANALYSIS ---------------------------------------------------

  #READ-IN OBSERVED CLIMATE DATA
  precip_obs_grid   <- read_csv("./data/precip_grid_mm.csv")
  precip_obs_ground <- read_csv("./data/precip_ground_mm.csv") %>%
    select(Date, Source = ID, value) %>%
    mutate(Source = paste0("G_", Source))

  #READ-IN STREAMFLOW DATA
  flow_data <- read_csv("data/flow_gage_cms.csv", skip=1) %>%
    mutate(Date = as.Date(paste(Year, Month,"1", sep="-")),
           Flow = Flow_obs * 86.4 * lubridate::days_in_month(Date) / area)

  # p <- ggplot(flow_data, aes(x = Date)) +
  #   theme_bw() +
  #   geom_line(aes(y=Flow_obs), color = "black") +
  #   labs(x = "", y = "cms"); p
  # p <- p + geom_line(aes(y=Flow_int), color = "blue"); p


  # CORRELATION MATRIX OF PRECIP +++++++++++++++++++++++++++++++++++++++++++++++
  df <- precip_obs_ground %>% spread(key = Source, value = value)
  #data <- precip_obs_grid %>% spread(key = Source, value = value)

  cormatrix <- df %>% select(-Date) %>% cor(., use="pairwise.complete.obs")
  corrplot(corr = cormatrix, type = "full", addCoef.col="black",
    cl.lim = c(0,1), cl.length = 5, tl.cex = 1, tl.pos = "lt",
    method = "shade", tl.col="black", tl.srt = 45, diag = TRUE)

  #data <- precip_obs_ground
  data <- precip_obs_grid

  # PRECIP - TIME-SERIES (MONTHLY)  ++++++++++++++++++++++++++++++++++++++++++++
  p <- ggplot(data, aes(x = Date, y = value)) +
    theme_bw(base_size = 11) +
    geom_line() +
    facet_wrap(~ Source, ncol = 1) +
    labs(x = "Date", y = "mm")
  #ggsave("precip_grid_tseries.png", p, width = 10, height = 8)

  p <- viz_tseries(data = data, x = "Date", y = "value", type = "Source") +
    labs(x = "Date", y = "mm")
  #ggsave("precip_grid_tseries.png", p, width = 10, height = 5)

  # PRECIP - TIME-SERIES (ANNUAL)  +++++++++++++++++++++++++++++++++++++++++++++
  precip_annual <- precip_obs_grid %>%
    mutate(Year = year(Date)) %>%
    group_by(Year, Source) %>%
    summarize(value = sum(value)) %>%
    filter(Source == "Princeton")

  p <- ggplot(precip_annual, aes(x = Year, y = value)) +
      theme_bw(base_size = 10) +
      geom_line() +
      scale_x_continuous(breaks = seq(1950,2000,5)) +
      stat_smooth(method = lm, se = FALSE) +
      labs(x = "Date", y = "mm")
  #ggsave("precip_Princeton_annual.png", p, width = 7, height = 4)

  #Mann-Kendal trend-test
  MannKendall(precip_annual$value)





  # PRECIP - MONTHLY - BOXPLOT  ++++++++++++++++++++++++++++++++++++++++++++++++
  p2 <- precip_obs_grid %>%
    mutate(Month = month(Date) %>% as.factor()) %>%
    filter(Source == "Princeton") %>%
    {ggplot(., aes(x = Month, y = value, group = Month)) +
        theme_bw(base_size = 10) +
        geom_boxplot() +
        stat_summary(fun.y=mean, colour="darkred", geom="point") +
        scale_x_discrete(labels = month.abb) +
        labs(x = "Date", y = "mm")}

  ggsave("precip_Princeton_seasonal.png", p2, width = 7, height = 4)


  # PRECIP - FLOW DURATION CURVE  ++++++++++++++++++++++++++++++++++++++++++++++
  p <- viz_FDC(data = data, x = "Date", y = "value", type = "Source") +
    labs(y = "mm")
  ggsave("precip_grid_fdc.png", p, width = 5, height = 5)

  # PRECIP - BOXPLOT  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  p <- viz_boxplot(data = data, x = "Date", y = "value", type = "Source") +
    labs(x = "Month", y = "Precip. (mm)");

  p <- p + scale_x_discrete(labels = month.abb)

  ggsave("precip_grid_boxplot.png", p, width = 7, height = 4)

  # PRECIP - SCATTER PLOTS +++++++++++++++++++++++++++++++++++++++++++++++++++++
  data <- precip_obs_grid %>% spread(key = Source, value = value) %>%
    slice(133:600) %>%
    left_join((precip_obs_ground %>% spread(key = Source, value = value)),
      by = "Date")

  #gridded-precip versus gridded-precip
  p1 <- viz_scatter(data, x ="Princeton", y= "GPCC", eq.scale = TRUE)
  p2 <- viz_scatter(data, x ="Princeton", y= "CRU", eq.scale =TRUE)
  p3 <- viz_scatter(data, x ="GPCC", y= "CRU", eq.scale = TRUE)

  p  <- plot_grid(p1, p2, p3, labels = c("a","b","c"))
  ggsave("precip_scatterplot_grid_vs_grid.png", p, width = 6, height = 6)

  #ground-precip versus gridded-precip
  p1 <- viz_scatter(data, x ="G_9339006", y= "Princeton", eq.scale = TRUE)
  p2 <- viz_scatter(data, x ="G_9339006", y= "GPCC", eq.scale =TRUE)
  p3 <- viz_scatter(data, x ="G_9339006", y= "CRU", eq.scale = TRUE)
  p  <- plot_grid(p1, p2, p3, labels = c("a","b","c"))
  ggsave("precip_scatterplot_ground_vs_grid.png", p, width = 6, height = 6)

  # FLOW VS PRECIP - SCATTER PLOTS +++++++++++++++++++++++++++++++++++++++++++++
  data <- precip_obs_grid %>% spread(key = Source, value = value) %>%
    left_join((precip_obs_ground %>% spread(key = Source, value = value)),
      by = "Date") %>% left_join(flow_obs, by = "Date")

  #Gridded-data vs observed flow
  p1 <- viz_scatter(data, x ="Flow_obs", y= "Princeton", eq.scale = F)
  p2 <- viz_scatter(data, x ="Flow_obs", y= "GPCC", eq.scale = F)
  p3 <- viz_scatter(data, x ="Flow_obs", y= "CRU", eq.scale = F)
  p  <- plot_grid(p1, p2, p3, labels = c("a","b","c"))
  ggsave("precip_scatterplot_flow_vs_grid.png", p, width = 6, height = 6)

  #Ground data vs observed flow
  p1 <- viz_scatter(data, x ="Flow_obs", y= "G_9339000", eq.scale = F)
  p2 <- viz_scatter(data, x ="Flow_obs", y= "G_9339002", eq.scale = F)
  p3 <- viz_scatter(data, x ="Flow_obs", y= "G_9339006", eq.scale = F)
  p4 <- viz_scatter(data, x ="Flow_obs", y= "G_9339023", eq.scale = F)
  p5 <- viz_scatter(data, x ="Flow_obs", y= "G_9339048", eq.scale = F)
  p6 <- viz_scatter(data, x ="Flow_obs", y= "G_9339049", eq.scale = F)
  p7 <- viz_scatter(data, x ="Flow_obs", y= "G_9339051", eq.scale = F)
  p8 <- viz_scatter(data, x ="Flow_obs", y= "G_9339055", eq.scale = F)
  p9 <- viz_scatter(data, x ="Flow_obs", y= "G_9339057", eq.scale = F)
  p10<- viz_scatter(data, x ="Flow_obs", y= "G_9339078", eq.scale = F)
  p  <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 4,
                  labels = c("a","b","c","d","e","f","g","h","i","j"))
  ggsave("precip_scatterplot_flow_vs_ground.png", p, width = 9, height = 12)


### WAVELET ANALYSIS ***********************************************************
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
  ANNUAL_PRCP <- forcings__annual %>% dplyr::filter(ID == 1) %>% .[["Prec"]]

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

  var_title <- c('s', 'Stochastic')
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


# CLIMATE PROJECTIONS  ---------------------------------------------------------

  #13 unique qualitative colors for each GCM group
  gcm_colors <- c('#a6cee3','#1f78b4','#b2df8a',
    '#33a02c','#fb9a99','#e31a1c',
    '#fdbf6f', '#ff7f00','#cab2d6','#6a3d9a',
    '#ffff99','#b15928', 'gray20',
    'steelblue',"royalblue")

  #GCM metadata
  gcm_gen <- read_csv("data/CMIP5_genealogy.csv", skip = 1) %>%
    select(Model = GCM, Family) %>%
    mutate(color = gcm_colors[Family])

  write_csv(gcm_proj, path = "./data/CMIP5_genealogy.csv")

  #BASED ON DATA FROM SUNGWOOK
  gcm_raw_data <- read_csv("./data/CMIP5_raw_deltaclim.csv", skip = 1)
  gcm_proj <- filter(gcm_raw_data, scenario != "hist") %>% na.omit() %>%
    select(Model = GCM, scenario, del_temp=temp, del_prec=prec) %>%
    mutate(del_prec = del_prec * 100) %>%
    mutate(scenario = factor(scenario,
      levels = c("rcp26", "rcp45", "rcp60", "rcp85"),
      labels = c("RCP 2.6", "RCP 4.5", "RCP 6.0", "RCP 8.5")))

  write_csv(gcm_proj, path = "CMIP5_deltaclim.csv")

  #PLOT CMIP5 PROJECTIONS (ALL)
  df <- gcm_proj #%>% filter(scenario %in% c("rcp26", "rcp45"))
  point_size <- 2.5

  p1 <- ggplot(df, aes(x = del_temp, y=del_prec, fill=Model, shape=scenario)) +
    theme_bw(base_size = 10) +
    labs(x = alabs[[1]], y = alabs[[2]]) +
    #Set geoms
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = point_size, stroke = 1) +
    #Set scales
    scale_shape_manual(name = "Scenario", values = c(21,22,23,24)) +
    scale_fill_manual(values = gcm_gen$color, labels = gcm_gen$Model) +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = c(-35,+55), breaks = atick$prec2)+
    #Set legend
    guides(
      fill = guide_legend(title.position = "top",
        override.aes = list(shape = 21,linetype = "dashed"),
            order = 1, ncol = 3, keyheight = 0.65),
      shape  = guide_legend(title.position = "top",
        order = 2, ncol = 1, keyheight = 0.70)) +
    theme(legend.direction = "horizontal",
      legend.position = "bottom",
      legend.box = "horizontal")


  ggsave(plot = p1, filename = "scatterpl_CMIP5.tiff", height = 7, width = 5)



  #Obtained from Brent Boehlert, Industrial Economics

  #145 Model projections: 600 months (50 years) x 4 variables in the order of:
  #1  Observed (based on the Princeton dataset)
  #23 CMIP5 hindcasts
  #56 CMIP3 GCM runs
  #43 CMIP5 GCM runs
  #22 CMIP5 GCM runs processed by Bruce Hewitson's group at UCT
  #(using pressure patterns to drive precip rather
  #than precip directly from the GCMs)

  #GCM metadata, forcing, and grouping information
  source("Codes/Mw_input.R")

  #13 unique qualitative colors for each GCM group
  gcm_colors <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
    '#fdbf6f', '#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 'gray20')

  gcm_meta <- read_csv("Inputs/gcm_metadata.csv", col_names = TRUE)

  gcm_meta %<>%
    left_join(gcm_groups, by = "Model") %>%
    mutate(Family = multipleReplace(Family,
      what = c(1,2,3,4,5,6,8,9,10,11,12,13,15), by= 1:13),
      Color  =  gcm_colors[Family])

  gcm_colorkey <- gcm_meta %>%
    filter(Ensemble == "CMIP5") %>%
    group_by(Model) %>% slice(1) %>%
    ungroup() %>%
    arrange(Family) %>%
    select(Model, Family, Color) %>%
    mutate(Model = as.factor(Model), Family = as.factor(Family))



  #GCM forcings
  forcings_gcm <- read_csv("./Inputs/forcing_gcm.csv", skip = 1)

  #Monthly and annual GCM data
  gcm_mon <- forcings_gcm
  gcm_yr  <- group_by(forcings_gcm, Year,Index) %>%
    summarize(prec=sum(Prec), temp=mean(Tavg)) %>% ungroup()

  # Climate statistics *** first define window period ***
  window_hist <- 1950:2000
  window_gcm  <- 2046:2050
  gcm_hist_window <- gcm_yr %>%
    filter(Year %in% window_hist) %>%
    left_join(gcm_meta %>% select(Index,Name,Ensemble,Scenario), by = "Index")
  gcm_proj_window <- gcm_yr %>%
    filter(Year %in% window_gcm) %>%
    left_join(gcm_meta %>% select(Index,Name,Ensemble,Scenario), by = "Index")

  #Absolute and relative means
  gcm_stats <- bind_rows(gcm_hist_window, gcm_proj_window) %>%
    group_by(Index) %>%
    summarize(mu_prec = mean(prec), mu_temp = mean(temp),
              sd_prec = sd(prec), sd_temp = sd(temp)) %>%
    mutate(del_prec = 100 + 100  * (mu_prec - mu_prec[1])/mu_prec[1],
           del_temp = mu_temp - mu_temp[1])

  #Append gcm meta information
  gcm_meta %<>% select(Index, Model, Ensemble, Scenario, Family, Color)
  gcm_stats %<>% left_join(gcm_meta, by = "Index")

  #Filter by model
  CMIP5_stats <- gcm_stats %>%
    filter(Ensemble == "CMIP5") %>%
    mutate(Scenario = factor(Scenario, levels = c("historical","rcp45","rcp85"),
                           labels = c("Historical","RCP 4.5","RCP 8.5")),
           Family = as.factor(Family)) %>%
    filter(Scenario != "Historical")

  #CMIP5_dat <- CMIP5_stats %>%
  #  mutate(del_prec = del_prec / 100)
  #write_csv(CMIP5_dat, "data_CMIP5.csv")

  df <- read_csv("GCM_full.csv") %>%
    filter(scenario != "hist") %>%
    na.omit() %>%
    select(Model = GCM, Scenario = scenario, del_temp = temp, del_prec = prec) %>%
    mutate(del_prec = del_prec * 100)

  range(df$del_temp)
  range(df$del_prec)

  #PLOT CMIP5 PROJECTIONS (ALL)
  df <- df
  point_size <- 2.5
  p1 <- ggplot(df, aes(x = del_temp, y=del_prec, fill=Model, shape=Scenario)) +
    theme_bw(base_size = 10) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 100, linetype = "dashed") +
    geom_point(size = point_size, stroke = 1) +
    #scale_shape_manual(name = "Scenario", values = c(21,22)) +
    scale_x_continuous(expand=c(0,0), limits = alim$temp, breaks = atick$temp) +
    scale_y_continuous(expand=c(0,0), limits = alim$prec, breaks = atick$prec,
      labels = seq(-30,50,10)) +
    #scale_fill_manual(values = gcm_colorkey$Color, labels = gcm_colorkey$Model) +
    guides(fill = guide_legend(
      override.aes = list(shape = 21,linetype = "dashed"),
      order = 1, ncol = 1, keyheight = 0.65),
      shape  = guide_legend(order = 2, ncol = 1, keyheight = 0.70)) +
    theme(legend.position = "bottom") +
    labs(x = alabs[[1]], y = alabs[[2]])

  #ggsave(plot = p1, filename = "scatterpl_CMIP5_2046_2050.png", height = 5, width = 6.5)

  #PLOT CMIP5 PROJECTIONS (GROUP MEANS)
  CMIP5_group_means <- CMIP5_stats %>%
    group_by(Family, Scenario) %>%
    summarize(del_prec = mean(del_prec), del_temp = mean(del_temp))

  df <- CMIP5_group_means

  p2 <- ggplot(df, aes(x = del_temp, y=del_prec)) +
    theme_bw(base_size = 10) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 100, linetype = "dashed") +
    geom_point(aes(fill=Family, shape=Scenario),
      data = df, size = point_size, stroke = 1) +
    scale_fill_manual(values = gcm_colorkey$Color, limits = gcm_colorkey$Family) +
    scale_shape_manual(name = "Scenario", values = c(21,22)) +
    scale_x_continuous(limits=range(0,4), breaks=aticks$x) +
    scale_y_continuous(limits=range(aticks$y1), labels = aticks$y2 , breaks=aticks$y1) +
    labs(x = alabs[[1]], y = alabs[[2]]) +
    guides(fill = FALSE)

  p2 <- p2 + stat_ellipse(level = 0.999, type = "norm")


  #ggsave(plot = p2, filename = "scatterpl_CMIP5m_2046_2050.png",
  #height = 5, width = 6.5)

  library(cowplot)
  p <- plot_grid(p1, p2, align = 'v', labels = c("a","b"), ncol = 1)
  ggsave('GCM_proj.png', p, height = 9, width = 6)


# DEMAND ANALYSIS ---------------------------------------------------------

  # # PROJECTIONS +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #
  # demand_proj <- data_frame(
  #   year = c(2015, 2035),
  #   "low"  = c(38,68),
  #   "medium" = c(38, 77),
  #   "high" = c(38,103)) %>%
  #   gather(key = scenario, value = value, -year) %>%
  #   mutate(scenario = factor(scenario, levels = c("low","medium","high")))
  #
  # demand_range <- data_frame(
  #   year = c(2015, 2035),
  #   min  = c(38,57), max = c(38, 114)) %>%
  #   gather(key = scenario, value = value, -year)
  #
  # demand_proj_yr <- demand_proj %>%
  #   group_by(scenario) %>%
  #   do(Demand = as.data.frame(approx(x=.$year, .$value, xout = 2015:2035))) %>%
  #   unnest(Demand) %>%
  #   select(scenario, year = x, value = y)
  #
  # demand_range_yr <- demand_range %>%
  #   group_by(scenario) %>%
  #   do(Demand = as.data.frame(approx(x=.$year, .$value, xout = 2015:2035))) %>%
  #   unnest(Demand) %>%
  #   select(scenario, year = x, value = y) %>%
  #   spread(key = scenario, value = value)
  #
  # df  <- demand_proj_yr
  # df2 <- demand_range_yr
  #
  # p <- ggplot(data = df, aes(x = year)) +
  #   theme_bw() +
  #   geom_ribbon(aes(min = min, max = max), data = df2, alpha = 0.3) +
  #   geom_line(aes(y = value, color = scenario), size = 0.8) +
  #   scale_y_continuous(limits = c(35,115), breaks = seq(30,120,10)) +
  #   labs(x = "Year", y = "Annual demand (MCM)")
  #
  # ggsave("demand_proj.png", p, height = 4, width = 7)
  #
  # #Develop demand projections
  # demand_proj <- demand_dat %>%
  #   group_by(Type, Scenario) %>%
  #   do(Demand = as.data.frame(
  #     approxExtrap(x=.$Year, .$value, xout=year_range))) %>%
  #   unnest(Demand) %>%
  #   select(Type, Scenario, Year = x, value = y) %>%
  #   mutate(value = value * 0.365) %>%
  #   group_by(Scenario, Type) %>%
  #   mutate(Dummy = ifelse(Year > 2035, 0, value)) %>%
  #   mutate(value2 = ifelse(Dummy == 0, max(Dummy), value)) %>%
  #   ungroup() %>%
  #   select(Type, Scenario, Year, value, value2)
  #
  # ggplot(data = demand_proj, aes(x = Year, color = Scenario)) +
  #   facet_wrap(~ Type, nrow = 2) +
  #   geom_line(aes(y = value), size = 0.8) +
  #   geom_line(aes(y = value2), linetype = "dashed", size = 0.8)
  #
  # mon_cf <- as.numeric(days_in_month(1:12)/365)
  #demand_proj %<>%
  #  expand.grid.df(., Mon = 1:12) %>%
  #  as_data_frame() %>%
  #  select(Type, Scenario, Year, Mon = y, value2) %>%
  #  mutate(value = value2 * mon_cf[Mon]) %>%
  #  arrange(Type, Scenario, Year, Mon)

  #write_csv(x = demand_proj, path = "./Inputs/Mw_DT_demand_proj.csv")



  #Projection range
  year_range <- 2015:2070
  demand_target <- 186000 * 365 * 10^-6
  demand_2015 <- 103000 * 365 * 10^-6

  # Scenarios by Consultant (cubic meters per day)
  demand_scn_cmd <- data_frame(
    Year = c(2015, 2018, 2024, 2030),
    Low  = c(103000,	109756, 133226, 161531),
    Med  = c(103000, 116166, 145892, 182519),
    High = c(103000, 128555, 173377, 233247))  %>%
    gather(key = scenario, value = value, -Year)

  demand_scn_mcm <- demand_scn_cmd %>% mutate(value = value * 365 * 10^-6)
  demand_scn <- demand_scn_mcm

  #Demand projections
  demand_scn_yr <- demand_scn %>%
    group_by(scenario) %>%
    do(Demand = as.data.frame(
      Hmisc::approxExtrap(x=.$Year, y = .$value, xout =  year_range))) %>%
    unnest(Demand) %>%
    select(scenario, Year = x, value = y) %>%
    group_by(scenario) %>%
    mutate(value2 = ifelse(Year < 2034, value, value[which(Year == 2034)])) %>%
    ungroup() %>%
    mutate(scenario = factor(scenario, levels = c("Low", "Med","High"))) %>%
    arrange(scenario)

  filter(demand_scn_yr, Year == 2034)

  #Delta change by scenario & grid
  demand_change <- demand_scn_yr %>%
    filter(Year == 2034) %>%
    mutate(ValScn = value,
           DelScn = (ValScn - demand_2015)/demand_2015 * 100,
           DelGrd = c(60, 120, 180),
           ValGrd = DelGrd/100 * demand_2015 + demand_2015)

  #Demand scenarios for stress test (mcm/yr)
  demand_grid <- data_frame(
    Year = c(2015, 2034),
    dem1  = c(demand_2015, 60),
    dem2  = c(demand_2015, 83),
    dem3  = c(demand_2015, 105)) %>%
    gather(key = scenario, value = value, -Year)

  demand_grid_yr <- demand_grid %>%
    group_by(scenario) %>%
    do(Demand = as.data.frame(
      Hmisc::approxExtrap(x=.$Year, y = .$value, xout =  year_range))) %>%
    unnest(Demand) %>%
    select(scenario, Year = x, value = y) %>%
    group_by(scenario) %>%
    mutate(value2 = ifelse(Year < 2034, value, value[which(Year == 2034)])) %>%
    ungroup()

  #Plot as ribbons
  demand_scn_df  <- demand_scn_yr %>% rename(plot_var = value2)
  demand_grid_df <- demand_grid_yr %>% rename(plot_var = value2)

  demand_range <- demand_grid_df %>%
    group_by(Year) %>%
    summarize(ymin = min(plot_var), ymax = max(plot_var))

  filter(demand_scn_df, Year == 2034)
  filter(demand_grid_df, Year == 2034)

  p1 <- ggplot(demand_scn_df, aes(x = Year)) +
    theme_bw() +
    geom_line(aes(color = scenario, y = plot_var), size = 0.8) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), data = demand_range, alpha = 0.2) +
    scale_x_continuous(limits=c(2015, 2070), breaks= c(2015, seq(2020,2070,10))) +
    #scale_y_continuous(limits = c(30, 70), breaks = seq(30,70,10)) +
    geom_vline(xintercept = 2020, linetype = "dashed") +
    geom_vline(xintercept = 2070, linetype = "dashed") +
    labs(x = "", y = "Demand (MCM)")

  #Plot as time-series
  p2 <- ggplot(mapping = aes(x = Year, y = plot_var)) +
    geom_line(aes(color = scenario), data = demand_scn_df, size = 1) +
    geom_line(aes(group = scenario), data = demand_grid_df, size = 1) +
    geom_vline(xintercept = 2020, linetype = "dashed") +
    geom_vline(xintercept = 2070, linetype = "dashed") +
    labs(x = "", y = "Demand (MCM)") #+
  #scale_y_continuous(limits = c(30, 70), breaks = seq(30,70,10))


  ggsave(plot = p, filename = "demand_proj.png", width = 7, height =4)

  df <- demand_grid_df %>% select(level = scenario, Year, value) %>%
    mutate(level = factor(level, levels = c("dem1", "dem2", "dem3"),
                          labels = c(1,2,3)))

  #write_csv(x = df, 'demand_levels_2034.csv')















