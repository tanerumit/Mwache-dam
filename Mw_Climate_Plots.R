
#==============================================================================#
# Mombasa Project: Simulation and Analysis                                     #
# M.Umit Taner, Jan 2015                                                       #
# Hydrologic model is based on abcd                                            #
# Water Resources Model application is based on WEAP                           #
#==============================================================================#


  source("Mw_Input_pars.R")
  library(corrplot)
  library(scales)

# MET_DATA_ANALYSIS ------------------------------------------------------------

  #Read-in data
  grid_prec_long <- read.csv("./inputs/grid_prec_long_mm.csv") %>%
    as_data_frame() %>% select(-X) %>% mutate(Date = as.Date(Date))
  obs_prec_long  <- read.csv("./inputs/obs_prec_long_mm.csv") %>%
    as_data_frame() %>% select(Date, Source = Station_ID, value) %>%
    mutate(Source = paste0("RGS_", Source), Date = as.Date(Date))
  #Read-in monthly streamflow time-series (cms) & convert to mm
  qobs_date <- seq.Date(as.Date("1976-06-1"), as.Date("1990-05-1"), by="month")
  q_obs_dat <- read.csv("inputs/Obs_monthly_flow_cms.csv", skip=1)
  q_obs_df  <- q_obs_dat %>%  tbl_df() %>%
    select(-Year, -Mon) %>%
    mutate(
      Date = qobs_date,
      Flow = Flow_obs * 86400 * days_in_month(qobs_date) / (catchment_area * 10^3))


  #Data in wide-format
  grid_prec_wide <- grid_prec_long %>% spread(key = Source, value = value)
  obs_prec_wide  <- obs_prec_long %>% spread(key = Source, value = value)

  #Correlation matrix of met-stations
  cor_matrix <- cor(obs_prec_wide[,-1], use="pairwise.complete.obs")
  corrplot(corr = cor_matrix, type = "upper",   addCoef.col="black",
    cl.lim = c(0,1), cl.length = 5,
    method = "shade", tl.col="black", tl.srt = 45, diag = TRUE)

  #Scatter plot of data by pairs
  dat_tbl <-  obs_prec_wide
  dat_x <- "RGS_9339002"
  dat_y <- "RGS_9339006"

  p <- GG_SCATTER(dat_x, dat_y, save.plot = FALSE)

  #Diagnostic plots (observed precip)
  plots <- pl_TimeSeries(date = as.Date(obs_prec_wide$Date),
                         data = grid_prec_wide[,-1], unit = "mm",
                         agg_type = "sum", titles = FALSE, Guides = TRUE)
  plots$Monthly
  plots$Annual
  plots$boxplot
  plots$Seasonal
  plots$fdc

  #Ground-data versus gridded-data
  df <- obs_prec_wide %>% select(Date, RGS_9339048) %>%
    left_join(grid_prec_wide, by = "Date") %>%
    left_join(q_obs_df, by ="Date")

  #Ground-data versus gridded-data
  GG_SCATTER(x=df[["RGS_9339048"]], y=df[["GPCC"]], save.plot = TRUE)
  GG_SCATTER(x=df[["RGS_9339048"]], y=df[["CRU"]], save.plot = TRUE)
  GG_SCATTER(x=df[["RGS_9339048"]], y=df[["Princeton"]], save.plot = TRUE)


  GG_SCATTER(x=df[["Princeton"]], y=df[["GPCC"]], save.plot = FALSE)
  GG_SCATTER(x=df[["Princeton"]], y=df[["CRU"]], save.plot = FALSE)
  GG_SCATTER(x=df[["GPCC"]], y=df[["CRU"]], save.plot = FALSE)


  for (i in 1:3) {

    data <- grid_prec_wide
    x_name <- "Observed streamflow"
    y_name <- colnames(data)[i+1]

    df <- data[,c(1,i+1)] %>%
      left_join(q_obs_df, by ="Date") %>%
      select(x = Flow_obs, y = get(y_name))

    r2 = paste("r^2 == ", format(summary(lm(df$y ~ df$x))$r.squared, digits = 3))

    p <- ggplot(df, aes(x,y)) +
      theme_bw(base_size = 6) +
      geom_point(shape=1, color = "blue") +
      annotate("text", x=0, y=600,
               label = r2, parse=T, hjust = 0) +
      labs(x = x_name, y = y_name) +
      theme(axis.title.x = element_text(size=8),
            axis.title.y = element_text(size=8))

    plot_name <- paste0("scatter_",x_name,"_vs_", y_name,".png")
    ggsave(plot_name, plot = p, height = 5, width = 5, unit = "cm")

    }



################################################################################

# WARM DATA (Climate Realizations) ---------------------------------------------

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

################################################################################

# TIME-SERIES ANALYSIS  --------------------------------------------------------

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
    geom_abline(intercept = mean(model$y), slope = sd(model$y)) + ggtitle('qqplot')

  p <- arrangeGrob(pl_1, arrangeGrob(pl_2, pl_3, ncol=2), ncol=1)

  #ACF ON HISTORIC TIME-SERIES
  lambda <- powerTransform(ANNUAL_PRCP)$lambda
  TRANSFORM_DATA <- bcPower(ANNUAL_PRCP,lambda)
  pl_historic_ACF_org  <- PLOT_ACF(ANNUAL_PRCP,lag=10, name.arg='Historical')[['ACF']]
  pl_historic_ACF_trns <- PLOT_ACF(TRANSFORM_DATA,lag=10, name.arg='Historical')[['ACF']]

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
    gws_df  <- data.frame(gws_period, gws_sig, gws_obs, gws_mean, gws_min, gws_max)

    #GLOBAL WAVELET SPECTRA (ALL Time-series)
    p <- ggplot(gws_df, aes(x=gws_period)) +
      theme_bw(base_size=12) +
      ggtitle(var_title[r]) +
      geom_ribbon(aes(ymin=gws_min, ymax=gws_max), alpha=0.15) +
      geom_line(mapping=aes(y=gws_mean), color='black', size=0.6) +
      geom_line(mapping=aes(y=gws_sig), color="red", linetype="dashed", size=0.6) +
      geom_line(mapping=aes(y=gws_obs), color="blue", size=0.6) +
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
      scale_x_continuous(expand=c(0,0), limits=c(0,100), breaks=seq(0,100,25), label=comma) +
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


  ##############################################################################





